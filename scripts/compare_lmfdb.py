#!/usr/bin/env python3
"""Compare the monogenic column of data/degree<n>.txt with the LMFDB's nf_fields.monogenic.

Both use the same values: 1 = monogenic, -1 = not monogenic, 0 = not computed. The LMFDB
column can also be null, which counts as 0 here (and is counted separately in the
summary). A field missing from our data files counts as 0.

For each degree this prints a summary line, and then lists
  * conflicts: ours 1 and LMFDB -1, or ours -1 and LMFDB 1. One of the two is wrong. For
    ours 1 the stored monogenic_coeffs is a generator that check_data.py verifies, so it is a
    certificate that the LMFDB value is wrong;
  * fields where only we know the answer (ours 1 or -1, LMFDB 0 or null), which the
    LMFDB could take from us;
  * fields where only the LMFDB knows the answer (LMFDB 1 or -1, ours 0), which are
    still to do here;
and any LMFDB values other than 1, -1, 0 and null. Only the first --show of each list are
printed; --write DIR writes the whole lists to files (see below). The exit status is 1 if
there is a conflict or an unexpected value, and 0 otherwise.

The LMFDB values come from the database (`from lmfdb import db`, or lmfdb-lite, as in
find_generators.py): one query per degree for the fields with monogenic != 0, plus
counts. Or they come from a file (--lmfdb FILE) with lines `label|monogenic` or
`label|monogenic|coeffs` (`\\N` for null; header lines are skipped), e.g. an export of
nf_fields; then the fields of the file are taken to be all the LMFDB fields of those
degrees.

With --min-disc and/or --max-disc only the fields with |disc| in that range are compared,
on both sides (the LMFDB query gets a disc_abs condition, and our rows are filtered by the
discriminant in their label), so a degree with millions of fields can be checked in part.

While it runs it prints progress lines on stdout, starting with the time elapsed: what it is
reading or counting, the number of LMFDB rows read so far every --progress-every seconds
(default 10), and a one-line result as each degree finishes. --quiet leaves them out.

Files written by --write DIR, with all degrees together:
  conflicts.txt          label|ours|lmfdb|monogenic_coeffs (ours)
  ours_only.txt          label|monogenic, for the fields only we know: an LMFDB update
  lmfdb_only_yes.txt     label|[coeffs] (the LMFDB polynomial), LMFDB 1 and ours 0: input
                         for find_generators.py --input
  lmfdb_only_no.txt      label|[coeffs], LMFDB -1 and ours 0: input for
                         find_generators.py --input --method local, or for Magma
  unexpected.txt         label|lmfdb, LMFDB values other than 1, -1, 0, null

Usage:
    sage -python scripts/compare_lmfdb.py
    sage -python scripts/compare_lmfdb.py --degrees 3-4 --show 50
    sage -python scripts/compare_lmfdb.py --degrees 6 --max-disc 10^4
    sage -python scripts/compare_lmfdb.py --write compare/
    python3 scripts/compare_lmfdb.py --lmfdb nf_fields_monogenic.txt
"""

import argparse
import glob
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import DATA_DIR, NULL, data_path, format_coeffs, label_key, parse_coeffs, read_data
from progress import fmt_time

UNKNOWN = (0, None)


class Log:
    """Progress lines on stdout, each with the time since the start, so a long run
    shows it is busy. say() prints at once; tick() at most every `every` seconds."""

    def __init__(self, every=10, quiet=False):
        self.t0 = self.last = time.time()
        self.every = every
        self.quiet = quiet

    def say(self, msg):
        self.last = time.time()
        if not self.quiet:
            print("[%s] %s" % (fmt_time(self.last - self.t0), msg), flush=True)

    def tick(self, msg):
        if time.time() - self.last >= self.every:
            self.say(msg)


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(out))


def parse_int(s):
    """An integer, also written as 10^6 or 1e6."""
    if "^" in s:
        b, e = s.split("^")
        return int(b) ** int(e)
    return int(float(s)) if "e" in s.lower() else int(s)


class DiscRange:
    """min_disc <= |disc| <= max_disc, either bound None for no bound."""

    def __init__(self, min_disc=None, max_disc=None):
        self.min, self.max = min_disc, max_disc

    def query(self):
        """The disc_abs condition for an LMFDB query, as a dict to add to the query."""
        cond = {}
        if self.min is not None:
            cond["$gte"] = self.min
        if self.max is not None:
            cond["$lte"] = self.max
        return {"disc_abs": cond} if cond else {}

    def __contains__(self, label):
        D = int(label.split(".")[2])
        return (self.min is None or D >= self.min) and (self.max is None or D <= self.max)

    def __str__(self):
        if self.min is None and self.max is None:
            return ""
        return " with %s|disc|%s" % ("%d <= " % self.min if self.min is not None else "",
                                     " <= %d" % self.max if self.max is not None else "")


def data_degrees():
    return sorted(int(os.path.basename(p)[len("degree"):-len(".txt")])
                  for p in glob.glob(os.path.join(DATA_DIR, "degree*.txt")))


class DBSource:
    """The LMFDB values from the database."""

    def __init__(self, with_coeffs, disc=None):
        self.disc = disc or DiscRange()
        try:
            from lmfdb import db
        except ImportError:
            from lmf import db
        self.table = db.nf_fields
        self.cols = ["label", "monogenic"] + (["coeffs"] if with_coeffs else [])

    def degree(self, n, log):
        """(known, counts): known is {label: (monogenic, coeffs or None)} for the fields
        with monogenic not 0 or null, counts is a Counter with keys 1, -1, 0, None (the
        number of fields with that value), 'other' and 'total'."""
        counts = Counter()
        base = dict(degree=n, **self.disc.query())
        log.say("degree %d: counting the LMFDB fields%s" % (n, self.disc))
        counts["total"] = self.table.count(base)
        counts[0] = self.table.count(dict(base, monogenic=0))
        log.say("degree %d: %d LMFDB fields, %d with monogenic = 0; reading the rest"
                % (n, counts["total"], counts[0]))
        known = {}
        for r in self.table.search(dict(base, monogenic={"$ne": 0}), self.cols):
            known[r["label"]] = (r["monogenic"], [int(c) for c in r["coeffs"]] if "coeffs" in r else None)
            log.tick("degree %d: %d LMFDB fields with monogenic != 0 read so far" % (n, len(known)))
        log.say("degree %d: %d LMFDB fields with monogenic != 0 read" % (n, len(known)))
        for mono, _ in known.values():
            counts[mono if mono in (1, -1) else "other"] += 1
        counts[None] = counts["total"] - counts[0] - counts[1] - counts[-1] - counts["other"]
        return known, counts


class FileSource:
    """The LMFDB values from a file of lines label|monogenic[|coeffs]."""

    def __init__(self, path, disc=None):
        disc = disc or DiscRange()
        self.fields = {}
        with open(path) as fh:
            for line in fh:
                cells = line.strip().split("|")
                if len(cells) < 2 or not cells[0][:1].isdigit():
                    continue  # blank or header line
                label, mono = cells[0], cells[1]
                if label not in disc:
                    continue
                coeffs = parse_coeffs(cells[2]) if len(cells) > 2 and cells[2] != NULL else None
                self.fields.setdefault(int(label.split(".")[0]), {})[label] = (
                    None if mono == NULL else int(mono), coeffs)

    def degrees(self):
        return sorted(self.fields)

    def degree(self, n, log):
        fields = self.fields.get(n, {})
        known = {label: v for label, v in fields.items() if v[0] not in UNKNOWN}
        counts = Counter()
        counts["total"] = len(fields)
        for mono, _ in fields.values():
            counts[mono if mono in (1, -1, 0, None) else "other"] += 1
        return known, counts


def compare(ours, known):
    """Compare {label: (monogenic, monogenic_coeffs)} (ours) with {label: (monogenic, coeffs)}
    (the LMFDB fields with monogenic not 0 or null). Returns a dict of lists of labels:
    agree_yes, agree_no, conflict, ours_only_yes, ours_only_no, lmfdb_only_yes,
    lmfdb_only_no, unexpected."""
    out = {k: [] for k in ["agree_yes", "agree_no", "conflict", "ours_only_yes", "ours_only_no",
                           "lmfdb_only_yes", "lmfdb_only_no", "unexpected"]}
    for label, (lm, _) in known.items():
        mine = ours.get(label, (0, None))[0]
        if lm not in (1, -1):
            out["unexpected"].append(label)
        elif mine == 0:
            out["lmfdb_only_yes" if lm == 1 else "lmfdb_only_no"].append(label)
        elif mine == lm:
            out["agree_yes" if lm == 1 else "agree_no"].append(label)
        else:
            out["conflict"].append(label)
    for label, (mine, _) in ours.items():
        if mine != 0 and label not in known:
            out["ours_only_yes" if mine == 1 else "ours_only_no"].append(label)
    for labels in out.values():
        labels.sort(key=label_key)
    return out


SUMMARY = [("ours 1", None), ("ours -1", None), ("LMFDB 1", None), ("LMFDB -1", None),
           ("LMFDB 0", None), ("LMFDB null", None), ("LMFDB fields", None), ("agree 1", "agree_yes"),
           ("agree -1", "agree_no"), ("conflicts", "conflict"), ("only ours 1", "ours_only_yes"),
           ("only ours -1", "ours_only_no"), ("only LMFDB 1", "lmfdb_only_yes"),
           ("only LMFDB -1", "lmfdb_only_no")]


def summary_row(n, ours, counts, result):
    mine = Counter(m for m, _ in ours.values())
    values = [mine[1], mine[-1], counts[1], counts[-1], counts[0], counts[None], counts["total"]]
    values += [len(result[key]) for _, key in SUMMARY[len(values):]]
    return [str(n)] + [str(v) for v in values]


def print_table(rows):
    header = ["degree"] + [name for name, _ in SUMMARY]
    widths = [max(len(r[i]) for r in [header] + rows) for i in range(len(header))]
    for r in [header] + rows:
        print("  ".join(c.rjust(w) for c, w in zip(r, widths)))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", help="only these degrees, e.g. 3-4 or 3,5 (default: those "
                                          "with a data file, and those in --lmfdb)")
    parser.add_argument("--lmfdb", metavar="FILE",
                        help="read the LMFDB values from this file instead of the database")
    parser.add_argument("--show", type=int, default=10, metavar="N",
                        help="print at most N fields of each list per degree (default 10; "
                             "conflicts and unexpected values are always all printed)")
    parser.add_argument("--write", metavar="DIR", help="write the full lists to files in DIR")
    parser.add_argument("--progress-every", type=float, default=10, metavar="SECONDS",
                        help="while reading from the LMFDB, print progress this often (default 10)")
    parser.add_argument("-q", "--quiet", action="store_true", help="no progress lines")
    parser.add_argument("--min-disc", type=parse_int, help="only fields with |disc| >= this, e.g. 10^4")
    parser.add_argument("--max-disc", type=parse_int, help="only fields with |disc| <= this, e.g. 10^4")
    args = parser.parse_args()
    disc = DiscRange(args.min_disc, args.max_disc)
    log = Log(args.progress_every, args.quiet)

    log.say("reading %s" % args.lmfdb if args.lmfdb else "connecting to the LMFDB")
    source = FileSource(args.lmfdb, disc) if args.lmfdb else DBSource(bool(args.write), disc)
    if args.degrees:
        degrees = parse_degrees(args.degrees)
    else:
        degrees = sorted(set(data_degrees()) | set(source.degrees() if args.lmfdb else []))
    files = {}
    if args.write:
        os.makedirs(args.write, exist_ok=True)
        for name in ["conflicts", "ours_only", "lmfdb_only_yes", "lmfdb_only_no", "unexpected"]:
            files[name] = open(os.path.join(args.write, name + ".txt"), "w")

    table, total, details = [], Counter(), []
    bad = Counter()
    for n in degrees:
        ours = {}
        if os.path.exists(data_path(n)):
            log.say("degree %d: reading %s" % (n, os.path.basename(data_path(n))))
            ours = {label: (mono, coeffs) for label, mono, coeffs in read_data(data_path(n)) if label in disc}
        known, counts = source.degree(n, log)
        if args.lmfdb and not counts["total"]:
            details.append("degree %d: no fields%s in %s, skipped" % (n, disc, args.lmfdb))
            continue
        result = compare(ours, known)
        log.say("degree %d: %d rows here, %d LMFDB fields; %d agree, %d conflicts, %d only here, "
                "%d only in the LMFDB" % (
                    n, len(ours), counts["total"], len(result["agree_yes"]) + len(result["agree_no"]),
                    len(result["conflict"]), len(result["ours_only_yes"]) + len(result["ours_only_no"]),
                    len(result["lmfdb_only_yes"]) + len(result["lmfdb_only_no"])))
        bad.update(conflicts=len(result["conflict"]), unexpected=len(result["unexpected"]))
        table.append(summary_row(n, ours, counts, result))
        total.update({i: int(v) for i, v in enumerate(table[-1][1:])})

        lines = []
        for label in result["conflict"]:
            mine, zk = ours[label]
            lines.append("  CONFLICT %s: ours %d, LMFDB %d%s" % (
                label, mine, known[label][0], ", generator " + format_coeffs(zk) if zk else ""))
        for label in result["unexpected"]:
            lines.append("  UNEXPECTED %s: LMFDB monogenic = %r" % (label, known[label][0]))
        for key, what in [("ours_only_yes", "only we know are monogenic"),
                          ("ours_only_no", "only we know are not monogenic"),
                          ("lmfdb_only_yes", "only the LMFDB knows are monogenic"),
                          ("lmfdb_only_no", "only the LMFDB knows are not monogenic")]:
            labels = result[key]
            if labels and args.show:
                more = ", ..." if len(labels) > args.show else ""
                lines.append("  %d fields %s: %s%s" % (len(labels), what, ", ".join(labels[:args.show]), more))
        if lines:
            details.append("degree %d:" % n)
            details += lines

        if files:
            for label in result["conflict"]:
                mine, zk = ours[label]
                files["conflicts"].write("%s|%d|%d|%s\n" % (label, mine, known[label][0],
                                                            format_coeffs(zk) if zk else NULL))
            for label in result["ours_only_yes"] + result["ours_only_no"]:
                files["ours_only"].write("%s|%d\n" % (label, ours[label][0]))
            for key in ["lmfdb_only_yes", "lmfdb_only_no"]:
                for label in result[key]:
                    coeffs = known[label][1]
                    files[key].write("%s|%s\n" % (label, format_coeffs(coeffs) if coeffs else NULL))
            for label in result["unexpected"]:
                files["unexpected"].write("%s|%r\n" % (label, known[label][0]))
        sys.stdout.flush()

    log.say("done")
    print()
    if len(table) > 1:
        table.append(["all"] + [str(total[i]) for i in range(len(SUMMARY))])
    print_table(table)
    if details:
        print()
        print("\n".join(details))
    for fh in files.values():
        fh.close()
    print()
    print("%d conflicts, %d unexpected LMFDB values%s" % (bad["conflicts"], bad["unexpected"],
                                                          ", for the fields" + str(disc) if str(disc) else ""))
    sys.exit(1 if bad["conflicts"] or bad["unexpected"] else 0)


if __name__ == "__main__":
    main()
