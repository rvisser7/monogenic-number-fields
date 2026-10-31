#!/usr/bin/env python3
"""Go through the LMFDB fields of the given degrees and find generators of O_K.

Fields come from the LMFDB database (`from lmfdb import db`, or lmfdb-lite), or from a
file with one `label|coeffs` line per field, e.g. `3.1.244.1|[-6, 1, 0, 1]`. Degree 1 and
2 fields are left out (see generate_quadratic.py).

Two methods:
  * magma (default for degree <= 4): IndexFormEquation, via magma/monogenic.m. Exact: it
    proves which fields are monogenic and finds the canonical generator, so its rows
    replace rows from the PARI search. Needs `magma` on the PATH (or --magma).
  * pari: brute-force search over small elements (see monogenic.gp), any degree. It can
    only show a field is monogenic, so fields where it finds nothing are not written.

The rows are merged into data/degree<n>.txt every --batch fields, so an interrupted run
keeps what it has done, and runs can be split with --shard or disc ranges and run in
parallel.

Needs cypari2: pip install cypari2, or run with sage -python. To get lmfdb-lite in Sage:
    sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"

Usage:
    sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6
    sage -python scripts/find_generators.py --degrees 3-4 --only-new --shard 0/4   # and 1/4, ...
    sage -python scripts/find_generators.py --degrees 5 --method pari -B 4
    sage -python scripts/find_generators.py --degrees 3 --input fields.txt --print
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import NULL, data_path, find_generator, format_coeffs, merge_data, parse_coeffs, read_data
from progress import Progress, fmt_time

MAGMA_SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "magma", "monogenic.m")


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(d for d in out if d >= 3))


def parse_int(s):
    """An integer, also written as 10^6 or 1e6."""
    if s is None:
        return None
    if "^" in s:
        b, e = s.split("^")
        return int(b) ** int(e)
    return int(float(s)) if "e" in s.lower() else int(s)


def db_rows(degree, args):
    """(total, iterator of (label, coeffs)) for the LMFDB fields of this degree."""
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    query = {"degree": degree}
    disc = {}
    if args.min_disc is not None:
        disc["$gte"] = args.min_disc
    if args.max_disc is not None:
        disc["$lte"] = args.max_disc
    if disc:
        query["disc_abs"] = disc
    total = db.nf_fields.count(query)
    rows = db.nf_fields.search(query, ["label", "coeffs"])
    return total, ((r["label"], [int(c) for c in r["coeffs"]]) for r in rows)


def file_rows(path, degree, args):
    rows = []
    with open(path) as f:
        for line in f:
            cells = line.strip().split("|")
            if len(cells) < 2 or not cells[0][:1].isdigit():
                continue  # blank line or header
            label = cells[0]
            n, _, disc_abs, _ = (int(x) for x in label.split("."))
            if n != degree:
                continue
            if args.min_disc is not None and disc_abs < args.min_disc:
                continue
            if args.max_disc is not None and disc_abs > args.max_disc:
                continue
            rows.append((label, parse_coeffs(cells[1])))
    return len(rows), iter(rows)


def run_magma(fields, magma):
    """Run magma/monogenic.m on [(label, coeffs)]; returns {label: row} for the fields
    it decided."""
    with tempfile.TemporaryDirectory() as tmp:
        infile, outfile = os.path.join(tmp, "fields.txt"), os.path.join(tmp, "rows.txt")
        with open(infile, "w") as f:
            for label, coeffs in fields:
                f.write("%s|%s\n" % (label, format_coeffs(coeffs)))
        try:
            proc = subprocess.run([magma, "-b", "infile:=" + infile, "outfile:=" + outfile, MAGMA_SCRIPT],
                                  capture_output=True, text=True)
        except FileNotFoundError:
            sys.exit("magma not found: give its path with --magma, or use --method pari")
        if proc.returncode != 0 or "error" in proc.stdout.lower():
            sys.exit("magma failed:\n" + proc.stdout + proc.stderr)
        rows = {}
        if os.path.exists(outfile):
            with open(outfile) as f:
                for line in f:
                    if line.strip():
                        label, mono, coeffs = line.strip().split("|")
                        rows[label] = (label, int(mono), None if coeffs == NULL else parse_coeffs(coeffs))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 3-4 or 3,5")
    parser.add_argument("--input", help="read fields from this file instead of the LMFDB")
    parser.add_argument("--min-disc", type=parse_int, help="smallest |disc|, e.g. 10^6")
    parser.add_argument("--max-disc", type=parse_int, help="largest |disc|, e.g. 10^8")
    parser.add_argument("--method", choices=["magma", "pari"],
                        help="default: magma for degree <= 4, pari above")
    parser.add_argument("--magma", default="magma", help="the magma executable")
    parser.add_argument("-B", "--bound", type=int, default=3,
                        help="pari: max |coordinate| on the integral basis (default 3)")
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th field, from the i-th")
    parser.add_argument("--only-new", action="store_true",
                        help="skip fields that already have a row in the data files")
    parser.add_argument("--batch", type=int, default=1000,
                        help="merge into the data files after this many fields (default 1000)")
    parser.add_argument("--print", action="store_true",
                        help="print the rows instead of merging them")
    parser.add_argument("--replace", action="store_true",
                        help="pari: overwrite rows already in the data files (magma always does)")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()

    shard, nshards = (int(x) for x in args.shard.split("/"))
    progress = Progress(every=args.progress_every)

    for n in parse_degrees(args.degrees):
        method = args.method or ("magma" if n <= 4 else "pari")
        replace = method == "magma" or args.replace
        t0 = time.time()
        total, rows = file_rows(args.input, n, args) if args.input else db_rows(n, args)
        done_before = set()
        if args.only_new and os.path.exists(data_path(n)):
            done_before = {row[0] for row in read_data(data_path(n))}
        found, merged = Counter(), Counter()
        batch, todo = {}, []
        checked = 0
        progress.start("degree %d, %s" % (n, method))

        def flush():
            if method == "magma" and todo:
                rows = run_magma(todo, args.magma)
                found.update(row[1] for row in rows.values())
                batch.update(rows)
                todo.clear()
            for row in batch.values():
                if args.print:
                    print("%s|%d|%s" % (row[0], row[1], NULL if row[2] is None else format_coeffs(row[2])),
                          flush=True)
            if batch and not args.print:
                merged.update(dict(zip(["added", "replaced", "kept"],
                                       merge_data(data_path(n), batch, replace))))
            batch.clear()

        for i, (label, coeffs) in enumerate(rows):
            if i % nshards != shard or label in done_before:
                continue
            checked += 1
            if method == "magma":
                todo.append((label, coeffs))
            else:
                gen = find_generator(coeffs, args.bound)
                if gen is not None:
                    batch[label] = (label, 1, gen)
                    found[1] += 1
            if checked % args.batch == 0:
                flush()
            progress(i + 1, total, ", %d monogenic" % found[1])
        flush()
        print("degree %d: %d fields checked, %d monogenic, %d not monogenic, %d undecided"
              % (n, checked, found[1], found[0], checked - found[1] - found[0])
              + ("; %d added, %d replaced, %d kept" % (merged["added"], merged["replaced"], merged["kept"])
                 if not args.print else "")
              + ", in " + fmt_time(time.time() - t0), file=sys.stderr if args.print else sys.stdout)


if __name__ == "__main__":
    main()
