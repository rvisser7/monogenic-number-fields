#!/usr/bin/env python3
"""Go through the LMFDB fields of the given degrees and find generators of O_K.

Fields come from the LMFDB database (`from lmfdb import db`, or lmfdb-lite), or from a
file with one `label|coeffs` line per field, e.g. `3.1.244.1|[-6,1,0,1]`. Degree 1 and
2 fields are left out (see generate_quadratic.py), except by --method polredabs.

Methods:
  * polredabs: the quick bulk pass. Keeps exactly the fields where the root of the
    LMFDB (polredabs) polynomial f already generates O_K, i.e. disc(f) = disc(K), and
    writes them with zk_coeffs = f, which is the canonical generator (see README). This
    costs one polynomial discriminant per field, so a whole degree takes about as long as
    reading it from the LMFDB. With --lmfdb-index only the fields whose LMFDB `index`
    column is 1 are read (faster still; each is still checked). Other fields are left
    for the methods below. By rule 1 these replace any other stored generator.
  * magma (default for degree <= 4): IndexFormEquation, via magma/monogenic.m. Exact: it
    proves which fields are monogenic and finds the canonical generator, so its
    generators replace larger ones from the PARI search (see below). Needs `magma` on
    the PATH (or --magma).
  * thue (degree 3 only): PARI solves the index form equation, a cubic Thue equation,
    completely. Exact like magma, and gives the same rows, without needing Magma.
  * factor (meant for degree 4): factors the index form over Z. For a quartic field with
    an imaginary quadratic subfield one factor is a positive definite binary quadratic
    form, which leaves finitely many cases, each a polynomial equation in one variable
    (see monogenic.gp). Exact like thue for those fields, both ways, without Magma. For a
    real quadratic subfield it uses the cubic resolvent to prove most non-monogenic fields
    (it never finds generators there). Other fields are skipped. E.g. for the quartic
    fields not decided yet:
        --degrees 4 --method factor --only-new
  * pari: brute-force search over small elements (see monogenic.gp), any degree. When it
    finds nothing, the tests in `local` are tried; fields where both fail are not written.
  * local: only the tests that prove a field is NOT monogenic (common index divisors,
    and the index form missing +-1 modulo small prime powers; see monogenic.gp). Writes
    rows with monogenic = -1 and nothing else, e.g. for the fields the search left over:
        --method local --only-new

Conflicts: when a new row says 1 and the stored row -1, or the other way round, the
generator is checked independently (PARI computes the ring of integers of the field it
defines and checks [O_K : Z[alpha]] = 1 and the discriminant), whatever --replace says.
If it passes, the 1 row wins, since a generator is a certificate and a -1 from Magma (or
anywhere else) can be a bug; if not, the -1 row wins. Each conflict is printed as
"CONFLICT ..." and appended to conflicts.txt. To look again for generators of the fields
marked -1, e.g. the ones Magma decided:
    --method pari --recheck-negative -B 4
and --recheck-generators reruns a method on the fields whose stored generator is not the
LMFDB polynomial, e.g. after a change to the canonical order.

LMFDB monogenic fields: --lmfdb-monogenic reads only the fields that the LMFDB marks
monogenic (nf_fields.monogenic = 1) and that have no generator here yet (no row, or a
row marked 0 or -1), e.g. to run the search on just those:
    --degrees 6 --method pari --lmfdb-monogenic -B 4
A field marked -1 here is a conflict with the LMFDB, settled as above if the search
finds a generator; add --only-new to leave those out. compare_lmfdb.py lists the same
fields (lmfdb_only_yes.txt, and the conflicts).

Generators: whatever the method, a new generator replaces a different stored one only if
it comes first in the canonical order (the polredabs root, then smallest T2 norm, then
smallest minimal polynomial; see README) and passes the same check. Each replacement is
printed as "BETTER ...". --replace makes any new generator that passes the check
replace the stored one.

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
    sage -python scripts/find_generators.py --degrees 3 --method thue --max-disc 10^5
    sage -python scripts/find_generators.py --degrees 5-8 --method local --only-new --show-proofs
    sage -python scripts/find_generators.py --degrees 4 --method factor --only-new
    sage -python scripts/find_generators.py --degrees 2-47 --method polredabs
    sage -python scripts/find_generators.py --degrees 6 --method polredabs --lmfdb-index
    sage -python scripts/find_generators.py --degrees 3-4 --method pari --recheck-negative -B 4
    sage -python scripts/find_generators.py --degrees 3 --method thue --recheck-generators
    sage -python scripts/find_generators.py --degrees 6 --method pari --lmfdb-monogenic -B 4
"""

import argparse
import os
import subprocess
import sys
import tempfile
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import (NULL, cubic_generator, data_path, factor_generator, find_generator, format_coeffs,
                       merge_data, obstruction, parse_coeffs, read_data, root_generates)
from progress import Progress, fmt_time

MAGMA_SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "magma", "monogenic.m")


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(d for d in out if d >= 2))


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
    if args.lmfdb_index:
        query["index"] = 1
    if args.lmfdb_monogenic:
        query["monogenic"] = 1
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
    parser.add_argument("--method", choices=["magma", "thue", "factor", "pari", "local", "polredabs"],
                        help="default: magma for degree <= 4, pari above")
    parser.add_argument("--lmfdb-index", action="store_true",
                        help="polredabs: only read the fields whose LMFDB index column is 1")
    parser.add_argument("--magma", default="magma", help="the magma executable")
    parser.add_argument("-B", "--bound", type=int, default=3,
                        help="pari: max |coordinate| on the integral basis (default 3)")
    parser.add_argument("--maxp", type=int, default=50,
                        help="pari, local: try local obstructions at primes up to this (default 50)")
    parser.add_argument("--max-points", type=parse_int, default=10**5,
                        help="pari, local: skip a modulus needing more evaluations (default 10^5)")
    parser.add_argument("--show-proofs", action="store_true",
                        help="print the proof found for each field marked -1 by pari, local or factor")
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th field, from the i-th")
    parser.add_argument("--only-new", action="store_true",
                        help="skip fields that already have a row in the data files")
    parser.add_argument("--lmfdb-monogenic", action="store_true",
                        help="only the fields the LMFDB marks monogenic that have no generator "
                             "in the data files (no row, or marked 0 or -1)")
    parser.add_argument("--recheck-negative", action="store_true",
                        help="only the fields marked -1 in the data files, e.g. to search "
                             "for generators of the fields Magma says are not monogenic")
    parser.add_argument("--recheck-generators", action="store_true",
                        help="only the fields marked 1 whose generator is not the LMFDB "
                             "polynomial, to replace it by a smaller one if there is one")
    parser.add_argument("--batch", type=int,
                        help="merge into the data files after this many fields "
                             "(default 1000, and 10^6 for polredabs)")
    parser.add_argument("--print", action="store_true",
                        help="print the rows instead of merging them")
    parser.add_argument("--replace", action="store_true",
                        help="let any verified new generator replace the stored one, not "
                             "just one that comes first in the canonical order (1 against -1 "
                             "is decided by checking the generator, with or without this)")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()

    shard, nshards = (int(x) for x in args.shard.split("/"))
    progress = Progress(every=args.progress_every)

    if args.lmfdb_index and args.method != "polredabs":
        sys.exit("--lmfdb-index is for --method polredabs")
    recheck = args.recheck_negative or args.recheck_generators
    if recheck and args.only_new:
        sys.exit("--recheck-negative and --recheck-generators exclude --only-new")
    if args.lmfdb_monogenic and args.input:
        sys.exit("--lmfdb-monogenic reads the LMFDB; for a file, use the lists from "
                 "compare_lmfdb.py --write (lmfdb_only_yes.txt)")
    if args.lmfdb_monogenic and args.recheck_generators:
        sys.exit("--lmfdb-monogenic skips the fields with a generator, so it excludes --recheck-generators")

    for n in parse_degrees(args.degrees):
        method = args.method or ("magma" if n <= 4 else "pari")
        if n == 2 and method != "polredabs":
            continue  # see generate_quadratic.py
        if method == "thue" and n != 3:
            sys.exit("--method thue is for degree 3 only")
        batch_size = args.batch or (10**6 if method == "polredabs" else 1000)
        # merge_data settles a 1 against a -1 by verifying the generator, and two
        # generators by the canonical order (so exact methods beat the search).
        replace = args.replace
        t0 = time.time()
        total, rows = file_rows(args.input, n, args) if args.input else db_rows(n, args)
        old = {}
        if (args.only_new or recheck or args.lmfdb_monogenic or method == "polredabs") \
                and os.path.exists(data_path(n)):
            old = {row[0]: row for row in read_data(data_path(n))}
        done_before = old if args.only_new else {}
        found, merged, skipped = Counter(), Counter(), Counter()
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
                merged.update(dict(zip(["added", "replaced", "kept", "conflicts"],
                                       merge_data(data_path(n), batch, replace))))
            batch.clear()

        for i, (label, coeffs) in enumerate(rows):
            if i % nshards != shard or label in done_before:
                continue
            if recheck:
                prev = old.get(label)
                if not (prev and (args.recheck_negative and prev[1] == -1 or
                                  args.recheck_generators and prev[1] == 1 and prev[2] != coeffs)):
                    continue
            if args.lmfdb_monogenic and old.get(label, (label, 0))[1] == 1:
                continue  # we already have a generator
            checked += 1
            if method == "magma":
                todo.append((label, coeffs))
            elif method == "polredabs" and root_generates(coeffs, label):
                found[1] += 1
                row = (label, 1, coeffs)
                prev = old.get(label)
                if prev == row:
                    skipped["already there"] += 1
                else:
                    # merge_data settles a -1 (a conflict) and another generator (rule 1)
                    batch[label] = row
            elif method == "thue":
                gen = cubic_generator(coeffs)
                batch[label] = (label, 1, gen) if gen is not None else (label, -1, None)
                found[1 if gen is not None else -1] += 1
            elif method == "factor":
                decided, gen = factor_generator(coeffs)
                if decided:
                    batch[label] = (label, 1, gen) if gen is not None else (label, -1, None)
                    found[1 if gen is not None else -1] += 1
                    if args.show_proofs and gen is None:
                        print("%s: not monogenic, factored index form" % label, file=sys.stderr, flush=True)
            elif method in ("pari", "local"):
                gen = find_generator(coeffs, args.bound) if method == "pari" else None
                proof = obstruction(coeffs, args.maxp, args.max_points) if gen is None else None
                if gen is not None:
                    batch[label] = (label, 1, gen)
                    found[1] += 1
                elif proof is not None:
                    batch[label] = (label, -1, None)
                    found[-1] += 1
                    if args.show_proofs:
                        print("%s: not monogenic, %s" % (label, proof), file=sys.stderr, flush=True)
            if checked % batch_size == 0:
                flush()
            progress(i + 1, total, ", %d monogenic, %d not" % (found[1], found[-1]))
        flush()
        print("degree %d: %d fields checked, %d monogenic, %d not monogenic, %d undecided"
              % (n, checked, found[1], found[-1], checked - found[1] - found[-1])
              + ("; %d added, %d replaced, %d kept" % (merged["added"], merged["replaced"], merged["kept"])
                 + ("; %d CONFLICTS with the stored rows (see conflicts.txt)" % merged["conflicts"]
                    if merged["conflicts"] else "")
                 if not args.print else "")
              + "".join("; %d %s" % (k, why) for why, k in sorted(skipped.items()))
              + ", in " + fmt_time(time.time() - t0), file=sys.stderr if args.print else sys.stdout)


if __name__ == "__main__":
    main()
