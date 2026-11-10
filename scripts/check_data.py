"""Validate every data/degree<n>.txt.

For each row with monogenic = 1, with f the polynomial monogenic_coeffs and alpha its root,
this checks that
  * f is monic and irreducible of degree n;
  * [O_L : Z[alpha]] = 1, where L = Q(alpha) and O_L is computed independently by PARI
    (nfdisc): the index is sqrt(disc(f) / disc(O_L));
  * disc(O_L) is the discriminant read off the label;
  * with --source db or --fields: L is the field with that label, by comparing
    polredabs(f) with the LMFDB polynomial (falling back to nfisisom). The
    discriminant alone does not decide this, e.g. 4.0.1040.1 and 4.0.1040.2. If the
    root of the LMFDB polynomial generates O_K, monogenic_coeffs must be that polynomial.
Rows with monogenic = -1 (not monogenic) or 0 (not computed) must have a null monogenic_coeffs.

With --cross-check B the PARI brute-force search (height <= B) is run on each field as an
independent check of the Magma results: it must find no generator for a field marked -1,
and no generator smaller (see monogenic.gp) than the stored one for a field marked 1.

With --canonical, each generator that is not the root of the polredabs polynomial must be
the canonical one: every element of O_K with T2 at most that of the stored generator is
checked (see canonicalize.py, which replaces the ones that are not). Fields with more than
--max-points such elements are counted, not errors.

With --reprove, PARI decides each field again where it can: cubic fields exactly by
solving the Thue equation (the row must match, generator included), quartic fields exactly
with --method factor (likewise), and for the other rows marked -1 it looks for a common
index divisor or local obstruction. A -1 row
without such a proof is counted but is not an error, since Magma may have proved it some
other way.

    python3 scripts/check_data.py
    sage -python scripts/check_data.py --source db --degrees 3-4 --cross-check 2
    sage -python scripts/check_data.py --source db --degrees 3 --reprove
    python3 scripts/check_data.py --fields fields.txt
    python3 scripts/check_data.py --degrees 5-7 --canonical
"""

import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import (DATA_DIR, canonical_generator, cubic_generator, is_polredabs, factor_generator, find_generator, obstruction, format_data, generator_key, is_smaller, label_disc,
                       format_coeffs, label_key, locked, pari, parse_coeffs, read_data, write_data)
from progress import Progress


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return set(out)


def db_fields(labels, batch=1000):
    """{label: coeffs} from the LMFDB for the given labels."""
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    out = {}
    for i in range(0, len(labels), batch):
        for r in db.nf_fields.search({"label": {"$in": labels[i:i + batch]}}, ["label", "coeffs"]):
            out[r["label"]] = [int(c) for c in r["coeffs"]]
    return out


def check_generator(label, n, coeffs, lmfdb_coeffs):
    """Errors for a row with monogenic = 1."""
    f = pari.Polrev(coeffs)
    if len(coeffs) != n + 1 or coeffs[-1] != 1 or not pari.polisirreducible(f):
        return ["%s: not a monic irreducible polynomial of degree %d" % (label, n)]
    errors = []
    D = pari.nfdisc(f)
    index = pari.sqrtint(pari.poldisc(f) / D)
    if index != 1:
        errors.append("%s: Z[alpha] has index %s in the ring of integers" % (label, index))
    if D != label_disc(label):
        errors.append("%s: monogenic_coeffs defines a field of discriminant %s, expected %s"
                      % (label, D, label_disc(label)))
    if lmfdb_coeffs is None:
        return errors
    g = pari.Polrev(lmfdb_coeffs)
    if pari.polredabs(f) != g and not pari.nfisisom(f, g):
        errors.append("%s: monogenic_coeffs defines a different field from the LMFDB polynomial" % label)
    elif pari.poldisc(g) == D and coeffs != lmfdb_coeffs:
        errors.append("%s: the root of the LMFDB polynomial generates O_K, so monogenic_coeffs should be %s"
                      % (label, format_coeffs(lmfdb_coeffs)))
    return errors


def cross_check(label, mono, coeffs, field_coeffs, bound):
    """Errors from comparing a row with the PARI brute-force search."""
    found = find_generator(field_coeffs, bound)
    if found is None:
        return []
    if mono == -1:
        return ["%s: marked not monogenic, but PARI finds the generator %s" % (label, format_coeffs(found))]
    if mono == 1 and is_smaller(generator_key(found), generator_key(coeffs)):
        return ["%s: PARI finds the generator %s, smaller than the stored one" % (label, format_coeffs(found))]
    return []


def reprove(label, n, mono, coeffs, field_coeffs):
    """(errors, proved): errors from deciding the field again with PARI, and whether a -1
    row of degree > 3 has a PARI proof."""
    decided, how = False, None
    if n == 3:
        decided, gen, how = True, cubic_generator(field_coeffs), "the Thue equation"
    elif n == 4:
        decided, gen = factor_generator(field_coeffs)
        how = "the factored index form"
    if decided:
        if gen is None and mono == 1:
            return ["%s: marked monogenic, but %s has no solutions" % (label, how)], False
        if gen is not None and mono == -1:
            return ["%s: marked not monogenic, but %s gives %s" % (label, how, format_coeffs(gen))], False
        if gen is not None and gen != coeffs:
            return ["%s: %s gives the generator %s, not the stored one" % (label, how, format_coeffs(gen))], False
        return [], True
    if mono == -1:
        return [], obstruction(field_coeffs) is not None
    return [], False


def check_file(path, fields, cross_bound, progress, do_reprove=False, canonical_points=None):
    errors = []
    n = int(os.path.basename(path)[len("degree"):-len(".txt")])
    rows = read_data(path)
    with open(path) as fh:
        if fh.read() != format_data(rows):
            errors.append("not in the standard format (header, sorted, no spaces in monogenic_coeffs); "
                          "rewrite it with --fix")
    labels = [row[0] for row in rows]
    if labels != sorted(labels, key=label_key) or len(set(labels)) != len(labels):
        errors.append("rows are not sorted or not unique")
    skipped = unproved = not_canonical = 0
    progress.start(os.path.basename(path))
    for i, (label, mono, coeffs) in enumerate(rows):
        if int(label.split(".")[0]) != n:
            errors.append("%s: wrong degree file" % label)
        elif mono in (-1, 0):
            if coeffs is not None:
                errors.append("%s: monogenic = %d but monogenic_coeffs given" % (label, mono))
        elif mono != 1 or coeffs is None:
            errors.append("%s: monogenic must be 1, -1 or 0, with monogenic_coeffs iff 1" % label)
        else:
            errors += check_generator(label, n, coeffs, fields.get(label))
            if canonical_points and not is_polredabs(coeffs):
                gen, proved = canonical_generator(coeffs, canonical_points)
                if not proved:
                    not_canonical += 1
                elif gen != coeffs:
                    errors.append("%s: not the canonical generator: %s comes first (fix with "
                                  "canonicalize.py)" % (label, format_coeffs(gen)))
        if (cross_bound or do_reprove) and mono == -1 and label not in fields:
            skipped += 1
        if do_reprove and label in fields and (n in (3, 4) and mono != 0 or mono == -1):
            errs, proved = reprove(label, n, mono, coeffs, fields[label])
            errors += errs
            unproved += mono == -1 and not proved and not errs
        if cross_bound and mono != 0:
            # the polredabs polynomial if we have it; otherwise monogenic_coeffs defines the field
            field_coeffs = fields.get(label, coeffs)
            if field_coeffs is not None:
                errors += cross_check(label, mono, coeffs, field_coeffs, cross_bound)
        progress(i + 1, len(rows), ", %d errors" % len(errors))
    return len(rows), errors, skipped, unproved, not_canonical


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", help="only these degrees, e.g. 3-4 or 3,5")
    parser.add_argument("--source", choices=["db"],
                        help="get the LMFDB polynomials from the database (lmfdb or lmfdb-lite)")
    parser.add_argument("--fields", help='or from a file with lines "label|[coeffs]"')
    parser.add_argument("--cross-check", type=int, metavar="B",
                        help="also run the PARI search up to height B on every field")
    parser.add_argument("--reprove", action="store_true",
                        help="decide the fields again with PARI where it can (see above)")
    parser.add_argument("--canonical", action="store_true",
                        help="also prove that each generator is the canonical one (see above)")
    parser.add_argument("--max-points", type=int, default=10**5,
                        help="--canonical: the most elements to check for one field (default 10^5)")
    parser.add_argument("--fix", action="store_true",
                        help="first rewrite each file in the standard format")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()
    degrees = parse_degrees(args.degrees) if args.degrees else None
    fields = {}
    if args.fields:
        with open(args.fields) as fh:
            for line in fh:
                if line.strip():
                    label, coeffs = line.strip().split("|")
                    fields[label] = parse_coeffs(coeffs)
    progress = Progress(every=args.progress_every)

    ok = True
    for path in sorted(glob.glob(os.path.join(DATA_DIR, "degree*.txt")), key=lambda p: (len(p), p)):
        if degrees is not None and int(os.path.basename(path)[len("degree"):-len(".txt")]) not in degrees:
            continue
        if args.fix:
            with locked(path):
                write_data(path, read_data(path))
        if args.source == "db":
            fields.update(db_fields([row[0] for row in read_data(path)]))
        count, errors, skipped, unproved, not_canonical = check_file(
            path, fields, args.cross_check, progress, args.reprove,
            args.max_points if args.canonical else None)
        print("%s: %d rows, %d errors" % (os.path.basename(path), count, len(errors))
              + ("; %d rows with monogenic = -1 not cross-checked, as that needs the LMFDB "
                 "polynomial (--source db or --fields)" % skipped if skipped else "")
              + ("; %d rows with monogenic = -1 have no PARI proof (not an error)" % unproved
                 if unproved else "")
              + ("; %d generators not proved canonical (more than --max-points elements)" % not_canonical
                 if not_canonical else ""))
        for e in errors:
            print("  " + e)
        ok = ok and not errors
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
