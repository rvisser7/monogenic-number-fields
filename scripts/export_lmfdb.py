"""Export the monogenic fields to one file for the LMFDB: every row with monogenic = 1 from
data/degree<n>.txt, with only the label and monogenic_coeffs columns,

    label|monogenic_coeffs
    text|numeric[]

    1.1.1.1|[0,1]
    2.0.3.1|[1,-1,1]
    ...

sorted by degree and then like the LMFDB (|disc|, signature, index). Rows with
monogenic = -1 or 0 are left out. Needs only Python (no PARI).

    python3 scripts/export_lmfdb.py --output monogenic_coeffs.txt
    python3 scripts/export_lmfdb.py --degrees 2-4 --output monogenic_coeffs_2-4.txt
    python3 scripts/export_lmfdb.py | gzip > monogenic_coeffs.txt.gz
"""

import argparse
import glob
import os
import re
import sys

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
HEADERS = ("label|monogenic|monogenic_coeffs\ntext|smallint|numeric[]\n\n",
           "label|monogenic|zk_coeffs\ntext|smallint|numeric[]\n\n")
EXPORT_HEADER = "label|monogenic_coeffs\ntext|numeric[]\n\n"
NULL = r"\N"


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return set(out)


def label_key(label):
    """Sort by degree, |disc|, signature, then index, as the LMFDB does."""
    n, r, D, i = (int(x) for x in label.split("."))
    return (n, D, r, i)


def monogenic_rows(path):
    """(label, coeffs) for the rows of a data file with monogenic = 1, coeffs as stored
    (e.g. "[1,-1,1]")."""
    with open(path) as fh:
        text = fh.read()
    if not text.startswith(HEADERS):
        raise ValueError("%s: bad header" % path)
    rows = []
    for lineno, line in enumerate(text.split("\n")[3:], 4):
        if not line:
            continue
        label, mono, coeffs = line.split("|")
        if mono != "1":
            continue
        n = int(label.split(".")[0])
        if coeffs == NULL or not re.fullmatch(r"\[-?\d+(,-?\d+)*\]", coeffs) \
                or coeffs.count(",") != n:
            raise ValueError("%s:%d: %s has monogenic = 1 but monogenic_coeffs %s is not a "
                             "polynomial of degree %d" % (path, lineno, label, coeffs, n))
        rows.append((label, coeffs))
    return rows


def export(paths, out):
    """Write the monogenic rows of the data files `paths` to the file object `out`.
    Returns the number of rows."""
    rows = []
    for path in paths:
        rows += monogenic_rows(path)
    rows.sort(key=lambda row: label_key(row[0]))
    labels = [label for label, _ in rows]
    if len(set(labels)) != len(labels):
        raise ValueError("a label occurs twice")
    out.write(EXPORT_HEADER)
    for label, coeffs in rows:
        out.write("%s|%s\n" % (label, coeffs))
    return len(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", help="only these degrees, e.g. 2-4 or 3,5 (default: all)")
    parser.add_argument("--output", "-o", help="the file to write (default: standard output)")
    parser.add_argument("--data-dir", default=DATA_DIR, help="where degree<n>.txt are")
    args = parser.parse_args()

    paths = glob.glob(os.path.join(args.data_dir, "degree*.txt"))
    paths = [p for p in paths if re.fullmatch(r"degree\d+\.txt", os.path.basename(p))]
    if args.degrees:
        degrees = parse_degrees(args.degrees)
        paths = [p for p in paths if int(os.path.basename(p)[6:-4]) in degrees]
    if not paths:
        sys.exit("no data files found in %s" % args.data_dir)

    if args.output:
        tmp = args.output + ".tmp%d" % os.getpid()
        try:
            with open(tmp, "w") as fh:
                count = export(paths, fh)
            os.replace(tmp, args.output)
        finally:
            if os.path.exists(tmp):
                os.remove(tmp)
        print("%d monogenic fields written to %s" % (count, args.output), file=sys.stderr)
    else:
        count = export(paths, sys.stdout)
        print("%d monogenic fields" % count, file=sys.stderr)


if __name__ == "__main__":
    main()
