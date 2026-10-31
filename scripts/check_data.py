"""Validate every data/degree<n>.txt.

For each row with monogenic = 1 this checks that zk_coeffs is a monic irreducible
polynomial of degree n whose discriminant equals the field discriminant read off the
label, so its root generates the maximal order of the field it defines. Rows with
monogenic = 0 must have a null zk_coeffs.

With --fields (lines "label|[polredabs coeffs]") it also checks with nfisisom that
the polynomial defines the field with that label, which the discriminant alone
does not decide (e.g. 4.0.1040.1 and 4.0.1040.2).
"""

import argparse
import glob
import os
import sys

from monogenic import DATA_DIR, label_disc, label_key, pari, parse_coeffs, read_data


def check_file(path, fields):
    errors = []
    n = int(os.path.basename(path)[len("degree"):-len(".txt")])
    rows = read_data(path)
    labels = [row[0] for row in rows]
    if labels != sorted(labels, key=label_key) or len(set(labels)) != len(labels):
        errors.append("rows are not sorted or not unique")
    for label, mono, coeffs in rows:
        if int(label.split(".")[0]) != n:
            errors.append("%s: wrong degree file" % label)
        elif mono == 0:
            if coeffs is not None:
                errors.append("%s: monogenic = 0 but zk_coeffs given" % label)
        elif mono != 1 or coeffs is None:
            errors.append("%s: monogenic must be 0 or 1, with zk_coeffs iff 1" % label)
        else:
            f = pari.Polrev(coeffs)
            if len(coeffs) != n + 1 or coeffs[-1] != 1 or not pari.polisirreducible(f):
                errors.append("%s: not a monic irreducible polynomial of degree %d" % (label, n))
            elif pari.poldisc(f) != label_disc(label):
                errors.append("%s: discriminant %s, expected %s" % (label, pari.poldisc(f), label_disc(label)))
            elif label in fields and not pari.nfisisom(f, pari.Polrev(fields[label])):
                errors.append("%s: polynomial defines a different field" % label)
    return len(rows), errors


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fields", help='file with lines "label|[coeffs]"')
    args = parser.parse_args()
    fields = {}
    if args.fields:
        with open(args.fields) as fh:
            for line in fh:
                if line.strip():
                    label, coeffs = line.strip().split("|")
                    fields[label] = parse_coeffs(coeffs)

    ok = True
    for path in sorted(glob.glob(os.path.join(DATA_DIR, "degree*.txt")), key=lambda p: (len(p), p)):
        count, errors = check_file(path, fields)
        print("%s: %d rows, %d errors" % (os.path.basename(path), count, len(errors)))
        for e in errors:
            print("  " + e)
        ok = ok and not errors
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
