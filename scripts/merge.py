"""Merge rows "label|monogenic|zk_coeffs" (e.g. from magma/monogenic.m run by hand)
into data/degree<n>.txt.

    python3 scripts/merge.py rows.txt [more_rows.txt ...] [--replace]

Existing rows are kept unless --replace is given, except that when one row says 1 and
the other -1, the generator is verified and wins if it passes ([O_K : Z[alpha]] = 1),
and the conflict is printed and logged to conflicts.txt (see monogenic.merge_data).
"""

import argparse
from collections import defaultdict

from monogenic import NULL, data_path, merge_data, parse_coeffs


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rows", nargs="+")
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    new = defaultdict(dict)
    for path in args.rows:
        with open(path) as fh:
            for line in fh:
                if line.strip():
                    label, mono, coeffs = line.strip().split("|")
                    row = (label, int(mono), None if coeffs == NULL else parse_coeffs(coeffs))
                    new[int(label.split(".")[0])][label] = row

    for n, rows in sorted(new.items()):
        added, replaced, kept, conflicts = merge_data(data_path(n), rows, args.replace,
                                                      source="merge.py " + " ".join(args.rows))
        print("degree%d.txt: %d added, %d replaced, %d kept" % (n, added, replaced, kept)
              + ("; %d CONFLICTS (see conflicts.txt)" % conflicts if conflicts else ""))


if __name__ == "__main__":
    main()
