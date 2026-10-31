"""Merge rows "label|monogenic|zk_coeffs" (from find_generators.py or
magma/monogenic.m) into data/degree<n>.txt.

    python3 scripts/merge.py rows.txt [more_rows.txt ...] [--overwrite]

Existing rows are kept unless --overwrite is given. Rows from Magma's index form
equation solver are exact (they consider all generators), so merge them with
--overwrite over rows from the brute-force search.
"""

import argparse
import os
from collections import defaultdict

from monogenic import DATA_DIR, NULL, parse_coeffs, read_data, write_data


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("rows", nargs="+")
    parser.add_argument("--overwrite", action="store_true")
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
        path = os.path.join(DATA_DIR, "degree%d.txt" % n)
        old = {row[0]: row for row in read_data(path)} if os.path.exists(path) else {}
        changed = 0
        for label, row in rows.items():
            if label not in old or (args.overwrite and old[label] != row):
                old[label] = row
                changed += 1
        write_data(path, old.values())
        print("degree%d.txt: %d rows added or changed, %d total" % (n, changed, len(old)))


if __name__ == "__main__":
    main()
