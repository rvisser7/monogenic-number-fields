"""Add quadratic fields to data/degree2.txt. Every quadratic field is monogenic, and the
root of the polredabs polynomial (x^2 - x + (1 - D)/4 or x^2 - D/4) generates O_K.

    python3 scripts/generate_quadratic.py --max-disc 1000
    python3 scripts/generate_quadratic.py --min-disc 10^6 --max-disc 2*10^6

The rows for the fundamental discriminants D with min-disc <= |D| <= max-disc are merged
into data/degree2.txt as by find_generators.py (see monogenic.merge_data), so the rows
already there, including fields outside the range, are kept.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import data_path, merge_data, pari


def parse_int(s):
    """An integer, also written as 10^6 or 2*10^6."""
    return int(pari(s))


def quadratic_rows(min_disc, max_disc):
    """{label: (label, 1, coeffs)} for the quadratic fields with min_disc <= |D| <= max_disc."""
    rows = {}
    for D in range(-max_disc, max_disc + 1):
        if abs(D) < min_disc or D in (0, 1) or not pari.isfundamental(D):
            continue
        f = pari.polredabs(pari("x^2 - (%d)" % D))
        label = "2.%d.%d.1" % (2 if D > 0 else 0, abs(D))
        rows[label] = (label, 1, [int(c) for c in pari.Vecrev(f)])
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--min-disc", type=parse_int, default=1, help="smallest |D| (default 1)")
    parser.add_argument("--max-disc", type=parse_int, default=1000, help="largest |D| (default 1000)")
    args = parser.parse_args()

    rows = quadratic_rows(args.min_disc, args.max_disc)
    added, replaced, kept, conflicts = merge_data(data_path(2), rows,
                                                  source="generate_quadratic.py " + " ".join(sys.argv[1:]))
    print("%d quadratic fields: %d added, %d replaced, %d kept" % (len(rows), added, replaced, kept)
          + ("; %d CONFLICTS (see conflicts.txt)" % conflicts if conflicts else ""))


if __name__ == "__main__":
    main()
