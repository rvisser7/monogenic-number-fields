"""Write data/degree2.txt. Every quadratic field is monogenic, and the root of the
polredabs polynomial (x^2 - x + (1 - D)/4 or x^2 - D/4) generates O_K.

    python3 scripts/generate_quadratic.py --max-disc 1000
"""

import argparse
import os

from monogenic import DATA_DIR, find_generator, pari, write_data


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--max-disc", type=int, default=1000)
    args = parser.parse_args()

    rows = []
    for D in range(-args.max_disc, args.max_disc + 1):
        if D in (0, 1) or not pari.isfundamental(D):
            continue
        f = pari.polredabs(pari("x^2 - (%d)" % D))
        label = "2.%d.%d.1" % (2 if D > 0 else 0, abs(D))
        rows.append((label, 1, find_generator([int(c) for c in pari.Vecrev(f)])))
    write_data(os.path.join(DATA_DIR, "degree2.txt"), rows)
    print("wrote %d quadratic fields" % len(rows))


if __name__ == "__main__":
    main()
