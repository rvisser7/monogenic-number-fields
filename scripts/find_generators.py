"""Brute-force search for generators of O_K with PARI.

Reads fields as lines "label|[a0, a1, ..., an]" (ascending polredabs coefficients,
as in the LMFDB column nf_fields.coeffs), or from the LMFDB database directly, and
writes a row "label|1|[...]" for each field where a generator was found. Fields
where nothing was found are listed on stderr: they may still be monogenic.

    python3 scripts/find_generators.py fields.txt -o rows.txt -B 4
    sage -python scripts/find_generators.py --source db --degree 3 --max-disc 10^6 -o rows.txt

Merge the rows into data/ with scripts/merge.py.
"""

import argparse
import sys

from monogenic import find_generator, format_coeffs, parse_coeffs


def fields_from_file(path):
    with open(path) as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                label, coeffs = line.strip().split("|")
                yield label, parse_coeffs(coeffs)


def fields_from_db(degree, max_disc):
    from lmfdb import db  # needs an LMFDB checkout on the path
    query = {"degree": degree}
    if max_disc is not None:
        query["disc_abs"] = {"$lte": max_disc}
    for rec in db.nf_fields.search(query, ["label", "coeffs"]):
        yield rec["label"], [int(c) for c in rec["coeffs"]]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", help='file with lines "label|[coeffs]"')
    parser.add_argument("--source", choices=["file", "db"], default="file")
    parser.add_argument("--degree", type=int)
    parser.add_argument("--max-disc", type=lambda s: int(eval(s.replace("^", "**"))))
    parser.add_argument("-B", "--bound", type=int, default=3,
                        help="max |coordinate| on the integral basis (default 3)")
    parser.add_argument("-o", "--output", help="output file (default stdout)")
    args = parser.parse_args()

    if args.source == "db":
        if args.degree is None:
            parser.error("--source db needs --degree")
        fields = fields_from_db(args.degree, args.max_disc)
    elif args.input:
        fields = fields_from_file(args.input)
    else:
        parser.error("give an input file or --source db")

    out = open(args.output, "w") if args.output else sys.stdout
    for label, coeffs in fields:
        gen = find_generator(coeffs, args.bound)
        if gen is None:
            print("not found: %s" % label, file=sys.stderr)
        else:
            out.write("%s|1|%s\n" % (label, format_coeffs(gen)))
            out.flush()


if __name__ == "__main__":
    main()
