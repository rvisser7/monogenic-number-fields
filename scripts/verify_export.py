"""Independently verify a file written by export_lmfdb.py (label|monogenic_coeffs rows).

This script uses only Sage (and PARI through Sage), and none of the code in this
repository, so that it can serve as an independent check of the data.

For each row "label|[c0,c1,...,cn]", with f = c0 + c1 x + ... + cn x^n, alpha a root of f
and K the LMFDB field with that label (defining polynomial g, discriminant disc(K)),
it checks that

  1. f is monic of degree n, where n is the degree in the label n.r.D.i;
  2. f is irreducible over Q;
  3. Q(alpha) is isomorphic to K: g has a root beta in Q[x]/(f). PARI (nfroots) finds
     beta, and Sage checks g(beta) = 0 in Q[x]/(f) exactly. Since g is irreducible of the
     same degree n, Q(beta) = Q(alpha), so this does not depend on trusting nfroots
     (if f = g there is nothing to do);
  4. disc(f) = disc(K), where disc(K) is disc_sign * disc_abs in the LMFDB and also
     (-1)^((n-r)/2) * D read off the label.

Since disc(f) = [O_K : Z[alpha]]^2 disc(K), 3 and 4 show that Z[alpha] = O_K, i.e. that
alpha generates the ring of integers of the field with that label.

The fields come from the LMFDB (the lmfdb package, or lmfdb-lite: `from lmf import db`),
or with --fields from a file of "label|[g0,...,gn]" lines (e.g. exported from nf_fields;
disc(K) is then the discriminant read off the label). With --no-fields, step 3 is skipped
and step 4 only compares with the label: that still shows that Z[alpha] is the full ring
of integers of Q(alpha) and that disc(Q(alpha)) matches the label, but not which of the
fields with that discriminant Q(alpha) is.

    sage -python scripts/verify_export.py monogenic_coeffs.txt
    sage -python scripts/verify_export.py monogenic_coeffs.txt --fields fields.txt
    sage -python scripts/verify_export.py monogenic_coeffs.txt --no-fields

Every failing row is printed as "label|reason"; the exit code is 1 if any row fails.
"""

import argparse
import sys
import time

from sage.all import QQ, ZZ, PolynomialRing, pari

R = PolynomialRing(ZZ, "x")
S = PolynomialRing(QQ, "x")


def parse_coeffs(s):
    s = s.strip()
    if not (s.startswith("[") and s.endswith("]")):
        raise ValueError("not a list: %s" % s)
    return [ZZ(c) for c in s[1:-1].split(",")]


def label_parts(label):
    n, r, D, i = (ZZ(x) for x in label.split("."))
    return n, r, D, i


def label_disc(label):
    n, r, D, _ = label_parts(label)
    return (-1) ** ((n - r) // 2) * D


def check(label, coeffs, field=None):
    """None if the row is correct, else the reason it is not. field is (g, disc(K)) with g
    the coefficients of the LMFDB polynomial, or None to skip the isomorphism check."""
    n, r, D, _ = label_parts(label)
    if not (0 <= r <= n and (n - r) % 2 == 0 and D > 0):
        return "bad label"
    f = R(coeffs)
    if f.degree() != n:
        return "degree %s, not %s" % (f.degree(), n)
    if f.leading_coefficient() != 1:
        return "f is not monic"
    if not f.is_irreducible():
        return "f is reducible"
    disc_f = f.discriminant()
    if disc_f != label_disc(label):
        return "disc(f) = %s, but the label says %s" % (disc_f, label_disc(label))
    if field is not None:
        g_coeffs, disc_K = field
        if disc_f != disc_K:
            return "disc(f) = %s, but disc(K) = %s in the LMFDB" % (disc_f, disc_K)
        g = R(g_coeffs)
        if g.degree() != n or not g.is_irreducible():
            return "the LMFDB polynomial is not irreducible of degree %s" % n
        if g == f:  # alpha is a root of g
            return None
        # PARI wants the variable of the field (y) to have lower priority than x.
        roots = pari.nfroots(pari(f).subst("x", "y"), pari(g))
        if len(roots) == 0:
            return "Q(alpha) is not isomorphic to K (g has no root in Q(alpha))"
        beta = S([QQ(c) for c in pari(roots[0]).lift().Vecrev()])
        if g.change_ring(QQ)(beta) % f.change_ring(QQ) != 0:
            return "Q(alpha) is not isomorphic to K (bad root of g)"
    return None


def read_export(path):
    """(label, coeffs string) for each row of an exported file."""
    with open(path) as fh:
        header = [next(fh).rstrip("\n") for _ in range(3)]
        if header != ["label|monogenic_coeffs", "text|numeric[]", ""]:
            raise ValueError("%s: bad header %s" % (path, header))
        for line in fh:
            line = line.strip()
            if line:
                label, coeffs = line.split("|")
                yield label, coeffs


def read_fields(path):
    """{label: (g, disc(K))} from a file of label|[g0,...,gn] lines."""
    fields = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and "|" in line and line[0].isdigit():
                label, g = line.split("|")[:2]
                fields[label] = (parse_coeffs(g), label_disc(label))
    return fields


def lmfdb_fields(labels):
    """{label: (g, disc(K))} from the LMFDB table nf_fields."""
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    out = {}
    for rec in db.nf_fields.search({"label": {"$in": list(labels)}},
                                   ["label", "coeffs", "disc_sign", "disc_abs"]):
        out[rec["label"]] = ([ZZ(c) for c in rec["coeffs"]],
                             ZZ(int(rec["disc_sign"])) * ZZ(int(rec["disc_abs"])))
    return out


def batches(rows, size):
    batch = []
    for row in rows:
        batch.append(row)
        if len(batch) == size:
            yield batch
            batch = []
    if batch:
        yield batch


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("export", help="the file written by export_lmfdb.py")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--fields", help="read the LMFDB polynomials from this file "
                                        "(label|[g0,...,gn] lines) instead of the LMFDB")
    group.add_argument("--no-fields", action="store_true",
                       help="skip the isomorphism check (no LMFDB polynomials)")
    parser.add_argument("--batch", type=int, default=2000, help="labels per LMFDB query")
    parser.add_argument("--progress-every", type=int, default=100000)
    args = parser.parse_args()

    fields = read_fields(args.fields) if args.fields else None
    checked = failed = 0
    seen = set()
    start = time.time()
    for batch in batches(read_export(args.export), args.batch):
        if not args.no_fields and fields is None:
            known = lmfdb_fields(label for label, _ in batch)
        else:
            known = fields
        for label, coeffs in batch:
            if label in seen:
                reason = "label occurs twice"
            elif not args.no_fields and label not in known:
                reason = "label not in the LMFDB"
            else:
                try:
                    reason = check(label, parse_coeffs(coeffs),
                                   None if args.no_fields else known[label])
                except Exception as err:  # a malformed row
                    reason = "error: %s" % err
            seen.add(label)
            checked += 1
            if reason is not None:
                failed += 1
                print("%s|%s" % (label, reason), flush=True)
            if args.progress_every and checked % args.progress_every == 0:
                print("... %d rows checked, %d failed, %.0f s" % (checked, failed,
                      time.time() - start), file=sys.stderr, flush=True)
    print("%d rows checked, %d failed%s" % (
        checked, failed, " (isomorphism to the LMFDB fields not checked)" if args.no_fields
        else ""), file=sys.stderr)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
