#!/usr/bin/env python3
"""Check the integral bases (nf_fields.zk) and the index column stored in the LMFDB.

For a field K = Q(a) of degree n with LMFDB polynomial f (coeffs) and discriminant D
(read off the label), the stored basis b_1, ..., b_n (zk, polynomials in a) is the ring
of integers O_K exactly when
  1. it has n elements, which are linearly independent;
  2. each b_i is an algebraic integer: its characteristic polynomial (of b_i in
     Q[a]/(f)) has integer coefficients;
  3. the lattice they span has discriminant D: disc(b_1, ..., b_n) = det(B)^2 disc(f) = D,
     where B is the matrix of the b_i on 1, a, ..., a^(n-1).
Given 2, the lattice is a sublattice of O_K, of some index m, so its discriminant is
m^2 D, and 3 says m = 1. The `index` column, [O_K : Z[a]], must equal sqrt(disc(f) / D). All of this is exact (PARI), and it does not need to
compute O_K. With --nfdisc it also checks that D is PARI's discriminant of K (nfdisc,
which factors disc(f), so it is slow for large discriminants).

E.g. if the LMFDB gives 1, a, a^2, a^3 for 4.0.19982592.7 with f = x^4 + 252x^2 + 8673,
then disc(f) = 2401^2 * 19982592, so that basis has index 2401 in O_K (it misses
(a^2 + 126)/49 = sqrt(3)).

Fields come from the LMFDB database (`from lmfdb import db`, or lmfdb-lite, as in
find_generators.py): coeffs and index from nf_fields, zk from nf_fields_extra (joined on
the label, --batch labels per query). Or they come from a file (--input) with lines
    label|[coeffs]|zk|index
where zk is a JSON list of strings in a, e.g.
    4.0.19982592.7|[8673,0,252,0,1]|["1","a","1/49*a^2 + 18/7","1/49*a^3 + 25/7*a"]|2401
and index may be \\N. Header lines are skipped.

Each problem is printed as "label: ..." (at most --show per degree) and, with --write
FILE, written as label|kind|message. The kinds are: format (zk missing, the wrong
length, or not readable), dependent, not_integral, disc (the basis is not all of O_K),
index (the index column disagrees), nfdisc. The exit status is 1 if there is a problem.

Usage:
    sage -python scripts/check_lmfdb_zk.py --degrees 4 --max-disc 10^6
    sage -python scripts/check_lmfdb_zk.py --degrees 2-47 --write zk_problems.txt
    sage -python scripts/check_lmfdb_zk.py --degrees 6 --shard 0/4    # and 1/4, ...
    python3 scripts/check_lmfdb_zk.py --input nf_fields_zk.txt
"""

import argparse
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_lmfdb import DiscRange, Log, parse_degrees, parse_int
from monogenic import NULL, label_disc, label_key, parse_coeffs, pari

KINDS = ["format", "dependent", "not_integral", "disc", "index", "nfdisc"]


def parse_element(e, f):
    """A basis element as stored (a string in a, e.g. '1/49*a^2 + 18/7', or a list of
    rational coefficients on 1, a, a^2, ...) -> a PARI polynomial in x of degree < n,
    reduced mod f."""
    if isinstance(e, (list, tuple)):
        g = pari.Polrev([pari(str(c)) for c in e]) if e else pari(0)
    else:
        g = pari.subst(pari(str(e)), "a", "x")
    if g.type() not in ("t_INT", "t_FRAC", "t_POL") or \
            (g.type() == "t_POL" and pari.variables(g) != pari("[x]")):
        raise ValueError("not a polynomial in a: %r" % (e,))
    return g % f if g.type() == "t_POL" else g


def check_zk(label, coeffs, zk, index=None, nfdisc=False):
    """Problems [(kind, message)] with the LMFDB basis zk (and index column) of the field
    with this label and polynomial (ascending coeffs). Empty if zk is an integral basis."""
    n = int(label.split(".")[0])
    D = label_disc(label)
    f = pari.Polrev(coeffs)
    out = []
    if nfdisc and pari.nfdisc(f) != D:
        out.append(("nfdisc", "PARI's disc(K) is %s, not %d as in the label" % (pari.nfdisc(f), D)))
    ind = pari.sqrtint(pari.poldisc(f) // D) if pari.poldisc(f) % D == 0 else None
    if ind is None or ind ** 2 * D != pari.poldisc(f):
        out.append(("disc", "disc(f) = %s is not a square times %d" % (pari.poldisc(f), D)))
        ind = None
    if index is not None and ind is not None and int(index) != int(ind):
        out.append(("index", "the index column is %s, but [O_K : Z[a]] = sqrt(disc(f) / D) = %s"
                    % (index, ind)))
    if zk is None:
        return out + [("format", "no zk")]
    if len(zk) != n:
        return out + [("format", "zk has %d elements, not %d" % (len(zk), n))]
    try:
        elts = [parse_element(e, f) for e in zk]
    except Exception as err:
        return out + [("format", "cannot read zk: %s" % err)]
    B = pari.matrix(n, n, [pari.polcoef(g, i, "x") for i in range(n) for g in elts])
    det = pari.matdet(B)
    if det == 0:
        return out + [("dependent", "the zk elements are linearly dependent")]
    bad = [i + 1 for i, g in enumerate(elts)
           if pari.denominator(pari.content(pari.charpoly(pari.Mod(g, f)))) != 1]
    if bad:
        out.append(("not_integral", "zk element%s %s (%s) %s not an algebraic integer"
                    % ("s" if len(bad) > 1 else "", ", ".join(map(str, bad)),
                       ", ".join(str(zk[i - 1]) for i in bad), "is" if len(bad) == 1 else "are")))
    disc = det ** 2 * pari.poldisc(f)
    if disc != D:
        msg = "the basis has discriminant %s, not %d" % (disc, D)
        if not bad and disc % D == 0 and pari.issquare(disc // D):
            msg += ", so it has index %s in O_K" % pari.sqrtint(disc // D)
        out.append(("disc", msg))
    return out


def get_db():
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    return db


def db_fields(n, disc, log, shard, batch=10000, db=None):
    """(label, coeffs, zk, index) for the LMFDB fields of degree n: coeffs and index from
    nf_fields, zk from nf_fields_extra (same label), fetched `batch` labels at a time."""
    db = db or get_db()
    query = dict(degree=n, **disc.query())
    total = db.nf_fields.count(query)
    log.say("degree %d: %d LMFDB fields%s" % (n, total, disc))
    rows = []

    def flush():
        zk = {r["label"]: r.get("zk") for r in
              db.nf_fields_extra.search({"label": {"$in": [r[0] for r in rows]}}, ["label", "zk"])}
        out = [(label, coeffs, zk.get(label), index) for label, coeffs, index in rows]
        rows.clear()
        return out

    for i, r in enumerate(db.nf_fields.search(query, ["label", "coeffs", "index"])):
        if i % shard[1] != shard[0]:
            continue
        rows.append((r["label"], [int(c) for c in r["coeffs"]], r.get("index")))
        if len(rows) >= batch:
            yield from flush()
    if rows:
        yield from flush()


def file_fields(path, degrees, disc, shard):
    """{n: [(label, coeffs, zk, index)]} from a file of label|[coeffs]|zk|index lines."""
    out = {}
    with open(path) as fh:
        lines = [line.rstrip("\n").split("|") for line in fh]
    i = 0
    for cells in lines:
        if len(cells) < 3 or not cells[0][:1].isdigit():
            continue  # blank or header line
        label = cells[0]
        n = int(label.split(".")[0])
        if (degrees and n not in degrees) or label not in disc:
            continue
        i += 1
        if (i - 1) % shard[1] != shard[0]:
            continue
        zk = None if cells[2] == NULL else json.loads(cells[2])
        index = int(cells[3]) if len(cells) > 3 and cells[3] not in (NULL, "") else None
        out.setdefault(n, []).append((label, parse_coeffs(cells[1]), zk, index))
    for rows in out.values():
        rows.sort(key=lambda r: label_key(r[0]))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", help="e.g. 3-4 or 3,5 (required with the database)")
    parser.add_argument("--input", metavar="FILE", help="read the fields from this file instead of the LMFDB")
    parser.add_argument("--min-disc", type=parse_int, help="only fields with |disc| >= this, e.g. 10^4")
    parser.add_argument("--max-disc", type=parse_int, help="only fields with |disc| <= this, e.g. 10^6")
    parser.add_argument("--nfdisc", action="store_true",
                        help="also check the label's discriminant with PARI's nfdisc (slow)")
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th field, from the i-th")
    parser.add_argument("--batch", type=int, default=10000,
                        help="labels per nf_fields_extra query (default 10000)")
    parser.add_argument("--show", type=int, default=20, metavar="N",
                        help="print at most N problems per degree (default 20)")
    parser.add_argument("--write", metavar="FILE", help="write every problem to FILE (label|kind|message)")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="print progress this often (default 30)")
    parser.add_argument("-q", "--quiet", action="store_true", help="no progress lines")
    args = parser.parse_args()
    if not args.input and not args.degrees:
        sys.exit("--degrees is required when reading from the LMFDB")
    shard = tuple(int(x) for x in args.shard.split("/"))
    disc = DiscRange(args.min_disc, args.max_disc)
    log = Log(args.progress_every, args.quiet)
    degrees = parse_degrees(args.degrees) if args.degrees else []
    from_file = file_fields(args.input, degrees, disc, shard) if args.input else None
    out = open(args.write, "w") if args.write else None
    any_problem = False
    for n in (sorted(from_file) if from_file is not None else degrees):
        rows = from_file[n] if from_file is not None else db_fields(n, disc, log, shard, args.batch)
        checked, kinds, bad, shown = 0, Counter(), 0, 0
        for label, coeffs, zk, index in rows:
            checked += 1
            problems = check_zk(label, coeffs, zk, index, args.nfdisc)
            if problems:
                bad += 1
                any_problem = True
            for kind, msg in problems:
                kinds[kind] += 1
                if shown < args.show:
                    print("%s: %s" % (label, msg), flush=True)
                    shown += 1
                if out:
                    out.write("%s|%s|%s\n" % (label, kind, msg))
            log.tick("degree %d: %d fields checked, %d with problems" % (n, checked, bad))
        print("degree %d: %d fields checked, %d with problems%s"
              % (n, checked, bad, "".join("; %d %s" % (kinds[k], k) for k in KINDS if kinds[k])), flush=True)
    if out:
        out.close()
    sys.exit(1 if any_problem else 0)


if __name__ == "__main__":
    main()
