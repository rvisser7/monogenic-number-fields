#!/usr/bin/env python3
"""Check the inessential primes (nf_fields.inessentialp) stored in the LMFDB.

The inessential primes of a number field K (also called common index divisors) are the
primes p that divide [O_K : Z[alpha]] for every alpha in O_K with K = Q(alpha), i.e. the
primes dividing the index of the field i(K) (the `index` column). This script recomputes
them independently of the LMFDB, with Hensel's criterion:

    p divides i(K)  iff  for some f, more primes of K above p have residue degree f
                         than there are monic irreducible polynomials of degree f over F_p.

(Only the residue degrees matter, not the ramification indices.) The number of monic
irreducible polynomials of degree f over F_p is (1/f) sum_{d | f} mu(d) p^(f/d), and the
residue degrees come from PARI's idealprimedec, on an order that is maximal at the primes
p < n = [K : Q] (only those can be inessential, since the count for f = 1 is p). So
e.g. 2 is inessential in the cubic field 3.3.961.1, where 2 splits completely and there
are only 2 linear polynomials over F_2.

For each field it checks:
  * inessentialp is the set of primes computed as above (kind `inessentialp`), or reports
    that it is missing (kind `missing`, null in the LMFDB, shown there as "Not computed");
  * the primes dividing the `index` column are the same primes, and `index` divides
    [O_K : Z[a]] = sqrt(disc(f) / disc(K)) for the LMFDB polynomial f (kind `index`; the
    exact power of each prime is not checked);
  * a field with an inessential prime is not monogenic, so neither the LMFDB `monogenic`
    column (kind `monogenic`) nor our data/degree<n>.txt (kind `ours`, unless --no-data)
    may say 1 for it. Every row with 1 in data/ has a verified generator, so `ours` would
    point at a bug in this script.

Fields come from the LMFDB database (`from lmfdb import db`, or lmfdb-lite, as in
find_generators.py), all columns from nf_fields. Or they come from a file (--input) with
lines
    label|[coeffs]|inessentialp|index|monogenic
where the last two columns may be left out or be \\N, and inessentialp is a list like [2]
or [] (\\N if missing), e.g.
    3.3.961.1|[8,-10,-1,1]|[2]|2|-1
Header lines are skipped.

Each problem is printed as "label: ..." (at most --show per degree) and, with --write FILE,
written as label|kind|message. With --update FILE, the computed inessentialp of every field
where the LMFDB value is wrong or missing goes to FILE, for the LMFDB editors, with header
label|inessentialp, text|integer[]. The exit status is 1 if there is a problem.

Usage:
    sage -python scripts/check_lmfdb_inessential.py --degrees 3-8 --max-disc 10^6
    sage -python scripts/check_lmfdb_inessential.py --degrees 2-47 --write inessential.txt --update inessential_update.txt
    sage -python scripts/check_lmfdb_inessential.py --degrees 6 --shard 0/4    # and 1/4, ...
    python3 scripts/check_lmfdb_inessential.py --input nf_fields_inessential.txt
"""

import argparse
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_lmfdb import DiscRange, Log, parse_degrees, parse_int
from monogenic import NULL, data_path, label_disc, label_key, parse_coeffs, pari, read_data

KINDS = ["missing", "inessentialp", "index", "monogenic", "ours"]


def mobius(n):
    out, p = 1, 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            out = -out
        p += 1
    return -out if n > 1 else out


def n_irreducible(p, f):
    """The number of monic irreducible polynomials of degree f over F_p."""
    return sum(mobius(d) * p ** (f // d) for d in range(1, f + 1) if f % d == 0) // f


def inessential_primes(coeffs):
    """The inessential primes (common index divisors) of K = Q[x]/(f), f given by its
    ascending coeffs, by Hensel's criterion (see the top of this file)."""
    n = len(coeffs) - 1
    ps = [int(p) for p in pari.primes([2, n - 1])] if n > 2 else []
    if not ps:
        return []
    nf = pari.nfinit([pari.Polrev(coeffs), ps])  # maximal at the primes in ps
    out = []
    for p in ps:
        degrees = Counter(int(P[3]) for P in pari.idealprimedec(nf, p))  # P[3] = P.f
        if any(k > n_irreducible(p, f) for f, k in degrees.items()):
            out.append(p)
    return out


def parse_primes(value):
    """inessentialp as stored: a list of integers (or strings), or a string such as
    '[2, 3]', '{2,3}' or '2,3'; None or \\N (not computed) -> None."""
    if value is None or value == NULL:
        return None
    if isinstance(value, str):
        value = [c for c in value.strip().strip("[]{}").split(",") if c.strip()]
    return sorted(int(p) for p in value)


def fmt(primes):
    return "{%s}" % ", ".join(map(str, primes)) if primes else "none"


def check_field(label, coeffs, inessentialp, index=None, monogenic=None, ours=None):
    """(problems, computed): problems [(kind, message)] with the LMFDB columns of the
    field with this label and polynomial (ascending coeffs), computed its inessential
    primes. inessentialp as stored (see parse_primes), index and monogenic the LMFDB
    values or None, ours our monogenic value (data/) or None."""
    cid = inessential_primes(coeffs)
    out = []
    stored = parse_primes(inessentialp)
    if stored is None:
        out.append(("missing", "inessentialp is not computed; it should be %s" % fmt(cid)))
    elif stored != cid:
        out.append(("inessentialp", "inessentialp is %s, but the inessential primes are %s"
                    % (fmt(stored), fmt(cid))))
    if index is not None:
        index = int(index)
        if index < 1:
            out.append(("index", "the index column is %d, not a positive integer" % index))
        else:
            primes = [int(p) for p in pari.factor(index)[0]] if index > 1 else []
            if primes != cid:
                out.append(("index", "the index column is %d, but the inessential primes are %s"
                            % (index, fmt(cid))))
            f, D = pari.Polrev(coeffs), label_disc(label)
            q = pari.poldisc(f) / D
            if q.type() == "t_INT" and pari.issquare(q) and int(pari.sqrtint(q)) % index != 0:
                out.append(("index", "the index column is %d, which does not divide [O_K : Z[a]] = %s"
                            % (index, pari.sqrtint(q))))
    if cid and monogenic is not None and int(monogenic) == 1:
        out.append(("monogenic", "the LMFDB says monogenic, but %s %s inessential"
                    % (fmt(cid), "is" if len(cid) == 1 else "are")))
    if cid and ours == 1:
        out.append(("ours", "data/ says monogenic, but %s %s inessential"
                    % (fmt(cid), "is" if len(cid) == 1 else "are")))
    return out, cid


def get_db():
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db
    return db


def db_fields(n, disc, log, shard, db=None):
    """(label, coeffs, inessentialp, index, monogenic) for the LMFDB fields of degree n."""
    db = db or get_db()
    query = dict(degree=n, **disc.query())
    log.say("degree %d: %d LMFDB fields%s" % (n, db.nf_fields.count(query), disc))
    cols = ["label", "coeffs", "inessentialp", "index", "monogenic"]
    for i, r in enumerate(db.nf_fields.search(query, cols)):
        if i % shard[1] == shard[0]:
            yield (r["label"], [int(c) for c in r["coeffs"]], r.get("inessentialp"),
                   r.get("index"), r.get("monogenic"))


def file_fields(path, degrees, disc, shard):
    """{n: [(label, coeffs, inessentialp, index, monogenic)]} from a file."""
    out = {}
    i = 0
    with open(path) as fh:
        for line in fh:
            cells = line.rstrip("\n").split("|")
            if len(cells) < 3 or not cells[0][:1].isdigit():
                continue  # blank or header line
            label = cells[0]
            n = int(label.split(".")[0])
            if (degrees and n not in degrees) or label not in disc:
                continue
            i += 1
            if (i - 1) % shard[1] != shard[0]:
                continue
            index, mono = [int(c) if c not in (NULL, "") else None for c in (cells[3:5] + ["", ""])[:2]]
            out.setdefault(n, []).append((label, parse_coeffs(cells[1]), cells[2], index, mono))
    for rows in out.values():
        rows.sort(key=lambda r: label_key(r[0]))
    return out


def our_data(n):
    """{label: monogenic} from data/degree<n>.txt ({} if there is no such file)."""
    path = data_path(n)
    return {label: mono for label, mono, _ in read_data(path)} if os.path.exists(path) else {}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", help="e.g. 3-8 or 3,5 (required with the database)")
    parser.add_argument("--input", metavar="FILE", help="read the fields from this file instead of the LMFDB")
    parser.add_argument("--min-disc", type=parse_int, help="only fields with |disc| >= this, e.g. 10^4")
    parser.add_argument("--max-disc", type=parse_int, help="only fields with |disc| <= this, e.g. 10^6")
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th field, from the i-th")
    parser.add_argument("--no-data", action="store_true", help="do not compare with data/degree<n>.txt")
    parser.add_argument("--show", type=int, default=20, metavar="N",
                        help="print at most N problems per degree (default 20)")
    parser.add_argument("--write", metavar="FILE", help="write every problem to FILE (label|kind|message)")
    parser.add_argument("--update", metavar="FILE",
                        help="write label|inessentialp with the right value for every wrong or missing one")
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
    update = open(args.update, "w") if args.update else None
    if update:
        update.write("label|inessentialp\ntext|integer[]\n\n")
    any_problem = False
    for n in (sorted(from_file) if from_file is not None else degrees):
        rows = from_file[n] if from_file is not None else db_fields(n, disc, log, shard)
        ours = {} if args.no_data else our_data(n)
        checked, kinds, bad, shown, with_cid = 0, Counter(), 0, 0, 0
        for label, coeffs, inessentialp, index, mono in rows:
            checked += 1
            problems, cid = check_field(label, coeffs, inessentialp, index, mono, ours.get(label))
            with_cid += bool(cid)
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
            if update and any(kind in ("missing", "inessentialp") for kind, _ in problems):
                update.write("%s|[%s]\n" % (label, ",".join(map(str, cid))))
            log.tick("degree %d: %d fields checked, %d with problems" % (n, checked, bad))
        print("degree %d: %d fields checked, %d with inessential primes, %d with problems%s"
              % (n, checked, with_cid, bad, "".join("; %d %s" % (kinds[k], k) for k in KINDS if kinds[k])),
              flush=True)
    for fh in (out, update):
        if fh:
            fh.close()
    sys.exit(1 if any_problem else 0)


if __name__ == "__main__":
    main()
