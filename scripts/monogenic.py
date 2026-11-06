"""Shared helpers: PARI interface and the data file format."""

import os
import sys
import time

import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**8, silent=True)
pari('read("%s")' % os.path.join(os.path.dirname(os.path.abspath(__file__)), "monogenic.gp"))

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
# merge_data appends a line here for every field where a new row says 1 and the stored
# row says -1, or the other way round (see merge_data).
CONFLICT_LOG = os.path.join(os.path.dirname(DATA_DIR), "conflicts.txt")
HEADER = "label|monogenic|zk_coeffs\ntext|smallint|numeric[]\n\n"
NULL = r"\N"


def parse_coeffs(s):
    """'[1, -1, 1]' or '1,-1,1' -> [1, -1, 1]"""
    return [int(c) for c in s.strip().strip("[]").split(",") if c.strip()]


def format_coeffs(coeffs):
    return "[" + ",".join(str(c) for c in coeffs) + "]"


def label_disc(label):
    """The discriminant of the field with LMFDB label n.r.D.i."""
    n, r, D, _ = (int(x) for x in label.split("."))
    return (-1) ** ((n - r) // 2) * D


def label_key(label):
    """Sort by degree, |disc|, signature, then index, as the LMFDB does."""
    n, r, D, i = (int(x) for x in label.split("."))
    return (n, D, r, i)


def root_generates(coeffs, label):
    """Whether a root of f (ascending coeffs) generates O_K, for K = Q[x]/(f) the field
    with this label: disc(f) = [O_K : Z[alpha]]^2 disc(K), so this is disc(f) = disc(K)."""
    return int(pari.poldisc(pari.Polrev(coeffs))) == label_disc(label)


def find_generator(coeffs, bound=3):
    """Ascending coefficients of the minimal polynomial of the canonical generator
    of O_K, K = Q[x]/(f) with f given by ascending coeffs, or None if the search up
    to height `bound` found nothing (which does not prove non-monogenicity)."""
    res = pari("mono_generator")(pari.nfinit(pari.Polrev(coeffs)), bound)
    if res == 0:
        return None
    return [int(c) for c in pari.Vecrev(res)]


def cubic_generator(coeffs):
    """For a cubic field, the same as find_generator, but exact: PARI solves the index
    form equation (a Thue equation) completely, so None means O_K is provably not
    monogenic."""
    res = pari("mono_cubic_generator")(pari.nfinit(pari.Polrev(coeffs)))
    if res == 0:
        return None
    return [int(c) for c in pari.Vecrev(res)]


def factor_generator(coeffs):
    """Decide monogenicity by factoring the index form (mono_factor_generator in
    monogenic.gp): exact for every totally complex quartic field (and for quartic fields
    with an imaginary quadratic subfield), and for other quartic fields it can prove
    non-monogenicity.
    Returns (True, generator coeffs) if O_K is monogenic, with the canonical generator,
    (True, None) if it is provably not, and (False, None) if undecided."""
    res = pari("mono_factor_generator")(pari.nfinit(pari.Polrev(coeffs)))
    if res == -1:
        return False, None
    if res == 0:
        return True, None
    return True, [int(c) for c in pari.Vecrev(res)]


def obstruction(coeffs, maxp=50, max_points=10**5):
    """A proof that O_K is not monogenic, K = Q[x]/(f), or None if none was found:
    "cid p" if p is a common index divisor, or "local m1 m2" if the index form takes
    neither the value 1 mod m1 nor the value -1 mod m2 (see monogenic.gp). Tries primes
    up to maxp, skipping moduli that need more than max_points evaluations."""
    res = pari("mono_obstruction")(pari.nfinit(pari.Polrev(coeffs)), maxp, max_points)
    if res == 0:
        return None
    return " ".join(str(c) for c in res)


def generator_key(coeffs):
    """The key by which the canonical generator is the smallest (see monogenic.gp),
    for the generator with minimal polynomial given by ascending coeffs."""
    return pari("mono_polkey")(pari.Polrev(coeffs))


def is_smaller(key1, key2):
    return pari.lex(key1, key2) < 0


def is_polredabs(coeffs):
    """Whether f (ascending coeffs) is its own polredabs polynomial. For a generator of
    O_K this says that it is the root of the LMFDB polynomial, which is rule 1 of the
    canonical choice (see README)."""
    f = pari.Polrev(coeffs)
    return pari.polredabs(f) == f


def more_canonical(coeffs, other):
    """Whether the generator with minimal polynomial `coeffs` comes before the one with
    `other` in the canonical order (README): the root of the polredabs polynomial first,
    then the polredabs order: smallest T2 norm, then smallest S(P). Both must generate
    O_K for the same field."""
    if coeffs == other:
        return False
    p1, p2 = is_polredabs(coeffs), is_polredabs(other)
    if p1 != p2:
        return p1
    return is_smaller(generator_key(coeffs), generator_key(other))


def read_data(path):
    """Rows (label, monogenic, coeffs or None) of a data file."""
    with open(path) as fh:
        lines = fh.read().split("\n")
    if "\n".join(lines[:3]) + "\n" != HEADER:
        raise ValueError("%s: bad header" % path)
    rows = []
    for line in lines[3:]:
        if not line:
            continue
        label, mono, coeffs = line.split("|")
        rows.append((label, int(mono), None if coeffs == NULL else parse_coeffs(coeffs)))
    return rows


def format_data(rows):
    """The contents of a data file with rows (label, monogenic, coeffs or None),
    sorted like the LMFDB."""
    return HEADER + "".join(
        "%s|%d|%s\n" % (label, mono, NULL if coeffs is None else format_coeffs(coeffs))
        for label, mono, coeffs in sorted(rows, key=lambda row: label_key(row[0])))


def write_data(path, rows):
    with open(path, "w") as fh:
        fh.write(format_data(rows))


def data_path(n):
    return os.path.join(DATA_DIR, "degree%d.txt" % n)


def verify_generator(coeffs, label):
    """Whether the root alpha of f (ascending coeffs) provably generates O_K, for K the
    field with this label: f is monic and irreducible of degree n, and PARI's ring of
    integers of L = Q(alpha) (nfdisc, which factors disc(f) completely) has discriminant
    disc(f) = disc(K), i.e. [O_L : Z[alpha]] = 1 and disc(L) = disc(K). This is
    independent of whatever produced f, so a row that passes it is a certificate."""
    n = int(label.split(".")[0])
    if coeffs is None or len(coeffs) != n + 1 or coeffs[-1] != 1:
        return False
    f = pari.Polrev(coeffs)
    if not pari.polisirreducible(f):
        return False
    D = label_disc(label)
    return pari.poldisc(f) == D and pari.nfdisc(f) == D


def log_conflict(message, log=None):
    """Print a conflict to stderr and append it, with the time, to the conflict log
    (default CONFLICT_LOG)."""
    print("CONFLICT " + message, file=sys.stderr, flush=True)
    log = log or CONFLICT_LOG
    if log:
        with open(log, "a") as fh:
            fh.write("%s %s\n" % (time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), message))


def merge_data(path, rows, replace=False, source=None, log=None):
    """Merge rows {label: (label, monogenic, coeffs)} into a data file.

    * A field without a row, or with monogenic = 0, gets the new row; a new row with
      monogenic = 0 never replaces a decided one.
    * If one row says 1 and the other -1 (a conflict: one of the methods has a bug, e.g.
      Magma's IndexFormEquation missing a solution), the generator is checked with
      verify_generator, whatever `replace` says. If it generates O_K, the 1 row is
      written (or kept), since that is a certificate; otherwise the -1 row is. Either
      way the conflict is printed and logged to `log` (default CONFLICT_LOG, conflicts.txt), with `source`
      (default: the command line) as the origin of the new row.
    * If both rows say 1 with different generators, the new one replaces the stored one
      only if it is more canonical (more_canonical: the polredabs root, or smaller T2
      norm, then smaller S(P)) and passes verify_generator; with replace,
      whenever it passes verify_generator. Each replacement is printed ("BETTER ..."); a
      more canonical generator that fails verification is a conflict and is logged.
    * Otherwise (both -1) the existing row is kept.

    Returns (added, replaced, kept, conflicts)."""
    old = {row[0]: row for row in read_data(path)} if os.path.exists(path) else {}
    if source is None:
        source = " ".join([os.path.basename(sys.argv[0])] + sys.argv[1:])
    name = os.path.basename(path)
    added = replaced = kept = conflicts = 0
    for label, row in rows.items():
        prev = old.get(label)
        if prev is None:
            added += 1
        elif prev == row or row[1] == 0:
            kept += 1
            continue
        elif prev[1] == 0:
            replaced += 1
        elif prev[1] != row[1]:
            conflicts += 1
            yes, no = (row, prev) if row[1] == 1 else (prev, row)
            ok = verify_generator(yes[2], label)
            log_conflict("%s: %s has %s, the new row (from %s) has %s; the generator %s, so stored %d"
                         % (label, name, format_row_value(prev), source, format_row_value(row),
                            "is verified ([O_K : Z[alpha]] = 1)" if ok else "FAILED verification",
                            1 if ok else -1), log)
            winner = yes if ok else no
            if winner is prev:
                kept += 1
                continue
            replaced += 1
        elif row[1] == 1:
            better = more_canonical(row[2], prev[2])
            if not (better or replace):
                kept += 1
                continue
            if not verify_generator(row[2], label):
                conflicts += 1
                log_conflict("%s: %s has %s, the new row (from %s) has %s, which comes first in the "
                             "canonical order but FAILED verification, so kept the stored one"
                             % (label, name, format_row_value(prev), source, format_row_value(row)), log)
                kept += 1
                continue
            print("%s %s: replaced %s by %s (from %s)"
                  % ("BETTER" if better else "REPLACED", label, format_coeffs(prev[2]),
                     format_coeffs(row[2]), source), file=sys.stderr, flush=True)
            replaced += 1
        else:
            kept += 1
            continue
        old[label] = row
    write_data(path, old.values())
    return added, replaced, kept, conflicts


def format_row_value(row):
    """'1 [coeffs]', '-1' or '0', for messages."""
    return "%d %s" % (row[1], format_coeffs(row[2])) if row[1] == 1 and row[2] is not None else str(row[1])
