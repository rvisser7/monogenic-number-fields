"""Shared helpers: PARI interface and the data file format."""

import os
import re

import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**8, silent=True)
pari('read("%s")' % os.path.join(os.path.dirname(os.path.abspath(__file__)), "monogenic.gp"))

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
HEADER = "label|monogenic|zk_coeffs\ntext|smallint|numeric[]\n\n"
NULL = r"\N"


def parse_coeffs(s):
    """'[1, -1, 1]' or '1,-1,1' -> [1, -1, 1]"""
    return [int(c) for c in s.strip().strip("[]").split(",") if c.strip()]


def format_coeffs(coeffs):
    return "[" + ", ".join(str(c) for c in coeffs) + "]"


def label_disc(label):
    """The discriminant of the field with LMFDB label n.r.D.i."""
    n, r, D, _ = (int(x) for x in label.split("."))
    return (-1) ** ((n - r) // 2) * D


def label_key(label):
    """Sort by degree, |disc|, signature, then index, as the LMFDB does."""
    n, r, D, i = (int(x) for x in label.split("."))
    return (n, D, r, i)


def find_generator(coeffs, bound=3):
    """Ascending coefficients of the minimal polynomial of the canonical generator
    of O_K, K = Q[x]/(f) with f given by ascending coeffs, or None if the search up
    to height `bound` found nothing (which does not prove non-monogenicity)."""
    res = pari("mono_generator")(pari.nfinit(pari.Polrev(coeffs)), bound)
    if res == 0:
        return None
    return [int(c) for c in pari.Vecrev(res)]


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


def write_data(path, rows):
    """Write rows (label, monogenic, coeffs or None), sorted like the LMFDB."""
    with open(path, "w") as fh:
        fh.write(HEADER)
        for label, mono, coeffs in sorted(rows, key=lambda row: label_key(row[0])):
            fh.write("%s|%d|%s\n" % (label, mono, NULL if coeffs is None else format_coeffs(coeffs)))


def data_path(n):
    return os.path.join(DATA_DIR, "degree%d.txt" % n)


def merge_data(path, rows, replace=False):
    """Merge rows {label: (label, monogenic, coeffs)} into a data file. Existing rows
    are kept unless replace is set. Returns (added, replaced, kept)."""
    old = {row[0]: row for row in read_data(path)} if os.path.exists(path) else {}
    added = replaced = kept = 0
    for label, row in rows.items():
        if label not in old:
            added += 1
        elif replace and old[label] != row:
            replaced += 1
        else:
            kept += 1
            continue
        old[label] = row
    write_data(path, old.values())
    return added, replaced, kept
