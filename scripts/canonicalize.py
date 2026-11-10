#!/usr/bin/env python3
"""Make sure the stored generators are the canonical ones, and replace those that are not.

The canonical generator of O_K (see README) is the root of the polredabs polynomial if it
generates O_K, and otherwise the generator alpha with smallest T2(alpha), then smallest
S(P) = (|a1|, a1, ..., |an|, an) for its minimal polynomial P. For each row marked 1 whose
generator alpha is not the polredabs root, this checks every element of O_K with T2 at most
T2(alpha) (mono_canonical in monogenic.gp: qfminim on the T2 form modulo Z), which proves
that alpha is canonical or finds the canonical generator. That replaces the stored one
through the usual merge (it must pass the independent check; printed as "BETTER ...").

This works offline: the stored generator defines the field, so it does not need the LMFDB.
For degree 3 and 4, rows where more than --max-points elements would have to be checked (e.g.
generators with huge coefficients) are decided with the exact methods instead (Thue, or
--method factor, which find every generator); in higher degree they are left as they are and
listed as "NOT PROVED CANONICAL label".

Usage:
    python3 scripts/canonicalize.py --degrees 5-7
    python3 scripts/canonicalize.py --degrees 6 --shard 0/4 --max-points 10^6
    python3 scripts/canonicalize.py --degrees 4 --max-disc 10^5 --print
"""

import argparse
import os
import sys
import time
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from monogenic import (canonical_generator, cubic_generator, data_path, factor_generator,
                       format_coeffs, is_polredabs, merge_data, pari, read_data)
from progress import Progress, fmt_time


def parse_degrees(s):
    out = []
    for part in s.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return sorted(set(d for d in out if d >= 3))


def parse_int(s):
    """An integer, also written as 10^6 or 2*10^6."""
    return int(pari(s))


def exact_generator(n, coeffs):
    """For degree 3 and 4, the canonical generator from the exact methods (Thue, factor),
    which find every generator: used when there are too many elements to check, e.g. for
    generators with huge coefficients. Returns (coeffs, proved)."""
    f = [int(c) for c in pari.Vecrev(pari.polredabs(pari.Polrev(coeffs)))]
    if n == 3:
        gen = cubic_generator(f)
        return (gen, True) if gen is not None else (coeffs, False)
    decided, gen = factor_generator(f)
    return (gen, True) if decided and gen is not None else (coeffs, False)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--degrees", required=True, help="e.g. 5-7 or 4,6")
    parser.add_argument("--min-disc", type=parse_int, help="smallest |disc|, e.g. 10^6")
    parser.add_argument("--max-disc", type=parse_int, help="largest |disc|, e.g. 10^8")
    parser.add_argument("--shard", default="0/1", help="i/k: only every k-th row, from the i-th")
    parser.add_argument("--max-points", type=parse_int, default=10**5,
                        help="the most elements to check for one field (default 10^5)")
    parser.add_argument("--batch", type=int, default=1000,
                        help="merge into the data file after this many rows (default 1000)")
    parser.add_argument("--print", action="store_true",
                        help="print the better generators instead of merging them")
    parser.add_argument("--progress-every", type=float, default=30, metavar="SECONDS",
                        help="when not on a terminal, how often to print progress (default 30)")
    args = parser.parse_args()

    shard, nshards = (int(x) for x in args.shard.split("/"))
    progress = Progress(every=args.progress_every)
    for n in parse_degrees(args.degrees):
        path = data_path(n)
        if not os.path.exists(path):
            continue
        t0 = time.time()
        rows = []
        for label, mono, coeffs in read_data(path):
            disc = int(label.split(".")[2])
            if mono != 1 or args.min_disc is not None and disc < args.min_disc \
                    or args.max_disc is not None and disc > args.max_disc:
                continue
            rows.append((label, coeffs))
        rows = rows[shard::nshards]
        count, batch, merged = Counter(), {}, Counter()

        def flush():
            if args.print:
                for label, _, gen in batch.values():
                    print("%s|1|%s" % (label, format_coeffs(gen)), flush=True)
            elif batch:
                merged.update(dict(zip(["added", "replaced", "kept", "conflicts"],
                                       merge_data(path, batch))))
            batch.clear()

        progress.start("degree %d" % n)
        for i, (label, coeffs) in enumerate(rows):
            if is_polredabs(coeffs):
                count["polredabs root"] += 1
            else:
                gen, proved = canonical_generator(coeffs, args.max_points)
                if not proved and n in (3, 4):
                    gen, proved = exact_generator(n, coeffs)
                if not proved:
                    count["not proved"] += 1
                    print("NOT PROVED CANONICAL %s: %s" % (label, format_coeffs(coeffs)),
                          file=sys.stderr, flush=True)
                elif gen != coeffs:
                    count["better"] += 1
                    batch[label] = (label, 1, gen)
                else:
                    count["canonical"] += 1
                if len(batch) >= args.batch:
                    flush()
            progress(i + 1, len(rows), ", %d better" % count["better"])
        flush()
        print("degree %d: %d generators: %d polredabs roots, %d proved canonical, %d replaced by a "
              "better one, %d not proved (more than --max-points elements)"
              % (n, len(rows), count["polredabs root"], count["canonical"], count["better"],
                 count["not proved"])
              + ("" if args.print else "; merged: %d replaced, %d kept" % (merged["replaced"], merged["kept"])
                 + ("; %d CONFLICTS (see conflicts.txt)" % merged["conflicts"] if merged["conflicts"] else ""))
              + ", in " + fmt_time(time.time() - t0), file=sys.stderr if args.print else sys.stdout)


if __name__ == "__main__":
    main()
