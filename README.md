# Generators of rings of integers for monogenic number fields

This repository computes, for monogenic number fields $K$ in the [LMFDB](https://www.lmfdb.org/NumberField/),
an explicit $\alpha \in \mathcal{O}_K$ with $\mathcal{O}_K = \mathbb{Z}[\alpha]$.
See [LMFDB issue #6783](https://github.com/LMFDB/lmfdb/issues/6783).

For many monogenic fields the root of the polredabs defining polynomial already generates $\mathcal{O}_K$,
but not always: for [3.1.244.1](https://www.lmfdb.org/NumberField/3.1.244.1), with polredabs polynomial
$x^3 + x - 6$, the root has index 2, while a root of $x^3 - x^2 - 4x + 6$ generates $\mathcal{O}_K$.

## Data

The data lives in [`data/`](data/), with one file per degree: `degree2.txt`, `degree3.txt`, ...
Each file is in psycodict upload format:

```
label|monogenic|zk_coeffs
text|smallint|numeric[]

2.0.3.1|1|[1,-1,1]
2.0.4.1|1|[1,0,1]
```

* `monogenic` is 1 if $\mathcal{O}_K$ is monogenic, -1 if it provably is not, and 0 if this has not
  been computed. A field only appears in a file once this is decided, so fields missing from a
  file have `monogenic = 0`.
* `zk_coeffs` are the coefficients (constant term first, as in `nf_fields.coeffs`) of the minimal
  polynomial of a generator $\alpha$, or `\N` (null) when `monogenic` is not 1.

Since $\alpha$ is only determined up to $\alpha \mapsto \pm\alpha + k$ and a choice of class,
the stored $\alpha$ is:

1. the root of the polredabs polynomial, if it generates $\mathcal{O}_K$ (so `zk_coeffs = coeffs`);
2. otherwise, the generator whose minimal polynomial $P$ comes first in the order polredabs uses:
   smallest $T_2(P) = \sum_i |\alpha_i|^2$ (over all $n$ roots), then smallest $S(P)$ in
   lexicographic order, where $S(P) = (|a_1|, a_1, \dots, |a_n|, a_n)$ for
   $P = x^n + a_1x^{n-1} + \dots + a_n$. (Polredabs compares $|\operatorname{disc}(P)|$ in between,
   but that is $|\operatorname{disc}(K)|$ for every generator of $\mathcal{O}_K$.)

Rule 1 is the special case of rule 2 where the polredabs root generates $\mathcal{O}_K$, since the
polredabs polynomial comes first in the same order among all defining polynomials of $K$; it is
kept as a separate rule because $T_2$ is compared after rounding to 6 decimals.

Every row is checked by `scripts/check_data.py`: for `monogenic = 1`, PARI computes the ring of
integers of the field the polynomial defines and checks that $[\mathcal{O}_K : \mathbb{Z}[\alpha]] = 1$,
that its discriminant matches the label, and (with `--source db`) that it is the LMFDB field. With
`--cross-check B`, the PARI search is run as an independent check of the Magma results: it must find
no generator for a field marked -1, and none smaller than the stored one for a field marked 1.

## Methods

* **Polredabs root** (`--method polredabs`, any degree): the quick bulk pass. The root $\alpha$ of
  the LMFDB polynomial $f$ generates $\mathcal{O}_K$ iff $\operatorname{disc}(f) = \operatorname{disc}(K)$,
  since $\operatorname{disc}(f) = [\mathcal{O}_K : \mathbb{Z}[\alpha]]^2 \operatorname{disc}(K)$. That costs
  one polynomial discriminant per field (about $10^5$ fields a second), and these rows are already
  canonical (rule 1 above), so run this first and leave the other fields to the methods below. By
  rule 1 these replace any other stored generator, and a field marked -1 is a conflict (see
  below). With `--lmfdb-index` it reads only
  the fields whose LMFDB `index` column is 1, which is faster still; each is still checked.
* **Magma** ([`magma/monogenic.m`](magma/monogenic.m)): `IndexFormEquation(O, 1)` returns all
  generators up to equivalence for fields of degree at most 4, so it decides monogenicity and finds
  the canonical $\alpha$ exactly. See the
  [Magma handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/418).
* **PARI Thue solver** (`--method thue`, degree 3 only): the index form of a cubic field is a
  binary cubic form, so $\mathcal{O}_K = \mathbb{Z}[\alpha]$ is a Thue equation $I(x, y) = \pm 1$,
  which PARI's `thue` solves completely and unconditionally. This decides monogenicity and finds
  the canonical $\alpha$ exactly, like Magma, and gives the same rows (checked on all 21,845 cubic
  fields with $\lvert D \rvert \le 10^5$).
* **Factoring the index form** (`--method factor`, meant for degree 4): the index form $I(x)$
  (see below) is a polynomial over $\mathbb{Z}$, so $I(x) = \pm 1$ forces each of its factors to
  be $\pm 1$. For a quartic field $K$ with a quadratic subfield $k$ the cubic resolvent is
  reducible and $I$ factors (Gaál, Pethő and Pohst): one factor is the norm form of $k$ on
  $\mathcal{O}_K / (\mathcal{O}_k + \mathbb{Z})$, a binary quadratic form. When $k$ is imaginary
  it is positive definite, so it equals 1 at only finitely many points (found with `qfminim`),
  and at each of them $I(x) = \pm 1$ is a polynomial equation in the one remaining coordinate,
  solved exactly with `nfroots`. This decides monogenicity and finds the canonical $\alpha$
  exactly, like Magma, for every quartic field with an imaginary quadratic subfield (the
  totally complex V4 fields, and the D4 fields whose quadratic subfield is imaginary); other
  fields are skipped, unless the root of their polynomial generates. E.g. 4.0.1197.2 has no
  common index divisor or local obstruction, but
  $I = (x_3^2 - x_3 x_4 + x_4^2) \cdot Q(x_2, x_3, x_4)$, and $Q = \pm 1$ has no solution on the
  six points where the first factor is 1, so it is not monogenic.

  When $k$ is real (any signature) it can still prove that $K$ is *not* monogenic. Write
  $I = c \prod_l (U - t_l V)$, with $t_l$ the roots of the cubic resolvent and $U$, $V$ two
  ternary quadratic forms with rational coefficients (Gaál, Pethő and Pohst; found numerically,
  then the identity is checked exactly). A rational root $t$ splits this as $I = c \, L \, M$,
  $L = U - tV$ linear and $M$ quadratic in $(U, V)$, and $I = \pm 1$ forces $L$ and $M$ to
  take one of two values each. That is a line meeting a conic in the $(U, V)$ plane: at most two
  points for each choice of signs, and they must be rational with the right denominators. Often
  there are none (e.g. 4.2.1600.1: $5Q = (5x_2^2 - 2B)^2 + 16N^2$, so $\lvert I \rvert \ge 16/5$);
  otherwise $U(x) = u$, $V(x) = v$ is tested modulo small prime powers.

  *Totally complex quartic fields, all of them.* Without a quadratic subfield the resolvent is
  irreducible and $c \, v^3 P(u/v) = \pm 1$ is a cubic Thue equation, which PARI's `thue` solves
  completely and unconditionally, so again only finitely many $(U, V)$ are allowed. If $K$ is
  totally complex, with $\alpha_1, \alpha_2$ and $\alpha_3, \alpha_4$ complex conjugate,
  $\phi_1 - \phi_2 = \lvert \alpha_1 - \bar\alpha_3 \rvert^2$ and
  $\phi_1 - \phi_3 = \lvert \alpha_1 - \alpha_3 \rvert^2$ are $\ge 0$ (with
  $\phi_1 = \alpha_1\alpha_2 + \alpha_3\alpha_4$ and so on), so some $U - sV$ is a definite
  ternary form, and $U(x) = u$, $V(x) = v$ leaves finitely many $x$ (`qfminim`). So for every
  totally complex quartic field this decides monogenicity and finds the canonical $\alpha$
  exactly, like Magma. E.g. 4.0.3636.1 (S4): the Thue equation has two solutions, and neither
  comes from any $x$, so it is not monogenic. For the other signatures it only proves
  non-monogenicity; monogenic fields whose root does not generate are left to the search.

  Checks (fields from `nflist`): it decides all 1,813 totally complex quartic fields with
  $\lvert D \rvert \le 30000$ (779 monogenic), giving the stored generator for every one in
  `data/`; of the 5,039 C4, V4 and D4 fields with $\lvert D \rvert \le 10^5$ (every signature)
  it decided 4,968 before the Thue step was added; of the 866 mixed and totally real S4
  fields with $\lvert D \rvert \le 20000$ it decides 770, and 91 of the 96 left are monogenic.
  It never contradicts a row in `data/` or the search, except that for 4.0.23517.2 and
  4.0.76432.2 it finds a more canonical generator than the stored one (same $T_2$, smaller
  $S(P)$).
* **PARI brute force** ([`scripts/monogenic.gp`](scripts/monogenic.gp)): search
  $\alpha = \sum_{i \ge 2} c_i \omega_i$ over the LLL-reduced integral basis with
  $\max |c_i| = 1, 2, \dots, B$, computing $[\mathcal{O}_K : \mathbb{Z}[\alpha]]$ as a determinant.
  This works in any degree but can only show that a field *is* monogenic, and the generator it
  finds is the smallest one it saw, which may not be the global minimum. Use it where Magma
  does not apply (degree $\ge 5$), and prefer Magma's rows where both exist. When it finds
  nothing, the tests below are tried.
* **Proofs of non-monogenicity** (`--method local`, any degree): with
  $\alpha = \sum x_i \omega_i$, $\mathcal{O}_K = \mathbb{Z}[\alpha]$ iff the index form
  $I(x) = \det(1, \alpha, \dots, \alpha^{n-1})$, of degree $d = n(n-1)/2$, equals $\pm 1$.
  So $K$ is not monogenic if $I(x) = 1$ and $I(x) = -1$ each have no solution modulo some integer:
  1. *Common index divisors*: by Hensel's criterion, a prime $p$ divides every
     $[\mathcal{O}_K : \mathbb{Z}[\alpha]]$ iff, for some $f$, more primes above $p$ have residue
     degree $f$ than there are monic irreducible polynomials of degree $f$ over $\mathbb{F}_p$
     (e.g. 2 splitting completely in a cubic field). Only $p < n$ can do this.
  2. *Local obstructions*: $I(x) \equiv \pm 1 \pmod{p^k}$ has no solution, found by running over
     the $x$ modulo $p^k$. Since $I(tx) = t^d I(x)$, one modulus per prime decides solvability
     over $\mathbb{Z}_p$: $p$ if $p \nmid d$ (and only $\gcd(d, p - 1) > 1$ can fail), $p^{v+1}$ if
     $v = v_p(d) > 0$ ($2^{v+2}$ for $p = 2$). For example $I = 5x^3 - 2y^3$ for
     $\mathbb{Q}(\sqrt[3]{20})$ is never $\pm 1$ mod 9. Primes up to `--maxp` (50) are tried, and
     moduli needing more than `--max-points` ($10^5$) evaluations are skipped.

  Any obstruction found is a proof, so these never mark a monogenic field -1 (checked on all the
  monogenic fields in `data/`), but they can miss: some fields are not monogenic for global
  reasons only, e.g. for cubic fields with $\lvert D \rvert \le 10^5$ they prove 3,863 of the
  11,290 non-monogenic fields. For degree 3 and 4 use Magma or the Thue solver to settle the rest.

| Degree | Status |
|---|---|
| all | polredabs roots that generate $\mathcal{O}_K$, from `--method polredabs` |
| 2 | done for $\lvert D \rvert \le 1000$ (`generate_quadratic.py`): every quadratic field is monogenic |
| 3 | exact (Magma, or PARI Thue solver) |
| 4 | exact with Magma; exact with PARI (`--method factor`) for all totally complex fields, and it proves most other non-monogenic fields |
| $\ge 5$ | PARI brute force for 1, common index divisors and local obstructions for -1 |

## Scripts

```
python3 scripts/generate_quadratic.py --max-disc 1000              # data/degree2.txt
sage -python scripts/find_generators.py --degrees 2-47 --method polredabs   # every field whose polredabs root generates O_K
sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6  # LMFDB fields, Magma for degree <= 4
sage -python scripts/find_generators.py --degrees 5 --method pari -B 4
sage -python scripts/find_generators.py --degrees 3 --input fields.txt --print   # "label|[coeffs]" lines
sage -python scripts/find_generators.py --degrees 3 --method thue --max-disc 10^5   # exact, no Magma
sage -python scripts/find_generators.py --degrees 5-8 --method local --only-new  # prove the rest are not monogenic
sage -python scripts/find_generators.py --degrees 4 --method factor --only-new   # exact for totally complex quartics
sage -python scripts/find_generators.py --degrees 6 --method pari --lmfdb-monogenic -B 4  # LMFDB says monogenic, no generator here yet
python3 scripts/merge.py rows.txt [--replace]                       # merge rows from a manual Magma run
python3 scripts/check_data.py [--fix]                               # validate every data/degree<n>.txt
sage -python scripts/check_data.py --source db --cross-check 2      # also against the LMFDB and the PARI search
sage -python scripts/check_data.py --source db --reprove            # redo degree 3 with Thue, degree 4 by factoring, find proofs for -1 rows
sage -python scripts/compare_lmfdb.py --write compare/              # compare with nf_fields.monogenic in the LMFDB
python3 -m unittest discover tests                                  # run the tests (need cypari2)
```

`find_generators.py` reads the fields from the LMFDB database (the `lmfdb` package, or
[lmfdb-lite](https://github.com/roed314/lmfdb-lite):
`sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"`)
and merges its results into `data/` every `--batch` fields. Runs can be split with `--shard i/k`
or discriminant ranges and run in parallel; `--only-new` skips fields that already have a row.

### Which generator is kept

A field keeps one generator, the first in the canonical order above. When any method (or
`merge.py`) produces a different generator for a field that already has one, the new one replaces
it only if it comes first in that order, i.e. it is the polredabs root, or has smaller $T_2$ norm
(then smaller $S(P)$), and it passes the independent check below. Each replacement is
printed as `BETTER ...`. So a smaller generator from Magma or the Thue solver replaces one from the
PARI search (which only searches up to height `-B`), but a later search never replaces a smaller one.
`--replace` lets any checked generator replace the stored one. To bring the stored generators up to
date after a change to the canonical order, rerun a method on just the fields whose generator is not
the LMFDB polynomial, e.g. `find_generators.py --degrees 3 --method thue --recheck-generators`.

### Conflicts: a generator always beats -1

Magma's `IndexFormEquation` (or any of the methods here) could have a bug and wrongly say a field
is not monogenic, but a generator is a certificate that anyone can check. So whenever a new row
says 1 and the stored row says -1, or the other way round, the merge (`find_generators.py` with
any method, and `merge.py`) checks the generator independently, with or without `--replace`: PARI
computes the ring of integers of the field it defines (`nfdisc`, factoring the discriminant
completely) and checks that $[\mathcal{O}_K : \mathbb{Z}[\alpha]] = 1$ and that the discriminant
matches the label. If it passes, the 1 row is written or kept; otherwise the -1 row is. Either way
the conflict is printed as `CONFLICT ...` and appended to `conflicts.txt`, with the time, both rows
and the command that produced the new one, so it can be reported (e.g. to Magma) or looked at
again. To search again for generators of every field marked -1:

```
sage -python scripts/find_generators.py --degrees 3-4 --method pari --recheck-negative -B 4
```

`check_data.py --cross-check B` does the same search but only reports, without changing `data/`.

## Comparing with the LMFDB

The LMFDB table `nf_fields` has its own `monogenic` column, with the same values (1, -1, 0, and
also null, treated as 0). `scripts/compare_lmfdb.py` compares it with `data/` and prints, for each
degree, how many fields each side marks 1 and -1, how many agree, and

* **conflicts**: one side says 1 and the other -1. If ours is 1, the stored `zk_coeffs` (checked by
  `check_data.py`) is a generator, so the LMFDB value is wrong; if ours is -1, one of the proofs
  needs a second look. The script exits with status 1 if there are any;
* fields where only we know the answer (an update for the LMFDB), and fields where only the LMFDB
  does (still to do here).

`--write DIR` writes the full lists: `conflicts.txt`, `ours_only.txt` (`label|monogenic`),
`lmfdb_only_yes.txt` and `lmfdb_only_no.txt` (`label|[coeffs]`, ready for
`find_generators.py --input`) and `unexpected.txt` (LMFDB values other than 1, -1, 0, null).
`--lmfdb FILE` reads `label|monogenic` lines from a file instead of the database. `--min-disc` and
`--max-disc` (e.g. `--degrees 6 --max-disc 10^4`) compare only the fields with $\lvert D \rvert$ in
that range, on both sides, for degrees with millions of fields.

To run the search on just the fields the LMFDB says are monogenic but that have no generator
here yet (no row, or a row marked 0 or -1), straight from the database:

```
sage -python scripts/find_generators.py --degrees 6 --method pari --lmfdb-monogenic -B 4
```

This adds `monogenic = 1` to the `nf_fields` query and skips the fields already marked 1. A field
marked -1 here is a conflict with the LMFDB; if the search finds a generator it is verified and
replaces the -1 as above (add `--only-new` to leave those fields out).
