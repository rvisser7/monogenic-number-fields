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
2. otherwise, the generator of smallest $T_2$ norm, ties broken by the lexicographically
   smallest minimal polynomial (highest coefficient first). Within $\pm\alpha + k$ this means the
   trace lies in $[0, n/2]$.

Every row is checked by `scripts/check_data.py`: for `monogenic = 1`, PARI computes the ring of
integers of the field the polynomial defines and checks that $[\mathcal{O}_K : \mathbb{Z}[\alpha]] = 1$,
that its discriminant matches the label, and (with `--source db`) that it is the LMFDB field. With
`--cross-check B`, the PARI search is run as an independent check of the Magma results: it must find
no generator for a field marked -1, and none smaller than the stored one for a field marked 1.

## Methods

* **Magma** ([`magma/monogenic.m`](magma/monogenic.m)): `IndexFormEquation(O, 1)` returns all
  generators up to equivalence for fields of degree at most 4, so it decides monogenicity and finds
  the canonical $\alpha$ exactly. See the
  [Magma handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/418).
* **PARI Thue solver** (`--method thue`, degree 3 only): the index form of a cubic field is a
  binary cubic form, so $\mathcal{O}_K = \mathbb{Z}[\alpha]$ is a Thue equation $I(x, y) = \pm 1$,
  which PARI's `thue` solves completely and unconditionally. This decides monogenicity and finds
  the canonical $\alpha$ exactly, like Magma, and gives the same rows (checked on all 21,845 cubic
  fields with $\lvert D \rvert \le 10^5$).
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
| 2 | done for $\lvert D \rvert \le 1000$ (`generate_quadratic.py`): every quadratic field is monogenic |
| 3 | exact (Magma, or PARI Thue solver) |
| 4 | exact with Magma |
| $\ge 5$ | PARI brute force for 1, common index divisors and local obstructions for -1 |

## Scripts

```
python3 scripts/generate_quadratic.py --max-disc 1000              # data/degree2.txt
sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6  # LMFDB fields, Magma for degree <= 4
sage -python scripts/find_generators.py --degrees 5 --method pari -B 4
sage -python scripts/find_generators.py --degrees 3 --input fields.txt --print   # "label|[coeffs]" lines
sage -python scripts/find_generators.py --degrees 3 --method thue --max-disc 10^5   # exact, no Magma
sage -python scripts/find_generators.py --degrees 5-8 --method local --only-new  # prove the rest are not monogenic
python3 scripts/merge.py rows.txt [--replace]                       # merge rows from a manual Magma run
python3 scripts/check_data.py [--fix]                               # validate every data/degree<n>.txt
sage -python scripts/check_data.py --source db --cross-check 2      # also against the LMFDB and the PARI search
sage -python scripts/check_data.py --source db --reprove            # redo degree 3 with Thue, find proofs for -1 rows
python3 -m unittest discover tests                                  # run the tests (need cypari2)
```

`find_generators.py` reads the fields from the LMFDB database (the `lmfdb` package, or
[lmfdb-lite](https://github.com/roed314/lmfdb-lite):
`sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"`)
and merges its results into `data/` every `--batch` fields. Runs can be split with `--shard i/k`
or discriminant ranges and run in parallel; `--only-new` skips fields that already have a row.
Rows from Magma or the Thue solver replace rows from the PARI search, since they are exact.
