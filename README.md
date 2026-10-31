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

* `monogenic` is 1 if $\mathcal{O}_K$ is monogenic and 0 if it provably is not.
  A field only appears in a file once this is decided.
* `zk_coeffs` are the coefficients (constant term first, as in `nf_fields.coeffs`) of the minimal
  polynomial of a generator $\alpha$, or `\N` (null) when `monogenic` is 0.

Since $\alpha$ is only determined up to $\alpha \mapsto \pm\alpha + k$ and a choice of class,
the stored $\alpha$ is:

1. the root of the polredabs polynomial, if it generates $\mathcal{O}_K$ (so `zk_coeffs = coeffs`);
2. otherwise, the generator of smallest $T_2$ norm, ties broken by the lexicographically
   smallest minimal polynomial (highest coefficient first). Within $\pm\alpha + k$ this means the
   trace lies in $[0, n/2]$.

Every row with `monogenic = 1` is certified by `scripts/check_data.py`: the polynomial's discriminant
equals the field discriminant (so its root generates the maximal order), and optionally `nfisisom`
against the LMFDB polynomial.

## Methods

* **Magma** ([`magma/monogenic.m`](magma/monogenic.m)): `IndexFormEquation(O, 1)` returns all
  generators up to equivalence for fields of degree at most 4, so it decides monogenicity and finds
  the canonical $\alpha$ exactly. See the
  [Magma handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/418).
* **PARI brute force** ([`scripts/monogenic.gp`](scripts/monogenic.gp)): search
  $\alpha = \sum_{i \ge 2} c_i \omega_i$ over the LLL-reduced integral basis with
  $\max |c_i| = 1, 2, \dots, B$, computing $[\mathcal{O}_K : \mathbb{Z}[\alpha]]$ as a determinant.
  This works in any degree but can only show that a field *is* monogenic, and the generator it
  finds is the smallest one it saw, which may not be the global minimum. Use it where Magma
  does not apply (degree $\ge 5$), and prefer Magma's rows where both exist.

| Degree | Status |
|---|---|
| 2 | done for $\lvert D \rvert \le 1000$ (`generate_quadratic.py`): every quadratic field is monogenic |
| 3, 4 | scripts ready (Magma, PARI) |
| $\ge 5$ | PARI brute force only |

## Scripts

```
python3 scripts/generate_quadratic.py --max-disc 1000              # data/degree2.txt
sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6  # LMFDB fields, Magma for degree <= 4
sage -python scripts/find_generators.py --degrees 5 --method pari -B 4
sage -python scripts/find_generators.py --degrees 3 --input fields.txt --print   # "label|[coeffs]" lines
python3 scripts/merge.py rows.txt [--replace]                       # merge rows from a manual Magma run
python3 scripts/check_data.py [--fields fields.txt]                 # validate every data/degree<n>.txt
python3 -m unittest discover tests                                  # run the tests (need cypari2)
```

`find_generators.py` reads the fields from the LMFDB database (the `lmfdb` package, or
[lmfdb-lite](https://github.com/roed314/lmfdb-lite):
`sage -pip install -U "lmfdb-lite[pgbinary] @ git+https://github.com/roed314/lmfdb-lite.git"`)
and merges its results into `data/` every `--batch` fields. Runs can be split with `--shard i/k`
or discriminant ranges and run in parallel; `--only-new` skips fields that already have a row.
Rows from Magma replace rows from the PARI search, since Magma's are exact.
