
## Monogenic data

Data on monogenic number fields is given in this folder.

There is one file per degree: `degree2.txt`, `degree3.txt`, ...
Each file is in psycodict upload format:

```
label|monogenic|monogenic_coeffs
text|smallint|numeric[]

2.0.3.1|1|[1,-1,1]
2.0.4.1|1|[1,0,1]
```

* `monogenic` is 1 if $\mathcal{O}_K$ is monogenic, -1 if it provably is not, and 0 if this has not
  been computed. A field only appears in a file once this is decided, so fields missing from a
  file have `monogenic = 0`.
* `monogenic_coeffs` are the coefficients (constant term first, as in `nf_fields.coeffs`) of the minimal
  polynomial of a generator $\alpha$, or `\N` (null) when `monogenic` is not 1.

Since $\alpha$ is only determined up to $\alpha \mapsto \pm\alpha + k$ and a choice of class,
the stored $\alpha$ is:

1. the root of the polredabs polynomial, if it generates $\mathcal{O}_K$ (so `monogenic_coeffs = coeffs`);
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
