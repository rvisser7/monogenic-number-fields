# Generators of rings of integers for monogenic number fields

This repository computes, for monogenic number fields $K$ in the [LMFDB](https://www.lmfdb.org/NumberField/),
an explicit $\alpha \in \mathcal{O}_K$ with $\mathcal{O}_K = \mathbb{Z}[\alpha]$.
See [LMFDB issue #6783](https://github.com/LMFDB/lmfdb/issues/6783).

For many monogenic fields the root of the polredabs defining polynomial already generates $\mathcal{O}_K$,
but not always: for [3.1.244.1](https://www.lmfdb.org/NumberField/3.1.244.1), with polredabs polynomial
$x^3 + x - 6$, the root has index 2, while a root of $x^3 - x^2 - 4x + 6$ generates $\mathcal{O}_K$.
