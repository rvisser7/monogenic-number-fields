\\ Search for a generator alpha of O_K, i.e. O_K = Z[alpha].
\\
\\ Elements are written in the integral basis nf.zk returned by nfinit, which
\\ for a polredabs polynomial is the integral basis the LMFDB stores (zk).
\\
\\ Canonical choice (shared with magma/monogenic.m): if the root of the
\\ defining polynomial generates O_K, alpha is that root. Otherwise alpha is
\\ the generator of smallest T2 norm, ties broken by the lexicographically
\\ smallest minimal polynomial, after normalising alpha as in mono_normalize.

\\ Index [O_K : Z[alpha]] for alpha given by its coordinates v on nf.zk.
mono_index(nf, v) =
{
  my(n = poldegree(nf.pol), a = nfbasistoalg(nf, v), p = Mod(1, nf.pol), M);
  M = matrix(n, n);
  for (j = 1, n, M[, j] = nfalgtobasis(nf, p); p *= a);
  abs(matdet(M));
}

\\ Replace alpha by the "smallest" of the elements +-alpha + k (all of which
\\ generate the same order): smallest trace, which ends up in [0, n/2], then
\\ smallest minimal polynomial in lexicographic order. This also minimises T2
\\ over the class, since T2(a + k) = T2(a) + 2k Tr(a) + n k^2.
mono_normalize(nf, a) =
{
  my(n = poldegree(nf.pol), t, b);
  a -= floor(trace(a) / n);
  t = trace(a);
  b = if (t, 1 - a, -a);
  if (lex([trace(b), Vec(minpoly(b))], [t, Vec(minpoly(a))]) < 0, b, a);
}

\\ Sort key: T2 norm (rounded, to compare floats safely), then minimal polynomial.
mono_key(nf, a) = [round(norml2(nfeltembed(nf, a)) * 10^6), Vec(minpoly(a))];

\\ Try alpha = x first, then all sum_{i>=2} c_i zk[i] with max |c_i| = 1, 2,
\\ ..., B. At the first height where generators exist, return the smallest
\\ one by mono_key, as its minimal polynomial (t_POL).
\\ Returns 0 if none is found: that does NOT prove K is non-monogenic.
mono_generator(nf, B = 3) =
{
  my(n = poldegree(nf.pol), best = 0, bestkey, v, a, key);
  if (mono_index(nf, nfalgtobasis(nf, variable(nf.pol))) == 1, return(nf.pol));
  for (h = 1, B,
    forvec (w = vector(n - 1, i, [-h, h]),
      if (vecmax(abs(w)) != h, next);
      my(k = 1); while (!w[k], k++); if (w[k] < 0, next);
      v = concat([0], w)~;
      if (mono_index(nf, v) != 1, next);
      a = mono_normalize(nf, nfbasistoalg(nf, v));
      key = mono_key(nf, a);
      if (!best || lex(key, bestkey) < 0, best = a; bestkey = key));
    if (best, break));
  if (best, minpoly(best), 0);
}
