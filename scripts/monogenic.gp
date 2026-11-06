\\ Search for a generator alpha of O_K, i.e. O_K = Z[alpha].
\\
\\ Elements are written in the integral basis nf.zk returned by nfinit, which
\\ for a polredabs polynomial is the integral basis the LMFDB stores (zk).
\\
\\ Canonical choice (shared with magma/monogenic.m): if the root of the
\\ defining polynomial generates O_K, alpha is that root. Otherwise alpha is
\\ the generator whose minimal polynomial P comes first in the polredabs order:
\\ smallest T2(P), then (|disc P| is disc(K) for every generator) smallest
\\ S(P) = (|a1|, a1, ..., |an|, an) for P = x^n + a1 x^(n-1) + ... + an. The
\\ polredabs root, when it generates O_K, is first in this order anyway.

\\ Index [O_K : Z[alpha]] for alpha given by its coordinates v on nf.zk.
mono_index(nf, v) =
{
  my(n = poldegree(nf.pol), a = nfbasistoalg(nf, v), p = Mod(1, nf.pol), M);
  M = matrix(n, n);
  for (j = 1, n, M[, j] = nfalgtobasis(nf, p); p *= a);
  abs(matdet(M));
}

\\ S(P) = [|a1|, a1, ..., |an|, an] for P = x^n + a1 x^(n-1) + ... + an, the
\\ last tie-break of polredabs.
mono_S(P) = my(v = Vec(P)); concat(vector(#v - 1, i, [abs(v[i + 1]), v[i + 1]]));

\\ Sort key, from the minimal polynomial g of alpha: T2(g) = sum |root|^2 over all
\\ n roots (rounded, to compare floats safely), then S(g).
mono_polkey(g) = concat([round(norml2(polroots(g)) * 10^6)], mono_S(g));

\\ The same key for alpha in nf. (Not from nfeltembed, which gives one embedding
\\ per pair of complex ones, so would count those once instead of twice.)
mono_key(nf, a) = mono_polkey(minpoly(a));

\\ Replace alpha by the smallest, by mono_key, of the elements +-alpha + k (all
\\ of which generate the same order). Since T2(a + k) = T2(a) + 2k Tr(a) + n k^2,
\\ T2 is smallest for k next to -Tr(a)/n, so there are four candidates.
mono_normalize(nf, a) =
{
  my(k = floor(trace(a) / poldegree(nf.pol)), best = 0, bestkey, key);
  foreach ([a - k, a - k - 1, -a + k, -a + k + 1], b,
    key = mono_key(nf, b);
    if (!best || lex(key, bestkey) < 0, best = b; bestkey = key));
  best;
}

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

\\ ---------------------------------------------------------------------------
\\ Proving that O_K is NOT monogenic.
\\
\\ For alpha = sum x_i zk[i], the index form I(x) = det(1, alpha, ..., alpha^(n-1)),
\\ with the powers written on nf.zk, is a form of degree d = n(n-1)/2 in
\\ x_2, ..., x_n with integer coefficients (x_1 drops out), and O_K = Z[alpha]
\\ iff I(x) = +-1. So K is not monogenic if I(x) = 1 and I(x) = -1 both have no
\\ solution modulo some integers m. Two tests of that kind:

\\ Number of monic irreducible polynomials of degree f over F_p.
mono_nirred(p, f) = sumdiv(f, e, moebius(e) * p^(f/e)) / f;

\\ (1) Common index divisors. By Hensel's criterion, p divides every index
\\ [O_K : Z[alpha]] iff for some f there are more primes above p of residue
\\ degree f than monic irreducible polynomials of degree f over F_p. Only
\\ primes p < n can do this. Returns the list of such p (empty if none).
mono_common_index_divisors(nf) =
{
  my(n = poldegree(nf.pol), res = List());
  forprime (p = 2, n - 1,
    my(dec = idealprimedec(nf, p), F = vector(n));
    for (i = 1, #dec, F[dec[i].f]++);
    for (f = 1, n, if (F[f] > mono_nirred(p, f), listput(res, p); break)));
  Vec(res);
}

\\ The modulus p^k that decides whether I(x) = s has a solution over Z_p, or 0
\\ if p can only be an obstruction as a common index divisor. Since I is
\\ homogeneous, I(t x) = t^d I(x), so the values of I on primitive x are unions
\\ of cosets of the d-th powers of units:
\\  * if p does not divide d, x -> x^d is onto 1 + pZ_p, so working mod p is
\\    enough, and if moreover gcd(d, p - 1) = 1, every unit is a d-th power and
\\    one unit value of I gives all of them;
\\  * if p | d, with v = v_p(d), the d-th powers contain 1 + p^(v+1) Z_p
\\    (1 + 2^(v+2) Z_2 for p = 2), so working mod p^(v+1) (2^(v+2)) is enough.
mono_local_modulus(p, d) =
{
  my(v = valuation(d, p));
  if (v, return(p^(v + 1 + (p == 2))));
  if (gcd(d, p - 1) > 1, p, 0);
}

\\ Multiplication matrices of nf.zk: T[i] * v is zk[i] * v, all on nf.zk.
mono_multables(nf) =
{
  my(n = poldegree(nf.pol), E = matid(n));
  vector(n, i, matconcat(vector(n, j, nfeltmul(nf, E[, i], E[, j]))));
}

\\ I(x) mod m, for x on nf.zk with x[1] = 0.
mono_indexform_mod(T, x, m) =
{
  my(n = #x, A = sum(i = 2, n, x[i] * T[i]) % m, M = matrix(n, n), v = x);
  M[1, 1] = 1;
  for (j = 2, n, M[, j] = v; if (j < n, v = (A * v) % m));
  matdet(M) % m;
}

\\ (2) Local obstructions. For the modulus m = p^k, returns [a, b] where a
\\ (resp. b) is 1 if I(x) = 1 (resp. -1) is solvable mod m, or -1 if more
\\ than maxpts points would have to be tried (then nothing is proved). By
\\ homogeneity it is enough to try x whose first coordinate prime to p is 1
\\ (those with all coordinates divisible by p have I(x) = 0 mod p): x has a
\\ value s mod m iff I(x) lies in s * (units mod m)^d.
mono_local_signs(nf, T, p, m, maxpts = 10^5) =
{
  my(n = poldegree(nf.pol), d = n * (n - 1) / 2, q = m / p, Dp = Map(), Dm = Map(),
     plus = 0, minus = 0, npts, x, val, t);
  npts = sum(j = 2, n, q^(j - 2) * m^(n - j));
  if (npts > maxpts, return([-1, -1]));
  for (u = 1, m - 1, if (u % p,
    t = lift(Mod(u, m)^d); mapput(Dp, t, 1); mapput(Dm, (-t) % m, 1)));
  for (j = 2, n,
    forvec (w = concat(vector(j - 2, i, [0, q - 1]), vector(n - j, i, [0, m - 1])),
      x = concat([0], concat(concat(p * w[1 .. j - 2], [1]), w[j - 1 .. n - 2]))~;
      val = mono_indexform_mod(T, x, m);
      if (!plus && mapisdefined(Dp, val), plus = 1);
      if (!minus && mapisdefined(Dm, val), minus = 1);
      if (plus && minus, return([1, 1]))));
  [plus, minus];
}

\\ Look for a proof that O_K is not monogenic. Returns 0 if none was found, or
\\ a certificate, one of
\\   ["cid", p]          p is a common index divisor;
\\   ["local", m1, m2]   I(x) = 1 has no solution mod m1, I(x) = -1 none mod m2;
\\ trying the primes p <= maxp. Any one modulus that fails is a proof, so this
\\ never claims a monogenic field is not.
mono_obstruction(nf, maxp = 50, maxpts = 10^5) =
{
  my(n = poldegree(nf.pol), d = n * (n - 1) / 2, cid, T, m, s, mplus = 0, mminus = 0);
  if (n < 3, return(0));
  cid = mono_common_index_divisors(nf);
  if (#cid, return(["cid", cid[1]]));
  T = mono_multables(nf);
  forprime (p = 2, maxp,
    m = mono_local_modulus(p, d);
    if (!m, next);
    s = mono_local_signs(nf, T, p, m, maxpts);
    if (s[1] == 0 && !mplus, mplus = m);
    if (s[2] == 0 && !mminus, mminus = m);
    if (mplus && mminus, return(["local", mplus, mminus])));
  0;
}

\\ (3) Cubic fields: the index form is a binary cubic form I(x2, x3), and PARI's
\\ Thue equation solver finds every solution of I = +-1 (unconditionally, with
\\ thueinit(, 1)), so this decides monogenicity exactly, like Magma's
\\ IndexFormEquation. Returns the canonical generator (as in mono_generator)
\\ as its minimal polynomial, or 0 if O_K is not monogenic.
mono_indexform(nf, v) =
{
  my(n = poldegree(nf.pol), a = nfbasistoalg(nf, v), p = Mod(1, nf.pol), M);
  M = matrix(n, n);
  for (j = 1, n, M[, j] = nfalgtobasis(nf, p); p *= a);
  matdet(M);
}

mono_cubic_generator(nf) =
{
  my(I = (x2, x3) -> mono_indexform(nf, [0, x2, x3]~), a, e, s, t, P, th, sols, best = 0, bestkey, al, key);
  if (poldegree(nf.pol) != 3, error("mono_cubic_generator: degree 3 only"));
  if (mono_index(nf, nfalgtobasis(nf, variable(nf.pol))) == 1, return(nf.pol));
  a = I(1, 0); e = I(0, 1);
  s = I(1, 1) - a - e;            \\ b + c
  t = I(1, -1) - a + e;           \\ c - b
  P = a * 'y^3 + (s - t) / 2 * 'y^2 + (s + t) / 2 * 'y + e;
  th = thueinit(P, 1);
  sols = concat(thue(th, 1), thue(th, -1));
  for (i = 1, #sols,
    al = mono_normalize(nf, nfbasistoalg(nf, [0, sols[i][1], sols[i][2]]~));
    key = mono_key(nf, al);
    if (!best || lex(key, bestkey) < 0, best = al; bestkey = key));
  if (best, minpoly(best), 0);
}

\\ (4) Factoring the index form. I(x) is a polynomial over Z, so if it factors,
\\ I(x) = +-1 forces every factor to be +-1. When some factors are semidefinite
\\ quadratic forms, this bounds x: for a quartic field K with a quadratic subfield
\\ k the cubic resolvent of K is reducible, and I has a factor N(x) that is (up to
\\ sign) the norm form of k on O_K / (O_k + Z), of rank 2 (Gaal, Petho, Pohst).
\\ If k is imaginary (K is then totally complex: every such V4 field, and the D4
\\ fields whose quadratic subfield is imaginary, but no C4 field), N is positive
\\ definite, so N(x) = 1 has finitely many solutions modulo the kernel, and on each
\\ of them I = +-1 is a polynomial equation in the one remaining coordinate, which
\\ nfroots solves exactly. This decides monogenicity and finds all generators.

\\ The index form I(x2, ..., xn) on nf.zk as a polynomial in the variables X.
mono_indexform_sym(nf, X) =
{
  my(n = poldegree(nf.pol), T = mono_multables(nf), A, M = matrix(n, n), v);
  A = sum(i = 2, n, X[i - 1] * T[i]);
  v = vectorv(n, i, i == 1);
  for (j = 1, n, M[, j] = v; v = A * v);
  matdet(M);
}

\\ A unimodular matrix whose last columns span the saturated lattice K (m x k).
mono_complete_basis(K) =
{
  my(m = #K~, k = #K, W, inK, notK);
  W = matsnf(K, 1)[1]^-1;
  inK = [j | j <- [1 .. m], matrank(matconcat([K, W[, j]])) == k];
  notK = [j | j <- [1 .. m], matrank(matconcat([K, W[, j]])) > k];
  if (#inK != k, error("mono_complete_basis: lattice not saturated"));
  matconcat(vector(m, i, W[, if (i <= m - k, notK[i], inK[i - m + k])]));
}

\\ Decide monogenicity by factoring the index form (see above). Returns the
\\ canonical generator (as in mono_generator) as its minimal polynomial, 0 if O_K is
\\ provably not monogenic, or -1 if the method does not apply (no semidefinite
\\ quadratic factors, or they leave more than one coordinate free).
mono_factor_generator(nf) =
{
  my(n = poldegree(nf.pol), X, Y, J, F, g, G, sg, S, Q = List(), r, free, U, IY, Gp, V, w, xw, p,
     sols = List(), best = 0, bestkey, al, key);
  if (n < 3, return(-1));
  if (mono_index(nf, nfalgtobasis(nf, variable(nf.pol))) == 1, return(nf.pol));
  X = vector(n - 1, i, varlower(Str("mono_x", i + 1)));
  Y = vector(n - 1, i, varlower(Str("mono_y", i + 1)));
  J = mono_indexform_sym(nf, X);
  F = factor(J)[, 1];
  S = matrix(n - 1, n - 1);
  for (i = 1, #F,
    g = F[i];
    if (substvec(g, X, 'mono_t * X) != 'mono_t^2 * g, next);
    G = simplify(matrix(n - 1, n - 1, a, b, deriv(deriv(g, X[a]), X[b]) / 2));
    sg = qfsign(G);
    if (sg[1] && sg[2], next);
    listput(Q, [g, if (sg[1], 1, -1)]);
    S += if (sg[1], G, -G));
  if (!#Q, return(-1));
  r = matrank(S);
  free = n - 1 - r;
  if (free > 1, return(-1));
  U = if (free, mono_complete_basis(matkerint(S)), matid(n - 1));
  Gp = (U~ * S * U)[1 .. r, 1 .. r];
  IY = substvec(J, X, Vec(U * Y~));
  V = qfminim(2 * Gp, 2 * #Q, , 2)[3];
  for (k = 1, #V, foreach ([V[, k], -V[, k]], w,
    if (w~ * Gp * w != #Q, next);
    if (!free,
      xw = U * w;
      if (abs(substvec(J, X, Vec(xw))) == 1, listput(sols, xw)),
      p = substvec(IY, Y[1 .. r], Vec(w));
      foreach ([1, -1], e,
        if (p == e, error("mono_factor_generator: index form constant on a line"));
        foreach (nfroots(, subst(p - e, Y[n - 1], 'z)), t,
          if (denominator(t) == 1, listput(sols, U * concat(w, [t]~))))))));
  for (i = 1, #sols,
    al = mono_normalize(nf, nfbasistoalg(nf, concat([0], sols[i])));
    key = mono_key(nf, al);
    if (!best || lex(key, bestkey) < 0, best = al; bestkey = key));
  if (best, minpoly(best), 0);
}
