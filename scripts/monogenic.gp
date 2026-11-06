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
\\ provably not monogenic, or -1 if the method does not decide (no semidefinite
\\ quadratic factors, or they leave more than one coordinate free). For quartic
\\ fields it then tries mono_resolvent_generator (below): exact for every totally
\\ complex quartic field, and for the others it can prove non-monogenicity.
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
  if (!#Q, return(if (n == 4, mono_resolvent_generator(nf, X, J), -1)));
  r = matrank(S);
  free = n - 1 - r;
  if (free > 1, return(if (n == 4, mono_resolvent_generator(nf, X, J), -1)));
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

\\ (5) Quartic fields with a quadratic subfield, in general (also real ones). Write
\\ alpha = sum x_i zk[i] with conjugates a1..a4, and phi_1 = a1 a2 + a3 a4,
\\ phi_2 = a1 a3 + a2 a4, phi_3 = a1 a4 + a2 a3. Each phi_k is a quadratic form in
\\ x = (x2, x3, x4), phi_k = A0 + A1 t_k + A2 t_k^2, where t_k are the roots of the
\\ cubic resolvent P of nf.pol and A0, A1, A2 are forms with rational coefficients,
\\ so phi_i - phi_j = (t_i - t_j) (U - t_l V) with U = A1 + (t_1 + t_2 + t_3) A2,
\\ V = A2, and I(x) = c * V^3 P(U/V) = c * prod_l (U - t_l V) (Gaal, Petho and
\\ Pohst). U and V are found numerically, then the identity is checked exactly, so
\\ nothing below depends on the numerics.
\\
\\ If K has a quadratic subfield, P has a rational root t, so I = c * L * M with
\\ L = U - t V and M = V^2 (P / (T - t))(U / V). As forms in x, L = lam * N and
\\ M = mu * Q with N and Q primitive over Z, so I(x) = +-1 forces N(x) = +-1 and
\\ Q(x) = +-1, i.e. L(U, V) = +-lam and M(U, V) = +-mu: a line and a conic in the
\\ (U, V) plane, which meet in at most two points per choice of signs. These must
\\ be rational, with U in (1/dU) Z and V in (1/dV) Z (dU, dV the denominators of U
\\ and V). E.g. for 4.2.1600.1 there are none, so it is not monogenic. Every
\\ rational root of P gives such conditions; the candidates must satisfy all.
\\ Candidates (u, v) left over are then tested locally: U(x) = u, V(x) = v must be
\\ solvable modulo small prime powers.
\\
\\ Returns 0 if this proves O_K is not monogenic, and -1 otherwise (it never finds
\\ generators: the equations U(x) = u, V(x) = v can have infinitely many solutions
\\ in the real case; leave those fields to the search).

\\ The rational content of a polynomial in several variables.
mono_rcontent(p) = if (type(p) == "t_POL", content(apply(mono_rcontent, Vec(p))), p);

\\ [U, V, P, c] as above, for alpha on nf.zk with x1 = 0 and the variables X, or 0
\\ if the exact check fails. J is the index form in X.
mono_resolvent_UV(nf, X, J) =
{
  my(f = nf.pol, v = variable(f), r, B, th, Vi, A1 = 0, A2 = 0, s, P, U, V, c,
     pairs = [[1, 2, 3, 4], [1, 3, 2, 4], [1, 4, 2, 3]], ap = (z) -> bestappr(real(z), 10^40));
  if (poldegree(f) != 4, return(0));
  localprec(200);
  r = polroots(f);
  B = vector(4, j, vector(4, i, subst(lift(nf.zk[j]), v, r[i])));
  th = vector(3, k, my(p = pairs[k]); r[p[1]] * r[p[2]] + r[p[3]] * r[p[4]]);
  Vi = matrix(3, 3, k, e, th[k]^(e - 1))^-1;
  for (j = 2, 4, for (l = j, 4,
    my(cv, rr);
    cv = vector(3, k, my(p = pairs[k]);
      sum(h = 0, 1, my(a = p[2 * h + 1], b = p[2 * h + 2]);
        if (j == l, B[j][a] * B[j][b], B[j][a] * B[l][b] + B[l][a] * B[j][b])));
    rr = Vi * cv~;
    A1 += ap(rr[2]) * X[j - 1] * X[l - 1];
    A2 += ap(rr[3]) * X[j - 1] * X[l - 1]));
  s = ap(vecsum(th));
  P = 'T^3 - s * 'T^2 + ap(th[1] * th[2] + th[1] * th[3] + th[2] * th[3]) * 'T - ap(th[1] * th[2] * th[3]);
  U = A1 + s * A2; V = A2;
  if (V == 0, return(0));
  c = simplify(J / (V^3 * subst(P, 'T, U / V)));
  if (type(c) != "t_INT" && type(c) != "t_FRAC", return(0));
  [U, V, P, c];
}

\\ Whether q1(x) = u and q2(x) = v (rational ternary forms in X) is solvable mod m,
\\ as a necessary condition (both sides scaled to integers).
mono_local_pair(q1, q2, X, u, v, m) =
{
  my(G = vector(2), t = vector(2), d, M = 2 * m);
  foreach ([[1, q1, u], [2, q2, v]], e,
    d = denominator(mono_rcontent(e[2])) * denominator(e[3]);
    G[e[1]] = simplify(matrix(3, 3, a, b, deriv(deriv(2 * d * e[2], X[a]), X[b])));
    t[e[1]] = 2 * d * e[3]);
  forvec (w = vector(3, i, [0, m - 1]),
    my(x = w~);
    if ((x~ * G[1] * x / 2 - t[1]) % M == 0 && (x~ * G[2] * x / 2 - t[2]) % M == 0, return(1)));
  0;
}

\\ A positive definite member G = sign * (U - s V) (or sign * V) of the pencil of
\\ quadratic forms with Gram matrices GU, GV, as [G, s, sign] (s = oo for V), or 0.
mono_definite_pencil(GU, GV) =
{
  my(rr = vecsort([real(t) | t <- polroots(matdet(GU - 'w * GV)), abs(imag(t)) < 1e-10]), ts);
  ts = if (#rr, concat([rr[1] - 1, rr[#rr] + 1], vector(#rr - 1, i, (rr[i] + rr[i + 1]) / 2)), [0]);
  foreach (concat([oo], ts), t,
    my(s = if (t == oo, oo, bestappr(t, 10^6)), G = if (s == oo, GV, GU - s * GV), sg = qfsign(G));
    if (sg[1] == 3, return([G, s, 1]));
    if (sg[2] == 3, return([-G, s, -1])));
  0;
}

\\ Decide a quartic field with the cubic resolvent (see (5) above). The pairs
\\ (u, v) = (U(x), V(x)) allowed by I(x) = +-1 are: with a rational root of P, the
\\ points where a line meets a conic; if P is irreducible, the solutions of the
\\ Thue equation c v^3 P(u / v) = +-1 (PARI's thue, unconditional with
\\ thueinit(, 1)), scaled to integers. If K is totally complex, with a1, a2 and a3,
\\ a4 complex conjugate, phi_1 - phi_2 = |a1 - conj(a3)|^2 and phi_1 - phi_3 =
\\ |a1 - a3|^2 are >= 0, so some U - s V is positive or negative definite; on it,
\\ U(x) = u, V(x) = v leaves finitely many x (qfminim), so this finds every
\\ generator. Otherwise each (u, v) is tested modulo small prime powers.
\\ Returns the canonical generator, 0 if O_K is provably not monogenic, -1 if
\\ undecided.
mono_resolvent_generator(nf, X, J, moduli = [2, 4, 8, 16, 32, 3, 9, 27, 5, 25, 7, 11, 13]) =
{
  my(R, U, V, P, c, rts, dU, dV, cands = 0, pen, sols = List(), best = 0, bestkey, al, key);
  if (#mono_common_index_divisors(nf), return(0));
  R = mono_resolvent_UV(nf, X, J);
  if (!R, return(-1));
  [U, V, P, c] = R;
  if (denominator(content(P)) != 1, return(-1));
  dU = denominator(mono_rcontent(U)); dV = denominator(mono_rcontent(V));
  rts = nfroots(, subst(P, 'T, 'z));
  if (#rts,
    foreach (rts, t,
      my(L = U - t * V, Mh = P / ('T - t), lam, mu, new = List());
      lam = mono_rcontent(L);
      mu = mono_rcontent(simplify(V^2 * subst(Mh, 'T, U / V)));
      foreach ([lam, -lam], l0, foreach ([mu, -mu], m0,
        my(q = simplify('w^2 * subst(Mh, 'T, (l0 + t * 'w) / 'w)) - m0);
        if (q == 0, return(-1));          \\ the line lies on the conic: no bound
        foreach (nfroots(, q / mono_rcontent(q)), v0,
          my(u0 = l0 + t * v0);
          if (denominator(u0 * dU) == 1 && denominator(v0 * dV) == 1, listput(new, [u0, v0])))));
      cands = if (cands === 0, Set(new), setintersect(cands, Set(new)))),
    my(D = lcm(dU, dV), rhs = D^3 / c, th, ts);
    if (denominator(rhs) != 1, return(0));
    th = thueinit(subst(P, 'T, 'y), 1);
    ts = concat(thue(th, rhs), thue(th, -rhs));
    cands = Set([[s[1] / D, s[2] / D] | s <- ts, denominator(s[1] / D * dU) == 1 && denominator(s[2] / D * dV) == 1]));
  if (!#cands, return(0));
  if (nf.sign[1] != 0,
    foreach (cands, uv,
      if ([1 | m <- moduli, !mono_local_pair(U, V, X, uv[1], uv[2], m)] == [], return(-1)));
    return(0));
  my(GU = simplify(matrix(3, 3, a, b, deriv(deriv(U, X[a]), X[b]) / 2)),
     GV = simplify(matrix(3, 3, a, b, deriv(deriv(V, X[a]), X[b]) / 2)));
  pen = mono_definite_pencil(GU, GV);
  if (!pen, return(-1));
  my(G = pen[1], d = denominator(G));
  foreach (cands, uv,
    my(val = pen[3] * if (pen[2] == oo, uv[2], uv[1] - pen[2] * uv[2]), W);
    if (val <= 0, next);
    W = qfminim(2 * d * G, floor(2 * d * val), , 0)[3];
    for (k = 1, #W, foreach ([W[, k], -W[, k]], x,
      if (substvec(U, X, Vec(x)) == uv[1] && substvec(V, X, Vec(x)) == uv[2]
          && mono_index(nf, concat([0], x)) == 1, listput(sols, x)))));
  for (i = 1, #sols,
    al = mono_normalize(nf, nfbasistoalg(nf, concat([0], sols[i])));
    key = mono_key(nf, al);
    if (!best || lex(key, bestkey) < 0, best = al; bestkey = key));
  if (best, minpoly(best), 0);
}
