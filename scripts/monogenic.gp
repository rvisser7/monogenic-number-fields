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
\\ The precision grows with the coefficients (roots up to about 2^e), so that
\\ generators with huge coefficients get an exact rounded T2 too.
mono_polkey(g) =
{
  my(e = vecmax(apply(c -> if (c, exponent(c), 0), Vec(g))));
  localbitprec(256 + 4 * e);
  concat([round(norml2(polroots(g)) * 10^6)], mono_S(g));
}

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

\\ The T2 form on O_K modulo Z: Q(x) = T2(a) - Tr(a)^2 / n for a = sum_{i>=2} x_i zk[i],
\\ which is the smallest T2(a + s) over all real s, as a (real, positive definite)
\\ Gram matrix in x_2, ..., x_n. (nf.zk[1] = 1, so x_1 only translates a.)
mono_t2mod1(nf) =
{
  my(n = poldegree(nf.pol), G = nf.t2, T);
  T = vector(n, i, nfelttrace(nf, matid(n)[, i]));
  G = G - T~ * T / n;
  G[2..n, 2..n];
}

\\ ---------------------------------------------------------------------------
\\ Fast index test for the searches below. Elements are a = sum_{i>=2} w[i-1] zk[i]
\\ (w a row vector of length n - 1; x_1 only translates a, which does not change
\\ Z[a]), and the index [O_K : Z[a]] = |det(1, a, ..., a^(n-1))| is computed with the
\\ integer multiplication tables of nf.zk only (no nfbasistoalg / nfalgtobasis).
\\
\\ Two quick tests reject almost all elements before that (see mono_isgen). First,
\\ whether a prime p divides the index only depends on w mod p (it says that the
\\ powers of a are dependent mod p), so for p = 2 a table of the 2^(n-1) classes
\\ mod 2 says which classes can contain generators at all, and the search by
\\ height only visits those. a mod 2 must generate O_K / 2 O_K as an F_2-algebra:
\\ usually most classes can (median 70-90% on samples in degrees 8-12), but in
\\ some fields only 15-30% can.

\\ det(1, a, ..., a^(n-1)) on nf.zk, for T = mono_multables(nf) (or T mod p, then the
\\ determinant mod p).
mono_det(T, w) =
{
  my(n = #T, A = w[1] * T[2], M = matrix(n, n), v = vectorv(n, i, i == 1));
  for (i = 2, n - 1, if (w[i], A += w[i] * T[i + 1]));
  for (j = 1, n, M[, j] = v; if (j < n, v = A * v));
  matdet(M);
}

\\ The residue classes mod p (as above) that can contain generators: a Vecsmall t
\\ with t[mono_class(w, p)] = 1 iff p does not divide the index of a.
mono_class(w, p) = fromdigits(w % p, p) + 1;
mono_restable(T, p) =
{
  my(n = #T, t = vectorsmall(p^(n - 1)), Tp = [Mod(1, p) * M | M <- T]);
  forvec (r = vector(n - 1, i, [0, p - 1]),
    if (mono_det(Tp, r) != 0, t[mono_class(r, p)] = 1));
  t;
}

\\ The table mod 2 if it is worth building for about ncand candidates (it costs
\\ 2^(n-1) determinants mod 2, so it is built for at most 2^20 classes), else 0.
mono_maybe_table(T, ncand) =
{
  my(n = #T);
  if (n - 1 > 20 || 2^(n - 1) > ncand / 4, return(0));
  mono_restable(T, 2);
}

\\ The data for mono_isgen: [T, table mod 2 or 0, Z, fq, target]. A second, cheaper
\\ test comes before the determinant: disc(minpoly(a)) = [O_K : Z[a]]^2 disc(K), so a
\\ generator has disc(charpoly(a)) = disc(K). This is checked modulo a prime q of 30
\\ bits, in the power basis where it only takes charpoly and poldisc over F_q (several
\\ times faster than mono_det): with a = Z w~ / d (Z the numerators of nf.zk[2..n] on
\\ the powers of the root, d their denominator), charpoly(d a) has discriminant
\\ d^(n(n-1)) disc(charpoly(a)). Elements passing it (all generators, and a fraction
\\ of about 1/q of the others) are then checked exactly with mono_det.
mono_searchdata(nf) =
{
  my(n = poldegree(nf.pol), T = mono_multables(nf), Z, d, q = 2^30);
  Z = matconcat(vector(n - 1, i, Colrev(lift(nf.zk[i + 1]), n)));
  d = denominator(Z);
  until (d % q && nf.disc % q, q = nextprime(q + 1));
  [T, 0, Mod(d * Z, q), Mod(1, q) * nf.pol, Mod(d, q)^(n * (n - 1)) * nf.disc];
}

\\ Build the table mod 2 in S if it is worth it for about ncand candidates.
mono_add_table(S, ncand) = if (!S[2], S[2] = mono_maybe_table(S[1], ncand)); S;

\\ Whether a (given by w) generates O_K: the table mod 2 (if any), the discriminant
\\ mod q, then the index.
mono_isgen(S, w) =
{
  if (S[2] && !S[2][mono_class(w, 2)], return(0));
  if (poldisc(charpoly(Mod(Polrev(S[3] * w~), S[4]))) != S[5], return(0));
  abs(mono_det(S[1], w)) == 1;
}

\\ The element a of nf for w.
mono_elt(nf, w) = nfbasistoalg(nf, concat([0], w)~);

\\ Prove that the generator alpha of O_K gives the canonical generator, or find the
\\ canonical one. If beta generates O_K and T2(beta) <= T2(alpha), then
\\ Q(beta) <= T2(beta) <= T2(alpha), so beta + k (some k in Z) is among the x with
\\ Q(x) <= T2(alpha) (qfminim, one of each pair +-x); each of those that generates
\\ O_K is normalized as in mono_normalize, which covers every +-beta + k with smallest
\\ T2. So the smallest by mono_key of all of them is the canonical generator.
\\ (mono_key is only computed for generators x with Q(x) at most the T2 of the best
\\ one so far: since T2 >= Q after normalizing, the others cannot beat it.)
\\ Returns [minimal polynomial of the canonical generator, 1], or, when more than
\\ maxpts vectors would have to be checked, [minimal polynomial of the best one
\\ known, 0] (not proved).
mono_canonical(nf, a, maxpts = 10^5) =
{
  my(best = mono_normalize(nf, a), bestkey, t, tb, G = mono_t2mod1(nf), S, V, w, b, key);
  bestkey = mono_key(nf, best);
  t = norml2(polroots(minpoly(best)));
  \\ a margin for rounding: mono_key compares T2 rounded to 10^-6
  V = qfminim(G, t * (1 + 1e-12) + 1e-6, maxpts + 1, 2)[3];
  if (#V > maxpts, return([minpoly(best), 0]));
  S = mono_add_table(mono_searchdata(nf), #V);
  tb = t;
  for (k = 1, #V,
    w = V[, k]~;
    if (!mono_isgen(S, w), next);
    if (qfeval(G, w~) > tb * (1 + 1e-12) + 1e-5, next);
    b = mono_normalize(nf, mono_elt(nf, w));
    key = mono_key(nf, b);
    if (lex(key, bestkey) < 0, best = b; bestkey = key; tb = key[1] / 10^6));
  [minpoly(best), 1];
}

\\ The root of nf.pol, if it generates O_K: then the result of mono_generator.
mono_root_generator(nf, maxpts) =
{
  if (mono_index(nf, nfalgtobasis(nf, variable(nf.pol))) != 1, return(0));
  \\ the root of a polredabs polynomial comes first (README); any other root is checked
  if (polredabs(nf.pol) == nf.pol, [nf.pol, 1],
      mono_canonical(nf, Mod(variable(nf.pol), nf.pol), maxpts));
}

\\ The search for a generator: try alpha = x first; then elements in the order of
\\ T2 (mono_t2_search); if that finds none within maxpts points, all
\\ sum_{i>=2} c_i zk[i] with max |c_i| = 1, 2, ..., B (mono_box_search; B = 0 skips
\\ it). A generator found this way need not be canonical, so it is passed to
\\ mono_canonical, which finds the canonical generator. Returns [minimal
\\ polynomial (t_POL), 1 if proved canonical, else 0], or 0 if no generator is
\\ found: that does NOT prove K is non-monogenic.
\\ With t2 = 0, only the box search is done (the order of earlier versions); the
\\ result is the same whenever it is proved canonical.
mono_generator(nf, B = 3, maxpts = 10^5, t2 = 1) =
{
  my(r, G, S, w = 0);
  if (r = mono_root_generator(nf, maxpts), return(r));
  G = mono_t2mod1(nf);
  S = mono_searchdata(nf);
  if (t2, w = mono_t2_search(S, G, maxpts));
  if (!w && B > 0, w = mono_box_search(S, G, B));
  if (w, mono_canonical(nf, mono_elt(nf, w), maxpts), 0);
}

\\ Search in the order of T2: enumerate the x with Q(x) <= t (Q as in mono_t2mod1)
\\ for growing t, until one of them generates O_K, and return the one with the
\\ smallest Q (as w), or 0 if there are more than maxpts points before that.
\\ The canonical generator has the smallest T2, so it is found as soon as t
\\ reaches about its T2, usually after a few hundred elements, without looking at
\\ the many elements of height 1 with a large T2 (in degree n there are
\\ 3^(n-1) / 2 of those). t starts at twice the minimum of Q and grows by the
\\ factor 2^(2/(n-1)) each round, which about doubles the number of points (they
\\ fill an ellipsoid of dimension n - 1), so at most about twice as many points as
\\ needed are enumerated; the points of the previous rounds are skipped.
mono_t2_search(S, G, maxpts) =
{
  my(n = #G + 1, t, t0 = 0, V, best = 0, bq, q, w);
  t = 2 * qfminim(G, , 1, 2)[2];
  while (1,
    V = qfminim(G, t, maxpts + 1, 2)[3];
    if (#V > maxpts, return(0));
    S = mono_add_table(S, #V);
    for (k = 1, #V,
      w = V[, k]~;
      q = qfeval(G, V[, k]);
      if (q < t0 * (1 - 1e-9), next);   \\ checked in an earlier round
      if (!mono_isgen(S, w), next);
      if (!best || q < bq, best = w; bq = q));
    if (best, return(best));
    t0 = t;
    t *= 2^(2 / (n - 1)));
}

\\ Search all sum_{i>=2} c_i zk[i] with max |c_i| = 1, 2, ..., B, stopping at the
\\ first height where generators exist, and return the one with the smallest Q
\\ (as w), or 0. Only the classes mod 2 that can contain generators are visited
\\ (see above): for each, every coordinate runs over the values in [-h, h] of
\\ its parity.
mono_box_search(S, G, B) =
{
  my(n = #G + 1, tab, cls, best = 0, bq, q);
  S = mono_add_table(S, 3^(n - 1));
  tab = S[2];
  \\ the classes mod 2 to visit, as vectors of 0s and 1s
  cls = if (tab, [s | s <- forvec_list(n - 1), tab[mono_class(s, 2)]], [0]);
  for (h = 1, B,
    foreach (cls, c,
      my(vals = if (c === 0, vector(n - 1, i, [-h .. h]),
                    vector(n - 1, i, [v | v <- [-h .. h], (v - c[i]) % 2 == 0])));
      if ([1 | u <- vals, !#u] != [], next);
      forvec (u = vector(n - 1, i, [1, #vals[i]]),
        my(w = vector(n - 1, i, vals[i][u[i]]), k = 1);
        if (vecmax(abs(w)) != h, next);
        while (k < n && !w[k], k++);
        if (k == n || w[k] < 0, next);
        if (!mono_isgen(S, w), next);
        q = qfeval(G, w~);
        if (!best || q < bq, best = w; bq = q)));
    if (best, return(best)));
  0;
}

\\ All vectors of 0s and 1s of length m.
forvec_list(m) = my(L = List()); forvec (r = vector(m, i, [0, 1]), listput(L, r)); Vec(L);

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

\\ The same from the polynomial f alone, without nfinit for most p: if p^2 does not
\\ divide disc(f), then p does not divide [O_K : Z[a]], and (Dedekind) the residue degrees
\\ of the primes above p are the degrees of the irreducible factors of f mod p. Only the
\\ other p need an order that is maximal at p (nfinit([f, [p]])).
mono_residue_degrees(f, p, d = poldisc(f)) =
{
  if (d % p^2, my(F = factormod(f, p)[, 1]); return(vector(#F, i, poldegree(F[i]))));
  [P.f | P <- idealprimedec(nfinit([f, [p]]), p)];
}

mono_cid_pol(f) =
{
  my(n = poldegree(f), d = poldisc(f), res = List());
  forprime (p = 2, n - 1,
    my(F = vector(n));
    foreach (mono_residue_degrees(f, p, d), k, F[k]++);
    for (k = 1, n, if (F[k] > mono_nirred(p, k), listput(res, p); break)));
  Vec(res);
}

\\ The exact index of the field at p: v_p(i(K)), i(K) = gcd of [O_K : Z[a]] over all a in
\\ O_K with K = Q(a), for nf maximal at p (e.g. nfinit([f, [p]])). Let rad be the product
\\ of the primes P above p, I_j = rad^j, and for an ideal I
\\ L(a, I) = v_p [O_K : Z_p[a] + I]. Then
\\  * L(a, I) only depends on a mod I, and L(a, I) <= v_p [O_K : Z[b]] for every
\\    b = a mod I (since Z_p[b] + I = Z_p[a] + I): a lower bound for the whole class;
\\  * if w = v_p [O_K : Z[a]] then p^w O_K lies in Z_p[a], so L(a, p^k O_K) = w for
\\    k > w, and if I lies in p^(w+1) O_K (I_j does once j >= e (w + 1), e the largest
\\    ramification index) every b = a mod I has Z_p[b] = Z_p[a] (b = a + p c with c in
\\    Z_p[a], and Nakayama), so the same index.
\\ So with U = the smallest index found so far (at first: the root of nf.pol, and 200
\\ random small elements), a depth-first search over the classes mod I_1, I_2, ... that
\\ only refines a class while its bound is < U finds an a of the smallest index: any a with
\\ index w < U has bound < U in every class it lies in, so its class mod I_j with
\\ j = e (w + 1) is reached, and the representative of that class has index w. The
\\ classes mod I_1 are taken up to a -> u a + t (u, t in Z, p not | u), which does not
\\ change the p-part of the index. Returns [v, 1], or [U, 0, lo] with lo <= v <= U if
\\ more than maxcls classes would be needed. This searches O_K as a whole, so the
\\ number of classes grows with the degree; mono_field_index_p below splits it up by the
\\ primes above p (and the tests check that the two agree).
\\ (mono_index_bound: L(a, p^k O_K), with a = w on nf.zk and T = mono_multables(nf).)
mono_index_bound(T, w, p, k) =
{
  my(n = #T, q = p^k, A = sum(i = 2, n, w[i] * T[i]) % q, M = matrix(n, n), v = vectorv(n, i, i == 1));
  for (j = 1, n, M[, j] = v; if (j < n, v = (A * v) % q));
  M = mathnfmodid(M, q);
  valuation(prod(i = 1, n, M[i, i]), p);
}

mono_field_index_global(nf, p, maxcls = 10^6) =
{
  my(n = poldegree(nf.pol), T, U, st = List(), cnt = 0, e, w, j, r, b, dec, emax, emin, rad, H, Q, prec);
  U = (valuation(poldisc(nf.pol), p) - valuation(nf.disc, p)) / 2;
  if (U == 0, return([0, 1]));
  dec = idealprimedec(nf, p);
  emax = vecmax([P.e | P <- dec]); emin = vecmin([P.e | P <- dec]);
  T = mono_multables(nf);
  \\ a quick upper bound first, from some random small elements (the search is exact
  \\ whatever U it starts with, but a small U prunes much more)
  my(seed = getrand());
  setrand(1);
  for (i = 1, 200, my(r0 = mono_index_bound(T, concat([0], vector(n - 1, k, random(5) - 2))~, p, U)); if (r0 < U, U = r0));
  setrand(seed);
  \\ the ideals I_j = rad^j (rad = the product of the P), as HNF matrices, and bases of
  \\ I_j / I_(j+1) (an F_p-space, since p I_j is in I_(j+1)), for j <= emax U
  rad = idealhnf(nf, idealfactorback(nf, dec));
  H = vector(emax * U + 2); H[1] = matid(n);
  for (k = 2, #H, H[k] = idealmul(nf, H[k - 1], rad));
  Q = vector(emax * U + 1, k,
    my(W = Mod(H[k]^-1 * H[k + 1], p), B = if (k == 1, Mat(matid(n)[, 1]), matrix(n, 0)), C);
    for (i = 1, n, C = concat(Mod(B, p), Mod(matid(n)[, i], p));
      if (matrank(concat(W, C)) > matrank(concat(W, Mod(B, p))), B = concat(B, matid(n)[, i])));
    H[k] * B);
  \\ level 1: the classes mod rad, on a basis of O_K / rad that starts with 1, up to
  \\ a -> u a + t (u, t in Z, p not | u): no 1 component, first nonzero coordinate 1 (or
  \\ a = 0 mod rad, which is not the same as a in p O_K if p ramifies)
  my(B1 = Q[1], m = #B1);
  if ((p^(m - 1) - 1) / (p - 1) > maxcls, return([U, 0, 0]));
  forvec (x = vector(m - 1, i, [0, p - 1]),
    my(i = 1);
    while (i < m && !x[i], i++);
    if (i == m || x[i] == 1, listput(st, [B1 * concat([0], x)~, 1, 0])));
  while (#st,
    e = st[#st]; listpop(st);
    if (e[3] >= U, next);
    if (cnt++ > maxcls,
      return([U, 0, min(U, vecmin(vector(#st + 1, i, if (i <= #st, st[i][3], e[3]))))]));
    [w, j] = [e[1], e[2]];
    prec = ceil(j / emin);
    r = mono_index_bound(T, w, p, U);
    if (r < U, U = r);
    \\ the bound for the class w + I_j: v_p [O_K : Z_p[w] + I_j]
    b = if (j < emax * U, mono_ideal_bound(T, w, p, H[j + 1], prec), U);
    if (b < U && j < #Q,
      forvec (x = vector(#Q[j + 1], i, [0, p - 1]), listput(st, [w + Q[j + 1] * x~, j + 1, b]))));
  [U, 1];
}

\\ v_p [O_K : Z_p[a] + I], for an ideal I (HNF matrix H) containing p^k O_K.
mono_ideal_bound(T, w, p, H, k) =
{
  my(n = #T, q = p^k, A = sum(i = 2, n, w[i] * T[i]) % q, M = matrix(n, n), v = vectorv(n, i, i == 1));
  for (j = 1, n, M[, j] = v; if (j < n, v = (A * v) % q));
  M = mathnfmodid(concat(M, H), q);
  valuation(prod(i = 1, n, M[i, i]), p);
}

\\ v_p(i(K)) as above, much faster: O_K (x) Z_p is the product of the O_P for the primes
\\ P above p, and for a in O_K let phi_P be the minimal polynomial over F_p of a mod P.
\\ If G is the set of P with the same phi, a_G the component of a in O_G = prod_{P in G}
\\ O_P, then the charpolys of the a_G are coprime mod p, so Z_p[a] = prod_G Z_p[a_G] and
\\ v_p [O_K : Z[a]] = sum_G v_p [O_G : Z_p[a_G]]. So
\\   v_p(i(K)) = min over the ways to give each P an irreducible phi_P with deg phi_P | f_P
\\               of sum over the groups G (same phi) of c(G, phi),
\\ c(G, phi) = the smallest v_p [O_G : Z_p[b]] over b in O_G with phi(b) in every P of G,
\\ found by the search above inside O_G (mono_group_index). The minimum over the ways is a
\\ dynamic program over the subsets of the P, one phi at a time. Simplifications, all exact:
\\  * if every P in G is unramified with f_P = d = deg phi (k of them), then
\\    c(G, phi) = d (k(k-1)/2 + g(k)), where g(k) is the smallest sum of v(x_i - x_j) over
\\    k distinct x_i in W = W(F_q), q = p^d (mono_split_index): the index only depends on
\\    the charpolys of the b_P, so (applying Frobenius to each O_P = W) all b_P = w + p y_P
\\    with w the Teichmuller lift of one root of phi, and then only the d conjugate pairs
\\    of roots that agree mod p count, each with v = 1 + v(y_P - y_Q). With k = 1: c = 0;
\\  * c(G, phi) >= the number of pairs in G (each pair P, Q has v_p(Res) >= 1) plus the
\\    number of P in G with deg phi < f_P (then p | [O_P : Z_p[b]]); c(G, phi) is not
\\    computed when this bound cannot improve the minimum;
\\  * c(G, phi) = c(G, psi) for psi(x) = u^d phi((x - t) / u), u, t in Z, p not | u (b -> u b + t);
\\  * a P such that there are more irreducible polynomials of degree f_P than other P that
\\    could use one (d | f_Q) can always get a phi of its own, at cost 0, and removing P
\\    from a group does not increase its cost, so such P are left out (repeatedly).
\\ Returns [v, 1], or [U, 0, lo] (lo <= v <= U) if a search needs more than maxcls classes
\\ or there are too many ways (more than maxcls subset pairs).
\\ g(k): the smallest sum over pairs of v(x_i - x_j) for k distinct x_i in W(F_q): split
\\ them by residue into q classes (sizes s_i), each pair in a class has v >= 1, so
\\ g(k) = min over the sizes of sum (s_i (s_i - 1) / 2 + g(s_i)) (s_i < k).
mono_split_index(q, k, memo = Map()) =
{
  my(r, h);
  if (k <= q, return(0));
  if (mapisdefined(memo, k, &r), return(r));
  \\ h[j + 1][s + 1] = best way to put s elements into j classes
  my(J = min(q, k));
  h = vector(J + 1, j, vector(k + 1, s, if (j == 1, if (s == 1, 0, oo), oo)));
  for (j = 2, J + 1,
    for (s = 0, k,
      for (t = 0, min(s, k - 1),
        my(prev = h[j - 1][s - t + 1]);
        if (prev < oo, h[j][s + 1] = min(h[j][s + 1], prev + t * (t - 1) / 2 + mono_split_index(q, t, memo))))));
  r = h[J + 1][k + 1];
  mapput(memo, k, r);
  r;
}

mono_irreducibles(p, d) =
{
  my(L = List());
  forvec (c = vector(d, i, [0, p - 1]),
    my(g = 'x^d + Pol(c, 'x));
    if (polisirreducible(Mod(1, p) * g), listput(L, g)));
  Vec(L);
}

\\ The representative of phi up to x -> u x + t (mod p), as a coefficient vector.
mono_phi_key(phi, p) =
{
  my(d = poldegree(phi), best);
  for (u = 1, p - 1,
    for (t = 0, p - 1,
      my(g = lift(Mod(u, p)^d * subst(Mod(1, p) * phi, 'x, ('x - t) / Mod(u, p))), k = Vec(g));
      if (!best || lex(k, best) < 0, best = k)));
  best;
}

\\ v_p [O_K : span(e a^j, j < nG) + (1 - e) O_K + p^k O_K] = v_p [O_G : Z_p[a_G] + p^k O_G].
mono_group_bound(T, e, C, nG, a, p, k) =
{
  my(n = #T, q = p^k, A = sum(i = 1, n, a[i] * T[i]) % q, M = matrix(n, nG), v = e % q, H);
  for (j = 1, nG, M[, j] = v; if (j < nG, v = (A * v) % q));
  H = mathnfmodid(concat(M, C % q), q);
  valuation(prod(i = 1, n, H[i, i]), p);
}

\\ c(G, phi) as above, or cap if it is >= cap; [c, number of classes, done].
mono_group_index(nf, T, dec, E, G, phi, p, cap, maxcls) =
{
  my(n = #T, e = vectorv(n), nG = 0, C, B, r, st = List(), cnt = 0, U = cap, ent, a, k, b, one = vectorv(n, i, i == 1));
  for (i = 1, #G, e += E[G[i]]; nG += dec[G[i]].e * dec[G[i]].f);
  C = sum(i = 1, n, (one - e)[i] * T[i]);
  \\ a basis of O_G mod p starting with e (translations by Z are left out)
  B = Mat(e); r = 1;
  for (j = 1, n, if (r == nG, break);
    my(B2 = concat(B, (sum(i = 1, n, e[i] * T[i]))[, j]));
    if (matrank(Mod(B2, p)) > r, B = B2; r++));
  if (r < nG, error("mono_group_index: no basis of O_G"));
  B = B[, 2 .. nG];
  if (p^nG > maxcls, return([cap, 0, 0]));
  \\ the classes mod p: b mod p O_G with phi(b) in every P of G
  forvec (c = vector(nG, i, [0, p - 1]),
    my(a0 = B * c[2 .. nG]~ + c[1] * e, y = vectorv(n), v = one);
    for (i = 0, poldegree(phi), y += polcoef(phi, i) * v; v = nfeltmul(nf, v, a0) % p);
    if (sum(i = 1, #G, nfeltval(nf, y % p, dec[G[i]]) > 0) == #G, listput(st, [a0, 1, 0])));
  while (#st,
    ent = st[#st]; listpop(st);
    if (ent[3] >= U, next);
    if (cnt++ > maxcls, return([U, cnt, 0]));
    [a, k] = [ent[1], ent[2]];
    r = mono_group_bound(T, e, C, nG, a, p, U);
    if (r < U, U = r);
    b = if (k < U, mono_group_bound(T, e, C, nG, a, p, k), U);
    if (b < U,
      forvec (c = vector(nG - 1, i, [0, p - 1]), listput(st, [a + p^k * B * c~, k + 1, b]))));
  [U, cnt, 1];
}

mono_field_index_p(nf, p, maxcls = 10^6) =
{
  my(n = poldegree(nf.pol), U0, dec, m, keep, F, T, K, E, fa, phis, dp, memo = Map(), done = 1, work = 0, lo);
  U0 = (valuation(poldisc(nf.pol), p) - valuation(nf.disc, p)) / 2;
  if (U0 == 0, return([0, 1]));
  dec = idealprimedec(nf, p); m = #dec;
  F = vector(m, i, dec[i].f);
  \\ leave out the P that can always get a polynomial of their own
  keep = vector(m, i, 1);
  my(avail = Map(), changed = 1);
  while (changed,
    changed = 0;
    for (i = 1, m, if (keep[i],
      my(demand = sum(j = 1, m, keep[j] && j != i && F[j] % F[i] == 0), have);
      if (!mapisdefined(avail, F[i], &have), have = mono_nirred(p, F[i]));
      if (have > demand, keep[i] = 0; mapput(avail, F[i], have - 1); changed = 1))));
  keep = [i | i <- [1 .. m], keep[i]];
  if (!#keep, return([0, 1]));
  T = mono_multables(nf);
  K = U0 + 1;
  \\ idempotents mod p^K
  fa = matconcat([Col(dec), Col(vector(m, i, dec[i].e * K))]);
  E = vector(m, i, nfalgtobasis(nf, idealchinese(nf, fa, vector(m, j, j == i))));
  \\ the polynomials: every degree dividing an f_P of a kept P
  phis = List();
  foreach (Set(concat(vector(#keep, i, divisors(F[keep[i]])))), d,
    foreach (mono_irreducibles(p, d), g, listput(phis, g)));
  \\ the kept P by type: the unramified P with the same f are interchangeable (all O_P are
  \\ W(F_(p^f))), so the dynamic program only counts how many of each type a group has
  my(tys = List(), INF = K + 1, R, NS, sv);
  for (i = 1, #keep, my(P = dec[keep[i]], t = 0);
    if (P.e == 1, for (j = 1, #tys, if (tys[j][1] == 1 && tys[j][2] == P.f, t = j; break)));
    if (t, listput(tys[t][3], keep[i]), listput(tys, [P.e, P.f, List([keep[i]])])));
  my(nt = #tys, cnt = vector(nt, t, #tys[t][3]));
  R = vector(nt); R[1] = 1; for (t = 2, nt, R[t] = R[t - 1] * (cnt[t - 1] + 1));
  NS = R[nt] * (cnt[nt] + 1);
  if (NS > maxcls, return([U0, 0, 1]));
  sv = vector(NS, i, my(x = i - 1); vector(nt, t, (x \ R[t]) % (cnt[t] + 1)));  \\ state i -> counts
  dp = vector(NS, i, INF); dp[1] = 0;
  foreach (phis, phi,
    my(d = poldegree(phi), nd = dp, key = mono_phi_key(phi, p), comp = vector(nt, t, tys[t][2] % d == 0));
    if (!vecmax(comp), next);
    for (S = 2, NS,
      my(s = sv[S]);
      forvec (tau = vector(nt, t, [0, if (comp[t], s[t], 0)]),
        my(k = vecsum(tau), rest, base, lb, c);
        if (!k, next);
        rest = S - sum(t = 1, nt, tau[t] * R[t]);
        base = dp[rest];
        if (base >= INF, next);
        lb = k * (k - 1) / 2 + sum(t = 1, nt, tau[t] * (tys[t][2] > d));
        if (base + lb >= nd[S], next);
        if (vecmin(vector(nt, t, !tau[t] || (tys[t][1] == 1 && tys[t][2] == d))),
          c = d * (k * (k - 1) / 2 + mono_split_index(p^d, k)),
          my(mkey = [tau, key], got, cap = min(INF, nd[S] - base));
          if (mapisdefined(memo, mkey, &got), c = got,
            my(G = concat(vector(nt, t, vector(tau[t], i, tys[t][3][i]))), res);
            res = mono_group_index(nf, T, dec, E, G, phi, p, cap, maxcls);
            work += res[2];
            if (!res[3], done = 0);
            c = res[1];
            \\ only an exact value (below the cap) can be reused
            if (res[3] && c < cap, mapput(memo, mkey, c))));
        nd[S] = min(nd[S], base + c);
        if (work > maxcls, return([U0, 0, 1]))));
    dp = nd);
  my(S0 = NS - 1);
  if (dp[S0 + 1] > U0, error("mono_field_index_p: found no a of index <= ", U0));
  if (!done, return([dp[S0 + 1], 0, 1]));
  [dp[S0 + 1], 1];
}

\\ v_p(i(K)) for nf maximal at p: mono_field_index_global (fast in low degree, and
\\ when p ramifies) and mono_field_index_p (fast when p is unramified: the closed form),
\\ in that order or the other way round, the second one only if the first gives up
\\ (more than maxcls classes). Returns [v, 1], or [U, 0, lo] with lo <= v <= U if both
\\ give up.
mono_field_index(nf, p, maxcls = 2 * 10^4) =
{
  my(unram = vecmax([P.e | P <- idealprimedec(nf, p)]) == 1, r, s);
  r = if (unram, mono_field_index_p(nf, p, maxcls), mono_field_index_global(nf, p, maxcls));
  if (r[2], return(r));
  s = if (unram, mono_field_index_global(nf, p, maxcls), mono_field_index_p(nf, p, maxcls));
  if (s[2], return(s));
  [min(r[1], s[1]), 0, max(if (#r > 2, r[3], 0), if (#s > 2, s[3], 0))];
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

\\ All x in Z^3 with U(x) = u, V(x) = v (U, V rational ternary forms in X, (u, v)
\\ != (0, 0)), or 0 if this fails (a degenerate case). Such x lie on the conic
\\ H = v U - u V = 0, so (qfsolve, qfparam) x = (lam / g) w(p, q) with
\\ w = M [p^2, pq, q^2]~, p, q coprime, g the content of w(p, q) and lam an
\\ integer. As adj(M) w = det(M) [p^2, pq, q^2]~, g divides det(M). Then u = U(x) (or
\\ v = V(x) if u = 0) becomes the quartic Thue equation
\\ U(w(p, q)) = u g^2 / lam^2, for the finitely many possible g and lam, which
\\ PARI's thue solves completely (thueinit(, 1)). This is how Gaal, Petho and
\\ Pohst solve index form equations in quartic fields.
mono_conic_solve(U, V, X, u, v) =
{
  my(T = if (u, [U, u], [V, v]), H = v * U - u * V, G, s0, M, w, F, dF, k = 0, Fk, tF, Del, out = List());
  G = simplify(matrix(3, 3, a, b, deriv(deriv(H, X[a]), X[b]) / 2));
  if (G == 0 || matdet(G) == 0, return(0));
  G /= content(G);
  s0 = qfsolve(G);
  if (type(s0) != "t_COL", return([]));    \\ no rational point on the conic
  M = qfparam(G, s0);
  Del = matdet(M);
  if (!Del, return(0));
  w = M * ['p^2, 'p * 'q, 'q^2]~;
  F = substvec(T[1], X, Vec(w));
  if (F == 0, return(0));
  dF = denominator(mono_rcontent(F)); F *= dF;
  while (substvec(F, ['p, 'q], [1, k]) == 0, k++);
  Fk = substvec(F, ['p, 'q], ['y, 1 + k * 'y]);   \\ q = q' + k p, dehomogenized at q' = 1
  if (poldegree(Fk) != 4 || !issquarefree(Fk), return(0));
  tF = thueinit(Fk, 1);
  fordiv (abs(Del), g,
    my(t = T[2] * g^2 * dF);
    if (denominator(t) != 1 || !t, next);
    fordiv (abs(t), l2,
      if (!issquare(l2), next);
      my(lam = sqrtint(l2));
      foreach (thue(tF, t / l2), pq,
        my(pp = pq[1], qq = pq[2] + k * pq[1], xx);
        if (gcd(pp, qq) != 1, next);
        xx = substvec(w, ['p, 'q], [pp, qq]) * lam / g;
        foreach ([xx, -xx], y,
          if (denominator(y) == 1 && substvec(U, X, Vec(y)) == u && substvec(V, X, Vec(y)) == v,
            listput(out, y))))));
  Set(out);
}

\\ Decide a quartic field with the cubic resolvent (see (5) above). The pairs
\\ (u, v) = (U(x), V(x)) allowed by I(x) = +-1 are: with a rational root of P, the
\\ points where a line meets a conic; if P is irreducible, the solutions of the
\\ Thue equation c v^3 P(u / v) = +-1 (PARI's thue, unconditional with
\\ thueinit(, 1)), scaled to integers. If K is totally complex, with a1, a2 and a3,
\\ a4 complex conjugate, phi_1 - phi_2 = |a1 - conj(a3)|^2 and phi_1 - phi_3 =
\\ |a1 - a3|^2 are >= 0, so some U - s V is positive or negative definite; on it,
\\ U(x) = u, V(x) = v leaves finitely many x (qfminim), so this finds every
\\ generator. Otherwise each (u, v) is first tested modulo small prime powers,
\\ and the ones left are solved exactly with mono_conic_solve. So this decides
\\ every quartic field with irreducible nf.pol, up to the degenerate cases where
\\ it returns -1. Returns the canonical generator, 0 if O_K is provably not
\\ monogenic, -1 if undecided.
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
  my(GU = simplify(matrix(3, 3, a, b, deriv(deriv(U, X[a]), X[b]) / 2)),
     GV = simplify(matrix(3, 3, a, b, deriv(deriv(V, X[a]), X[b]) / 2)));
  pen = if (nf.sign[1] == 0, mono_definite_pencil(GU, GV), 0);
  if (pen,
    my(G = pen[1], d = denominator(G));
    foreach (cands, uv,
      my(val = pen[3] * if (pen[2] == oo, uv[2], uv[1] - pen[2] * uv[2]), W);
      if (val <= 0, next);
      W = qfminim(2 * d * G, floor(2 * d * val), , 0)[3];
      for (k = 1, #W, foreach ([W[, k], -W[, k]], x,
        if (substvec(U, X, Vec(x)) == uv[1] && substvec(V, X, Vec(x)) == uv[2], listput(sols, x))))),
    foreach (cands, uv,
      if ([1 | m <- moduli, !mono_local_pair(U, V, X, uv[1], uv[2], m)] != [], next);
      my(S = mono_conic_solve(U, V, X, uv[1], uv[2]));
      if (S === 0, return(-1));
      foreach (S, x, listput(sols, x))));
  sols = [x | x <- sols, mono_index(nf, concat([0], x)) == 1];
  for (i = 1, #sols,
    al = mono_normalize(nf, nfbasistoalg(nf, concat([0], sols[i])));
    key = mono_key(nf, al);
    if (!best || lex(key, bestkey) < 0, best = al; bestkey = key));
  if (best, minpoly(best), 0);
}
