// Decide monogenicity with Magma's index form equation solver, and find the
// canonical generator alpha of O_K when there is one.
//
//   magma -b infile:=fields.txt outfile:=rows.txt magma/monogenic.m
//
// infile has lines "label|[a0, a1, ..., an]" (ascending polredabs coefficients).
// For each field where the root of the polynomial generates O_K, or where
// IndexFormEquation applies (degree <= 4; in higher degree Magma needs a known
// power basis already), a row is appended to outfile:
//   label|1|[minimal polynomial of alpha, ascending]   if O_K is monogenic
//   label|0|\N                                         if it provably is not
// Other fields are listed on stdout for the brute-force search instead.
//
// The choice of alpha matches scripts/monogenic.gp: the root of the defining
// polynomial if it generates O_K, otherwise the generator of smallest T2 norm
// (normalised among +-alpha + k), ties broken by the minimal polynomial.

ZZ := Integers();
R<x> := PolynomialRing(ZZ);

function LexLess(u, v)
    for i in [1..Min(#u, #v)] do
        if u[i] ne v[i] then
            return u[i] lt v[i];
        end if;
    end for;
    return #u lt #v;
end function;

// Descending coefficients of the minimal polynomial, like PARI's Vec(minpoly(a)).
function MinpolyVec(a)
    return Reverse([ZZ | c : c in Coefficients(MinimalPolynomial(a))]);
end function;

// The smallest of +-a + k: smallest trace, then smallest minimal polynomial.
function Normalize(a, n)
    a := a - Floor(ZZ!Trace(a) / n);
    t := ZZ!Trace(a);
    b := t ne 0 select 1 - a else -a;
    tb := ZZ!Trace(b);
    if tb lt t or (tb eq t and LexLess(MinpolyVec(b), MinpolyVec(a))) then
        return b;
    end if;
    return a;
end function;

// Sort key [T2 rounded to 6 decimals] cat minimal polynomial, as mono_key in PARI.
function Key(a, K)
    t2 := &+[Abs(c)^2 : c in Conjugates(K!a)];
    return [Round(10^6 * t2)] cat MinpolyVec(a);
end function;

function CoeffString(f)
    return "[" cat Join([IntegerToString(ZZ!c) : c in Coefficients(f)], ", ") cat "]";
end function;

// Returns the row to write, or "" when IndexFormEquation does not apply.
function ProcessField(label, coeffs)
    f := R!coeffs;
    n := Degree(f);
    O := MaximalOrder(f);
    K := NumberField(O);
    if Discriminant(f) eq Discriminant(O) then
        return label cat "|1|" cat CoeffString(f);
    end if;
    if n gt 4 then
        return "";
    end if;
    sols := IndexFormEquation(O, 1);
    if #sols eq 0 then
        return label cat "|0|\\N";
    end if;
    best := Normalize(sols[1], n);
    bestkey := Key(best, K);
    for s in sols[2..#sols] do
        a := Normalize(s, n);
        k := Key(a, K);
        if LexLess(k, bestkey) then
            best := a;
            bestkey := k;
        end if;
    end for;
    g := R!MinimalPolynomial(best);
    assert Discriminant(g) eq Discriminant(O);
    return label cat "|1|" cat CoeffString(g);
end function;

if not assigned infile or not assigned outfile then
    print "usage: magma -b infile:=fields.txt outfile:=rows.txt magma/monogenic.m";
    quit;
end if;

F := Open(infile, "r");
while true do
    line := Gets(F);
    if IsEof(line) then
        break;
    end if;
    if #line eq 0 or line[1] eq "#" then
        continue;
    end if;
    parts := Split(line, "|");
    label := parts[1];
    coeffs := [ZZ | c : c in eval parts[2]];
    row := ProcessField(label, coeffs);
    if row eq "" then
        print "not handled:", label;
    else
        Write(outfile, row);
    end if;
end while;
quit;
