// blocks.m: the code 480-16-26 built a second time with Magma's group arithmetic, and the Brauer sandwich
// k_ss <= k <= k_C of Proposition 5 of README.md (h0, h2 separately). Optional (needs Magma; about 1 s).
// Run from two-block-and-amplified-codes/:   magma -b src/blocks.m
//   k_ss = 2 sum_i kappa_i(a) kappa_i(b),  k_C = 2 sum_{i,j} c_ij kappa_i(a) kappa_j(b),
// kappa_i(x) = dim ker(x | S_i) on the absolutely irreducible K[G]-modules S_i (K = GF(16) splits G), (c_ij) the
// Cartan matrix; h0 = dim R/(aR + Rb), h2 = dim R/(Ra + bR), R = F_2[G], k = h0 + h2.
SetQuitOnError(true);
SetColumns(0);
F := GF(2);
K := GF(16);

function ReadRows(fn)
    L := [ s : s in Split(Read(fn), "\n") | #s gt 0 ];
    hdr := StringToIntegerSequence(L[1]);
    return hdr[1], [ StringToIntegerSequence(L[1 + i]) : i in [1..hdr[2]] ];
end function;

t0 := Cputime();
Fr<x, y, z> := FreeGroup(3);
Q := quo< Fr | x^5, y^4, y*x*y^-1*x, z^12, (z, x), (z, y) >;
f, G := CosetAction(Q, sub< Q | >);          // regular permutation representation, degree |G|
X := f(Q.1); Y := f(Q.2); Zg := f(Q.3);
D := sub< G | X, Y >;
printf "G: order %o, IdentifyGroup %o; <x,y>: order %o, IdentifyGroup %o; G = <z> x <x,y>: %o\n",
    #G, IdentifyGroup(G), #D, IdentifyGroup(D),
    #(sub< G | Zg > meet D) eq 1 and Zg in Centre(G) and #G eq Order(Zg) * #D;
assert #G eq 240 and IdentifyGroup(G) eq <240, 40> and IdentifyGroup(D) eq <20, 1>;

// element z^c x^i y^j has index 20 c + 4 i + j (data/codes/480-16-26.json)
E := [ Id(G) : t in [1..240] ];
for c in [0..11] do for i in [0..4] do for j in [0..3] do
    E[20*c + 4*i + j + 1] := Zg^c * X^i * Y^j;
end for; end for; end for;
assert #Set(E) eq 240;
ix := func< g | Index(E, g) - 1 >;
a := [ Id(G), Zg^8 * X^3 * Y^2, Zg^11 * Y^3 ];
b := [ Id(G), Zg^2 * X * Y^2, Y ];
HX := [ { ix(s^-1 * c) : s in a } join { 240 + ix(c * t^-1) : t in b } : c in E ];
HZ := [ { ix(c * t) : t in b } join { 240 + ix(s * c) : s in a } : c in E ];
n, PX := ReadRows("instances/c480.hx");
_, PZ := ReadRows("instances/c480.hz");
printf "rows from Magma's multiplication = instances/c480.hx, .hz (row by row): %o, %o\n",
    HX eq [ Set(r) : r in PX ], HZ eq [ Set(r) : r in PZ ];
assert HX eq [ Set(r) : r in PX ] and HZ eq [ Set(r) : r in PZ ];
MX := Matrix(F, 240, 480, [ [ (q in r) select 1 else 0 : q in [0..479] ] : r in HX ]);
MZ := Matrix(F, 240, 480, [ [ (q in r) select 1 else 0 : q in [0..479] ] : r in HZ ]);
printf "rank H_X = %o, rank H_Z = %o, k = %o, H_X H_Z^T = 0: %o\n", Rank(MX), Rank(MZ),
    480 - Rank(MX) - Rank(MZ), IsZero(MX * Transpose(MZ));

// h0, h2 and the sandwich (regular representation: u -> u * g on the points of G)
N := #G;
pos := AssociativeArray();
for i in [1..N] do pos[E[i]] := i; end for;
La := ZeroMatrix(F, N, N); Rb := ZeroMatrix(F, N, N); Ra := ZeroMatrix(F, N, N); Lb := ZeroMatrix(F, N, N);
for u in [1..N] do
    for s in a do La[pos[s*E[u]], u] +:= 1; Ra[pos[E[u]*s], u] +:= 1; end for;
    for t in b do Rb[pos[E[u]*t], u] +:= 1; Lb[pos[t*E[u]], u] +:= 1; end for;
end for;
h0 := N - Rank(HorizontalJoin(La, Rb));
h2 := N - Rank(HorizontalJoin(Ra, Lb));
k := 2*N - Rank(HorizontalJoin(La, Rb)) - Rank(VerticalJoin(Rb, La));
printf "h0 = %o, h2 = %o, k = h0 + h2 = %o\n", h0, h2, k;
assert k eq h0 + h2 and k eq 16;
Ms := [ ChangeRing(M, K) : M in AbsolutelyIrreducibleModules(G, K) ];
r := #Ms;
kap := function(M, S)
    rep := Representation(M);
    return Dimension(M) - Rank(&+[ rep(g) : g in S ]);
end function;
ka := [ kap(M, a) : M in Ms ];
kb := [ kap(M, b) : M in Ms ];
Pims := [ ChangeRing(P, K) : P in ProjectiveIndecomposableModules(G, K) ];
C := ZeroMatrix(Integers(), r, r);
for i in [1..r] do
    hd := [ j : j in [1..r] | IsIsomorphic(ChangeRing(Pims[i]/JacobsonRadical(Pims[i]), K), Ms[j]) ];
    assert #hd eq 1;
    for S in CompositionFactors(Pims[i]) do
        j := [ jj : jj in [1..r] | IsIsomorphic(ChangeRing(S, K), Ms[jj]) ];
        assert #j eq 1;
        C[j[1], hd[1]] +:= 1;
    end for;
end for;
assert C eq Transpose(C);
printf "%o absolutely irreducible modules over GF(16), dimensions %o\n", r, [ Dimension(M) : M in Ms ];
for i in [ i : i in [1..r] | ka[i] gt 0 and kb[i] gt 0 ] do
    printf "simple module %o: dimension %o, kappa(a) = %o, kappa(b) = %o, Cartan entry c_ii = %o, c_ij = 0 for j <> i: %o\n",
        i, Dimension(Ms[i]), ka[i], kb[i], C[i, i], forall{ j : j in [1..r] | j eq i or C[i, j] eq 0 };
end for;
printf "simple modules with kappa(a) > 0: %o; with kappa(b) > 0: %o\n",
    #[ i : i in [1..r] | ka[i] gt 0 ], #[ i : i in [1..r] | kb[i] gt 0 ];
kss := 2 * &+[ ka[i]*kb[i] : i in [1..r] ];
kC := 2 * &+[ C[i, j]*ka[i]*kb[j] : i in [1..r], j in [1..r] ];
printf "k_ss = %o <= h0, h2 <= k_C/2, k_ss <= k = %o <= k_C = %o\n", kss, k, kC;
assert kss le 2*h0 and 2*h0 le kC and kss le 2*h2 and 2*h2 le kC;
printf "time %o s\n", Cputime(t0);
exit;
