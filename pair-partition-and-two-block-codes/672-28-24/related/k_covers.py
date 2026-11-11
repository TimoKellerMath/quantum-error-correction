#!/usr/bin/env python3
"""k of two-block group algebra codes related to [[672,28,24]] (README.md, Part B, Lemma 8).

    python3 672-28-24/related/k_covers.py    (prints k_covers.out; a few seconds, standard library only)

1. The weight-6 pair of Qian and Li (arXiv:2608.08996) over G' = Z84 x|_29 Z4 = <x, s | x^84 = s^4 = 1,
   s x s^-1 = x^29>: A = e + x s + x^3, B = s^3 + x^2 + x^4, over G'/K with K = <x^42> (n = 336, their
   [[336,12,20]]) and over G' itself (K trivial, n = 672), with A acting from the left (H_X = [L(A) | R(B)]) and
   from the right (H_X = [R(A) | L(B)]).
2. [[672,28,24]] = C(G; a, b) with G = <x, y | x^12 = y^28 = 1, y x y^-1 = x^5>, a = 1 + x^2 y^9 + x^10 y^11,
   b = 1 + x^10 y^22 + x^7 y^27 and H_X = [L(a) | R(b)] (the convention of the board entry), and its quotients
   C(G/<z>; a, b) by the central involutions z of G.
The code over G/K (K central) has one X-check and one Z-check per coset cK; coefficients are reduced mod 2.
k = n - rank H_X - rank H_Z over GF(2), after checking H_X H_Z^T = 0.
"""


def semidirect(m, q, r):
    """Z_m x|_r Z_q: elements (i, j) = x^i y^j, (i, j)(k, l) = (i + r^j k, j + l)."""
    assert pow(r, q, m) == 1 % m
    elts = [(i, j) for i in range(m) for j in range(q)]

    def mul(g, h):
        return ((g[0] + pow(r, g[1], m) * h[0]) % m, (g[1] + h[1]) % q)

    def inv(g):
        jj = (-g[1]) % q
        return ((-pow(r, jj, m) * g[0]) % m, jj)
    for g in elts:
        assert mul(g, inv(g)) == (0, 0)
    return elts, mul, inv


def mono(g):
    t = []
    for var, e in zip("xy", g):
        if e:
            t.append(var if e == 1 else f"{var}^{e}")
    return " ".join(t) if t else "e"


def rank(rows):
    piv = {}
    for v in rows:
        while v:
            p = v.bit_length() - 1
            if p not in piv:
                piv[p] = v
                break
            v ^= piv[p]
    return len(piv)


def k_code(group, a, b, K, left_a=True):
    """(n, k) of C(G/K; a, b); K a list of central elements forming a subgroup."""
    elts, mul, inv = group
    for z in K:
        assert all(mul(z, g) == mul(g, z) for g in elts)
    coset, reps = {}, []
    for g in elts:
        if g not in coset:
            for z in K:
                coset[mul(g, z)] = len(reps)
            reps.append(g)
    N = len(reps)
    HX, HZ = [], []
    for c in reps:
        if left_a:   # X-check: s^-1 c (block 1), c t^-1 (block 2); Z-check: c t (block 1), s c (block 2)
            rx = [coset[mul(inv(s), c)] for s in a] + [N + coset[mul(c, inv(t))] for t in b]
            rz = [coset[mul(c, t)] for t in b] + [N + coset[mul(s, c)] for s in a]
        else:        # a acting from the right, b from the left
            rx = [coset[mul(c, inv(s))] for s in a] + [N + coset[mul(inv(t), c)] for t in b]
            rz = [coset[mul(t, c)] for t in b] + [N + coset[mul(c, s)] for s in a]
        vx = vz = 0
        for qb in rx:
            vx ^= 1 << qb
        for qb in rz:
            vz ^= 1 << qb
        HX.append(vx)
        HZ.append(vz)
    assert all(bin(u & v).count("1") % 2 == 0 for u in HX for v in HZ)
    return 2 * N, 2 * N - rank(HX) - rank(HZ)


def main():
    E = (0, 0)
    lines = []
    # Qian-Li: x = (1, 0), s = (0, 1) in Z84 x|_29 Z4
    QL = semidirect(84, 4, 29)
    A = [E, (1, 1), (3, 0)]
    B = [(0, 3), (2, 0), (4, 0)]
    for left in (True, False):
        conv = "H_X = [L(A) | R(B)]" if left else "H_X = [R(A) | L(B)]"
        for name, K in (("K = <x^42>", [E, (42, 0)]), ("K trivial", [E])):
            n, k = k_code(QL, A, B, K, left)
            lines.append(f"Qian-Li A, B over Z84 x|_29 Z4, {name}, {conv}: n = {n}, k = {k}")
    # [[672,28,24]]: x = (1, 0), y = (0, 1) in Z12 x|_5 Z28
    G = semidirect(12, 28, 5)
    elts, mul, _ = G
    a = [E, (2, 9), (10, 11)]
    b = [E, (10, 22), (7, 27)]
    centre = [g for g in elts if all(mul(g, h) == mul(h, g) for h in elts)]
    invol = [g for g in centre if g != E and mul(g, g) == E]
    lines.append(f"G = Z12 x|_5 Z28: centre of order {len(centre)}, central involutions "
                 + ", ".join(mono(z) for z in invol))
    n, k = k_code(G, a, b, [E])
    lines.append(f"[[672,28,24]] = C(G; a, b): n = {n}, k = {k}")
    for z in invol:
        n, k = k_code(G, a, b, [E, z])
        lines.append(f"C(G/<{mono(z)}>; a, b): n = {n}, k = {k}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
