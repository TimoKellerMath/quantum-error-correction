#!/usr/bin/env python3
"""Rebuild the three codes of this directory from their definitions and check them.

    python3 src/rebuild.py [--board DIR]

Run from the directory two-block-and-amplified-codes/ (Python 3.8 or later, standard library only).

480-16-26 (Part A of README.md): the group G = C12 x Dic5 from its multiplication rule (group axioms by Light's
test), the code C(G; a, b) from a and b, compared with data/codes/480-16-26.json row by row. Checked: n, k by
ranks, CSS orthogonality, row and column weights, the parity statement (the all-one vector is the sum of the rows
of H_X and of H_Z), both witnesses, connectivity of the Tanner graph; the instances of the exact searches
(instances/c480_16.*, instances/c480v.*: the same rows after the stated relabelling, the same witnesses) and
their automorphism generators (each maps the rows of H_X and of H_Z onto themselves; the generated group has
order 240 and two orbits of 240 qubits); the tropidist problems of src/tropidist_inputs.py (H the check matrix,
L a basis of the logical classes); the decomposition of k over the blocks of F_2[<z^4 x>] (all 16 logical qubits
in one block); the local ring O/(N(a), N(b)) at that block (length 4, Hilbert function (1,1,1,1)) and the residue
ranks of the Clifford symbols (Proposition 3); k of the quotients and of the cover (Proposition 4).

860-130-10, 463-34-12 (Part B): the check matrices are rebuilt from the base entry in data/bases/ (a copy of
the board's codes/<base>.json) by Definition 6 of README.md and compared with data/codes/<slug>.json row by
row. Checked: the base (k by ranks, witnesses), that the base rows are independent, CSS orthogonality, n, k by
ranks, n = 5 n_base - k_base, k = 2 k_base, the check weight w and its formula, qubit degrees, both witnesses
(nontrivial, weight 2 d_base, equal to base witness (x) (e_0 + e_1)), connectivity, the tropidist problems, and
the instances of the exact searches (instances/amplified/: relabelled bases, the amplified codes, automorphism
generators).

--board DIR: DIR is a checkout of github.com/unitaryfoundation/qldpc-challenge. Then data/bases/*.json are
compared with DIR/codes/<base>.json, and, where DIR/codes/<slug>.json exists, the submitted codes with it.

Exit status 0 if every check passes. Every line of the output starts with PASS, FAIL or INFO.
"""
import argparse
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AMPLIFIED = ["860-130-10", "463-34-12"]
FAILURES = []


def report(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        FAILURES.append(what)
    return ok


def info(what):
    print("INFO " + what)


def path(*parts):
    return os.path.join(ROOT, *parts)


def load(*parts):
    with open(path(*parts), encoding="utf-8") as f:
        return json.load(f)


def popcount(v):
    return bin(v).count("1")


def vec(support):
    v = 0
    for q in support:
        v |= 1 << q
    return v


class Span:
    """Row space over GF(2), rows as Python integers, reduced by leading bit."""

    def __init__(self, rows=()):
        self.piv = {}
        for r in rows:
            self.add(r)

    def reduce(self, v):
        while v:
            p = self.piv.get(v.bit_length() - 1)
            if p is None:
                return v
            v ^= p
        return 0

    def add(self, v):
        v = self.reduce(v)
        if v:
            self.piv[v.bit_length() - 1] = v
        return bool(v)

    def __contains__(self, v):
        return self.reduce(v) == 0

    def __len__(self):
        return len(self.piv)


def rank(rows):
    return len(Span(rows))


def connected(n, rows):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for r in rows:
        for q in r[1:]:
            a, b = find(r[0]), find(q)
            if a != b:
                parent[a] = b
    return len({find(q) for q in range(n)}) == 1


def css_logical(w, own, other_rows):
    """w is a nontrivial logical of its type: even overlap with every opposite check, not in the own row space."""
    v = vec(w)
    return all(popcount(v & r) % 2 == 0 for r in other_rows) and v not in own


def css_params(n, X, Z):
    """X, Z lists of supports: (orthogonal, rank X, rank Z, k)."""
    MX, MZ = [vec(r) for r in X], [vec(r) for r in Z]
    orth = all(popcount(x & z) % 2 == 0 for x in MX for z in MZ)
    rx, rz = rank(MX), rank(MZ)
    return orth, rx, rz, n - rx - rz


def read_inst(fn):
    """Instance format of the search programs: 'n m', then one check per line (0-based qubits)."""
    with open(path(fn), encoding="ascii") as f:
        n, m = map(int, f.readline().split())
        rows = [list(map(int, ln.split())) for ln in f if ln.strip()]
    assert len(rows) == m
    return n, rows


def read_ints(fn):
    with open(path(fn), encoding="ascii") as f:
        return [int(t) for t in f.read().split()]


def read_gens(fn, n):
    """Generators in exactA's format (count, then permutations) or exactB's ('#' comment lines, then
    permutations); the image of qubit q is in position q."""
    with open(path(fn), encoding="ascii") as f:
        lines = [ln.split() for ln in f if ln.strip() and not ln.startswith("#")]
    if lines and len(lines[0]) == 1:
        cnt, lines = int(lines[0][0]), lines[1:]
        assert cnt == len(lines), fn
    gens = [list(map(int, ln)) for ln in lines]
    assert all(sorted(g) == list(range(n)) for g in gens), fn
    return gens


def rowset(rows):
    return {frozenset(r) for r in rows}


def preserves(g, rows):
    s = rowset(rows)
    return {frozenset(g[q] for q in r) for r in rows} == s


def generated(gens, n, limit=100000):
    """All elements of the permutation group generated by gens (breadth-first closure)."""
    e = tuple(range(n))
    seen, todo = {e}, [e]
    while todo:
        p = todo.pop()
        for g in gens:
            h = tuple(g[p[q]] for q in range(n))
            if h not in seen:
                seen.add(h)
                todo.append(h)
                if len(seen) > limit:
                    return seen
    return seen


def orbits(gens, n):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for g in gens:
        for q in range(n):
            a, b = find(q), find(g[q])
            if a != b:
                parent[a] = b
    sizes = {}
    for q in range(n):
        sizes[find(q)] = sizes.get(find(q), 0) + 1
    return sorted(sizes.values())


# ------------------------------------------------------------------------------------------- Part A: the group

class Group:
    """C_m x H with H = <x, y | x^p, y^J, y x y^-1 = x^-1> (p = 5: J = 4 gives Dic5, J = 2 gives D5; p = 1 kills
    x). Element z^c x^i y^j is the triple (c, i, j), c mod m, i mod p, j mod J, with index (p c + i) J + j; for
    C12 x Dic5 this is 20 c + 4 i + j. Product: (c,i,j)(c',i',j') = (c+c', i + (-1)^j i', j+j'), as
    y^j x^i' y^-j = x^((-1)^j i')."""

    def __init__(self, m, p, J):
        self.m, self.p, self.J = m, p, J
        self.els = [(c, i, j) for c in range(m) for i in range(p) for j in range(J)]
        self.order = len(self.els)
        self.index = {g: t for t, g in enumerate(self.els)}

    def norm(self, g):
        c, i, j = g
        return (c % self.m, i % self.p, j % self.J)

    def mul(self, g, h):
        (c, i, j), (c2, i2, j2) = g, h
        return self.norm((c + c2, i + (i2 if j % 2 == 0 else -i2), j + j2))

    def inv(self, g):
        c, i, j = g
        return self.norm((-c, i if j % 2 else -i, -j))

    def idx(self, g):
        return self.index[self.norm(g)]


def two_block(G, a, b):
    """H_X = (L_a | R_b), H_Z = (R_b^T | L_a^T): X-check c has block-1 qubits s^-1 c (s in a) and block-2 qubits
    c t^-1 (t in b); Z-check c has block-1 qubits c t (t in b) and block-2 qubits s c (s in a). Qubit (block e, g)
    has index e |G| + idx(g); checks in the order of idx(c)."""
    N = G.order
    HX = [sorted([G.idx(G.mul(G.inv(s), c)) for s in a] + [N + G.idx(G.mul(c, G.inv(t))) for t in b])
          for c in G.els]
    HZ = [sorted([G.idx(G.mul(c, t)) for t in b] + [N + G.idx(G.mul(s, c)) for s in a]) for c in G.els]
    return 2 * N, HX, HZ


def check_group(G):
    z, x, y, e = (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 0, 0)
    gens = [z, x, y]
    # Light's associativity test: (g h) s = g (h s) for all g, h and all s in a generating set
    light = all(G.mul(G.mul(g, h), s) == G.mul(g, G.mul(h, s)) for g in G.els for h in G.els for s in gens)
    seen, todo = {e}, [e]
    while todo:
        g = todo.pop()
        for s in gens:
            h = G.mul(g, s)
            if h not in seen:
                seen.add(h)
                todo.append(h)
    ident = all(G.mul(e, g) == g == G.mul(g, e) for g in G.els)
    inverses = all(G.mul(g, G.inv(g)) == e == G.mul(G.inv(g), g) for g in G.els)

    def power(g, k):
        h = e
        for _ in range(k):
            h = G.mul(h, g)
        return h

    def order(g):
        k, h = 1, g
        while h != e:
            h, k = G.mul(h, g), k + 1
        return k

    rel = (power(x, 5) == e and power(y, 4) == e and power(z, 12) == e
           and G.mul(G.mul(y, x), G.inv(y)) == G.inv(x) and all(G.mul(z, g) == G.mul(g, z) for g in G.els))
    dic = [g for g in G.els if g[0] == 0]
    invol = [g for g in dic if g != e and G.mul(g, g) == e]
    centre = [g for g in G.els if all(G.mul(g, h) == G.mul(h, g) for h in G.els)]
    report(light and len(seen) == 240 and ident and inverses and rel and invol == [(0, 0, 2)]
           and [order(z), order(x), order(y)] == [12, 5, 4] and len(centre) == 24,
           "480-16-26: G = C12 x Dic5 is a group of order 240 (Light's test, identity, inverses), generated by "
           "z, x, y with z^12 = x^5 = y^4 = 1, y x y^-1 = x^-1, z central; <x, y> of order 20 has the single "
           "involution y^2 (so it is Dic5); |Z(G)| = 24")


A_ELEMS = [(0, 0, 0), (8, 3, 2), (11, 0, 3)]      # a = 1 + z^8 x^3 y^2 + z^11 y^3
B_ELEMS = [(0, 0, 0), (2, 1, 2), (0, 0, 1)]       # b = 1 + z^2 x y^2 + y


def rep_index(q):
    """Labelling of instances/c480_16.*: index 120 e + 2 (5 c + i) + j2 of z^c x^i (y^2)^j2 y^e in each block."""
    blk, t = divmod(q, 240)
    e, rest = divmod(t, 120)
    ci, j2 = divmod(rest, 2)
    c, i = divmod(ci, 5)
    return 240 * blk + 20 * c + 4 * i + 2 * j2 + e


def check_480():
    G = Group(12, 5, 4)
    check_group(G)
    c = load("data", "codes", "480-16-26.json")
    n, HX, HZ = two_block(G, A_ELEMS, B_ELEMS)
    report(n == c["n"] == 480 and HX == c["checks"]["X"] and HZ == c["checks"]["Z"],
           "480-16-26: rebuilt check matrices equal data/codes/480-16-26.json row by row")
    orth, rx, rz, k = css_params(n, HX, HZ)
    report(orth and rx == rz == 232 and k == 16 == c["k"], f"480-16-26: CSS, rank H_X = {rx}, rank H_Z = {rz}, "
                                                            f"k = {k}")
    colw = {sum(q in set(r) for r in rows) for rows in (HX, HZ) for q in range(n)}
    report({len(r) for r in HX + HZ} == {6} and colw == {3},
           "480-16-26: every check has weight 6 (3 + 3), every qubit lies in 3 X-checks and 3 Z-checks")
    ones = (1 << n) - 1
    sx = sz = 0
    for r in HX:
        sx ^= vec(r)
    for r in HZ:
        sz ^= vec(r)
    report(sx == ones and sz == ones, "480-16-26: the sum of the rows of H_X and the sum of the rows of H_Z are the "
                                      "all-one vector, so every logical operator has even weight (Lemma 1)")
    MX, MZ = [vec(r) for r in HX], [vec(r) for r in HZ]
    wx, wz = c["witness"]["X"], c["witness"]["Z"]
    report(len(wx) == len(wz) == 26 and css_logical(wx, Span(MX), MZ) and css_logical(wz, Span(MZ), MX),
           "480-16-26: both witnesses are nontrivial logicals of weight 26 (X: in ker H_Z, not in rowsp H_X; Z: "
           "in ker H_X, not in rowsp H_Z)")
    report(connected(n, HX + HZ), "480-16-26: Tanner graph connected")
    sX, sZ = rowset(HX), rowset(HZ)
    # instance of the first exact searches (labelling rep_index)
    perm = [rep_index(q) for q in range(n)]
    report(sorted(perm) == list(range(n)), "480-16-26: rep_index is a bijection")
    _, ix = read_inst("instances/c480_16.hx")
    _, iz = read_inst("instances/c480_16.hz")
    report({frozenset(perm[q] for q in r) for r in ix} == sX and {frozenset(perm[q] for q in r) for r in iz} == sZ
           and len(ix) == len(iz) == 240, "480-16-26: instances/c480_16.{hx,hz} are the code in the labelling "
                                          "rep_index (same row sets)")
    for side, w in (("X", wx), ("Z", wz)):
        report(sorted(perm[q] for q in read_ints(f"instances/c480_16.w{side}")) == sorted(w),
               f"480-16-26: instances/c480_16.w{side} is the {side}-witness of the entry")
    # instance of the second exact searches (seeded relabelling q -> p[q])
    p = read_ints("instances/c480v.perm")
    _, vx = read_inst("instances/c480v.hx")
    _, vz = read_inst("instances/c480v.hz")
    report(sorted(p) == list(range(n)) and rowset(vx) == {frozenset(p[q] for q in r) for r in HX}
           and rowset(vz) == {frozenset(p[q] for q in r) for r in HZ} and len(vx) == len(vz) == 240,
           "480-16-26: instances/c480v.{hx,hz} are the code relabelled by q -> p[q], p = instances/c480v.perm")
    for side, w in (("X", wx), ("Z", wz)):
        report(read_ints(f"instances/c480v.w{side}") == sorted(p[q] for q in w),
               f"480-16-26: instances/c480v.w{side} is the relabelled {side}-witness")
    for inst, fx, fz, files in (("c480_16", ix, iz, ("c480_16.gens", "c480_16.gensB")),
                                ("c480v", vx, vz, ("c480v.gensA", "c480v.gensB"))):
        for fn in files:
            gens = read_gens(f"instances/{fn}", n)
            ok = all(preserves(g, fx) and preserves(g, fz) for g in gens)
            grp = generated(gens, n)
            report(ok and len(grp) == 240 and orbits(gens, n) == [240, 240],
                   f"480-16-26: instances/{fn}: {len(gens)} generators, each maps the rows of H_X and of H_Z of "
                   f"{inst} onto themselves; group of order {len(grp)}, qubit orbits {orbits(gens, n)}")
    check_tropidist_inputs("480-16-26", n, HX, HZ, 16)
    check_blocks(G)
    check_local_ring()
    check_tower()


def check_tropidist_inputs(slug, n, HX, HZ, k):
    """The problems written by src/tropidist_inputs.py: side Z is (H_X, a basis of the X-logical classes), side X
    mirrored; tropidist's least weight of e with H e = 0 and L e != 0 is then d_Z resp. d_X."""
    sys.path.insert(0, os.path.join(ROOT, "src"))
    import tropidist_inputs as T
    MX, MZ = [vec(r) for r in HX], [vec(r) for r in HZ]
    for side, H, comm, stab in (("Z", MX, MZ, MX), ("X", MZ, MX, MZ)):
        nn, Hm, Lm = T.problem(slug, side)
        sp = Span(stab)
        ok = nn == n and Hm == H and len(Lm) == k
        ok &= all(all(popcount(v & r) % 2 == 0 for r in comm) for v in Lm)
        for v in Lm:
            ok &= sp.add(v)
        report(ok, f"{slug}: tropidist problem, side {side}: H = H_{'X' if side == 'Z' else 'Z'} row by row, L = {k} "
                   f"{'X' if side == 'Z' else 'Z'}-logicals independent modulo the stabilizers (all logical classes)")


# ---------------------------------------------------------------------- Part A: blocks, local ring, quotients

def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def pdivmod(a, b):
    q = 0
    db = b.bit_length()
    while a.bit_length() >= db:
        s = a.bit_length() - db
        q ^= 1 << s
        a ^= b << s
    return q, a


def pinv(a, f):
    """inverse of a modulo f in F_2[t] (extended Euclid)."""
    r0, r1, s0, s1 = f, pdivmod(a, f)[1], 0, 1
    while r1 != 1:
        q, r = pdivmod(r0, r1)
        r0, r1, s0, s1 = r1, r, s1, s0 ^ pmul(q, s1)
        assert r1 != 0, "not invertible"
    return pdivmod(s1, f)[1]


FACTORS_15 = {"t+1": 0b11, "t^2+t+1": 0b111, "t^4+t+1": 0b10011, "t^4+t^3+1": 0b11001,
              "t^4+t^3+t^2+t+1": 0b11111}


def check_blocks(G):
    """k of C(G;a,b) e for the primitive idempotents e of F_2[<t>], t = z^4 x of order 15 (central idempotents of
    F_2[G], since conjugation by y maps t to t^4 and fixes every idempotent of F_2[t]/(t^15+1))."""
    P = (1 << 15) | 1
    prod = 1
    for f in FACTORS_15.values():
        prod = pmul(prod, f)
    report(prod == P, "480-16-26: t^15 + 1 = (t+1)(t^2+t+1)(t^4+t+1)(t^4+t^3+1)(t^4+t^3+t^2+t+1) over F_2")
    t = (4, 1, 0)
    tpow = [(0, 0, 0)]
    for _ in range(14):
        tpow.append(G.mul(tpow[-1], t))
    report(G.mul(tpow[-1], t) == (0, 0, 0) and len(set(tpow)) == 15
           and G.mul(G.mul((0, 0, 1), t), G.inv((0, 0, 1))) == tpow[4], "480-16-26: t = z^4 x has order 15 and "
                                                                       "y t y^-1 = t^4")
    N = G.order
    rmul = [[G.idx(G.mul(g, h)) for h in G.els] for g in G.els]      # rmul[g][h] = index of g h

    def alg(elems):
        v = 0
        for g in elems:
            v ^= 1 << G.idx(g)
        return v

    def amul(u, v):
        r = 0
        us = [i for i in range(N) if u >> i & 1]
        vs = [i for i in range(N) if v >> i & 1]
        for i in us:
            for j in vs:
                r ^= 1 << rmul[i][j]
        return r

    a, b = alg(A_ELEMS), alg(B_ELEMS)
    idem = {}
    total = 0
    for name, f in FACTORS_15.items():
        cof = pdivmod(P, f)[0]
        e = pdivmod(pmul(cof, pinv(cof, f)), P)[1]
        idem[name] = e
        E = 0
        for kk in range(15):
            if e >> kk & 1:
                E ^= 1 << G.idx(tpow[kk])
        eg = [amul(E, 1 << i) for i in range(N)]
        dim = rank(eg)
        d1 = rank([amul(a, u) for u in eg] + [amul(u, b) for u in eg])          # (u, v) -> a u + v b on eR^2
        d2 = rank([amul(u, b) | (amul(a, u) << N) for u in eg])                # w -> (w b, a w) on eR
        kk_ = 2 * dim - d1 - d2
        total += kk_
        info(f"480-16-26: block of {name} (t = z^4 x): dim e F_2[G] = {dim}, k = {kk_}")
        idem[name + " k"] = kk_
    s = 0
    for name in FACTORS_15:
        s ^= idem[name]
    orth = all(pdivmod(pmul(idem[p1], idem[p2]), P)[1] == (idem[p1] if p1 == p2 else 0)
               for p1 in FACTORS_15 for p2 in FACTORS_15)
    report(s == 1 and orth and total == 16 and idem["t^4+t+1 k"] == 16,
           "480-16-26: the five idempotents are orthogonal with sum 1, and all 16 logical qubits lie in the block of "
           "t^4 + t + 1, i.e. at the Galois orbit O_1 of the characters chi of order 15 with chi(z^4 x) a root of "
           "t^4 + t + 1 (Proposition 3)")


# F_16 = F_2[g]/(g^4 + g + 1)
def fmul(a, b):
    r = pmul(a, b)
    for s in (6, 5, 4):
        if r >> s & 1:
            r ^= 0b10011 << (s - 4)
    return r


def fpow(a, k):
    r = 1
    for _ in range(k):
        r = fmul(r, a)
    return r


def finv(a):
    return next(b for b in range(1, 16) if fmul(a, b) == 1)


def fbasis(rows):
    """a basis (reduced echelon form) over F_16 of the span of vectors given as lists."""
    rows = [list(r) for r in rows]
    rk, col = 0, 0
    ncol = len(rows[0]) if rows else 0
    while rk < len(rows) and col < ncol:
        piv = next((i for i in range(rk, len(rows)) if rows[i][col]), None)
        if piv is None:
            col += 1
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        iv = finv(rows[rk][col])
        rows[rk] = [fmul(iv, v) for v in rows[rk]]
        for i in range(len(rows)):
            if i != rk and rows[i][col]:
                f = rows[i][col]
                rows[i] = [u ^ fmul(f, v) for u, v in zip(rows[i], rows[rk])]
        rk += 1
        col += 1
    return rows[:rk]


def frank(rows):
    """rank over F_16 of vectors given as lists."""
    return len(fbasis(rows))


def check_local_ring():
    """Proposition 3: O = K[A_2], K = F_16, A_2 = <z^3> x <y^2> = C4 x C2 (basis element z^3u y^2v at 2u + v);
    chi the character of A_2' = <z^4> x <x> with chi(z^4) = zeta_3 = g^5, chi(x) = zeta_5 = g^3, extended to
    F_2[A] -> O by z^c x^i y^2j -> chi(z^4c x^i) z^(9c) y^2j (z^(4c) the odd part, z^(9c) the 2-part of z^c).
    a = alpha + y beta, b = gamma + delta y with alpha, beta, gamma, delta in F_2[A]; N(a) = alpha tau(alpha) +
    y^2 beta tau(beta), N(b) = gamma tau(gamma) + y^2 delta tau(delta), tau(x) = x^-1."""
    g = 0b10
    z3, z5 = fpow(g, 5), fpow(g, 3)
    report(fpow(g, 15) == 1 and fpow(g, 5) != 1 and fpow(g, 3) != 1 and fpow(z3, 3) == 1 and fpow(z5, 5) == 1,
           "480-16-26: g (g^4 = g + 1) has order 15 in F_16; zeta_3 = g^5, zeta_5 = g^3")

    def chi_elem(c, i, j2, inv=False):
        """image of z^c x^i y^(2 j2) in O as an 8-vector over K (inv: the character chi^tau)."""
        v = [0] * 8
        u = (9 * c) % 12 // 3
        val = fmul(fpow(z3, c % 3), fpow(z5, (-i if inv else i) % 5))
        v[2 * u + j2] = val
        return v

    def omul(p, q):
        r = [0] * 8
        for s in range(8):
            if p[s]:
                for t_ in range(8):
                    if q[t_]:
                        u = ((s // 2) + (t_ // 2)) % 4
                        w = (s % 2 + t_ % 2) % 2
                        r[2 * u + w] ^= fmul(p[s], q[t_])
        return r

    def oadd(p, q):
        return [u ^ v for u, v in zip(p, q)]

    def image(elems, inv=False):
        v = [0] * 8
        for (c, i, j2) in elems:
            v = oadd(v, chi_elem(c, i, j2, inv))
        return v

    def tau(elems):
        return [(c, -i % 5, j2) for (c, i, j2) in elems]

    def amul_A(p, q):        # product in F_2[A] of lists of (c, i, j2), A abelian
        out = {}
        for (c, i, j) in p:
            for (c2, i2, j2) in q:
                key = ((c + c2) % 12, (i + i2) % 5, (j + j2) % 2)
                out[key] = out.get(key, 0) ^ 1
        return [k for k, v in out.items() if v]

    # a = 1 + z^8 x^3 y^2 + z^11 y^3: alpha = 1 + z^8 x^3 y^2, beta = z^11 y^2 (z^11 y^3 = y z^11 y^2)
    # b = 1 + z^2 x y^2 + y: gamma = 1 + z^2 x y^2, delta = 1
    alpha, beta = [(0, 0, 0), (8, 3, 1)], [(11, 0, 1)]
    gamma, delta = [(0, 0, 0), (2, 1, 1)], [(0, 0, 0)]
    y2 = [(0, 0, 1)]
    Na = amul_A(alpha, tau(alpha))
    for e in amul_A(y2, amul_A(beta, tau(beta))):
        Na = [x for x in Na if x != e] if e in Na else Na + [e]
    Nb = amul_A(gamma, tau(gamma))
    for e in amul_A(y2, amul_A(delta, tau(delta))):
        Nb = [x for x in Nb if x != e] if e in Nb else Nb + [e]
    na, nb = image(Na), image(Nb)

    def residue(v):
        r = 0
        for c in v:
            r ^= c
        return r

    # Clifford symbols (README.md, "Where the 16 logical qubits sit"): P_a = [[alpha(chi), y^2 beta(chi^tau)],
    # [beta(chi), alpha(chi^tau)]], P_b = [[gamma(chi), y^2 delta(chi)], [delta(chi^tau), gamma(chi^tau)]]; residues
    # modulo the maximal ideal (the residue of y^2 is 1)
    Ra = [[residue(image(alpha)), residue(image(beta, True))], [residue(image(beta)), residue(image(alpha, True))]]
    Rb = [[residue(image(gamma)), residue(image(delta))], [residue(image(delta, True)), residue(image(gamma, True))]]
    det_a = fmul(Ra[0][0], Ra[1][1]) ^ fmul(Ra[0][1], Ra[1][0])
    det_b = fmul(Rb[0][0], Rb[1][1]) ^ fmul(Rb[0][1], Rb[1][0])
    report(frank(Ra) == 1 and frank(Rb) == 1 and det_a == residue(na) == 0 and det_b == residue(nb) == 0,
           "480-16-26: at chi the residue matrices of the Clifford symbols of a and b have rank 1, and their "
           "determinants are the residues of N(a)(chi), N(b)(chi), which vanish")
    basis = []
    for s in range(8):
        e = [0] * 8
        e[s] = 1
        basis.append(e)
    ideal = [omul(na, e) for e in basis] + [omul(nb, e) for e in basis]
    dimI = frank(ideal)
    one = basis[0]
    m1 = [oadd(e, one) for e in basis[1:]]                 # augmentation ideal, spanned by h - 1
    # Hilbert function h_i = dim (m^i + I) - dim (m^(i+1) + I), for i = 0, 1, ... until m^(i+1) lies in I
    hil, power, prev = [], m1, 8
    for _ in range(10):
        d = frank(ideal + power)
        hil.append(prev - d)
        prev = d
        if d == dimI:
            break
        power = fbasis([omul(p, q) for p in power for q in m1])
    length = 8 - dimI
    report(length == 4 and hil == [1, 1, 1, 1],
           f"480-16-26: O/(N(a)(chi), N(b)(chi)) has length {length} and Hilbert function {tuple(hil)}, so it is "
           "K[eps]/(eps^4), a curvilinear point; with [F_4 : F_2] = 2 and the two characters chi, chi^tau of the "
           "tau-orbit this gives k = 2 * 2 * 4 = 16 (Proposition 3)")


def check_tower():
    """Proposition 4: the same supports over C_m x Dic5 (z-exponents modulo m; for m = 24 the same exponents),
    over G/<x> = C12 x C4 and over G/<y^2> = C12 x D5."""
    out = []
    for m, p, J, lab in ((3, 5, 4, "C3 x Dic5"), (6, 5, 4, "C6 x Dic5"), (12, 5, 4, "C12 x Dic5"),
                         (24, 5, 4, "C24 x Dic5"), (4, 5, 4, "C4 x Dic5 = G/<z^4>"), (12, 1, 4, "C12 x C4 = G/<x>"),
                         (12, 5, 2, "C12 x D5 = G/<y^2>")):
        G = Group(m, p, J)
        n, HX, HZ = two_block(G, A_ELEMS, B_ELEMS)
        orth, rx, rz, k = css_params(n, HX, HZ)
        assert orth
        out.append((lab, n, k))
        info(f"480-16-26: two-block code with the same supports over {lab}: n = {n}, k = {k}")
    ks = [k for _, _, k in out]
    report(ks == [4, 8, 16, 16, 0, 0, 8], "480-16-26: k = 4, 8, 16, 16 over C_m x Dic5 for m = 3, 6, 12, 24; k = 0 "
                                          "over G/<z^4> and G/<x>; k = 8 over G/<y^2> (Proposition 4)")


# ------------------------------------------------------------------------------------- Part B: amplified codes

def amplify(n, hx, hz, N):
    """Definition 6 of README.md: [LGCE26, Eqs. (14), (15)] with G_X = G_Z = (1, ..., 1) of length N.

    Qubits: central (q, a) -> N q + a; side t_j -> N n + j for the j-th Z-check (register Q_2 (x) A_0);
    side s_i -> N n + m_Z + i for the i-th X-check (register Q_0 (x) A_2)."""
    mX, mZ = len(hx), len(hz)
    xcols = [[] for _ in range(n)]
    zcols = [[] for _ in range(n)]
    for i, r in enumerate(hx):
        for q in r:
            xcols[q].append(i)
    for j, r in enumerate(hz):
        for q in r:
            zcols[q].append(j)
    HX = [sorted([N * q + a for q in hx[i]] + [N * n + mZ + i]) for i in range(mX) for a in range(N)]
    HX += [sorted([N * q + a for a in range(N)] + [N * n + j for j in zcols[q]]) for q in range(n)]
    HZ = [sorted([N * q + a for q in hz[j]] + [N * n + j]) for j in range(mZ) for a in range(N)]
    HZ += [sorted([N * q + a for a in range(N)] + [N * n + mZ + i for i in xcols[q]]) for q in range(n)]
    return N * n + mX + mZ, HX, HZ


def check_base(b):
    n = b["n"]
    X = [vec(r) for r in b["checks"]["X"]]
    Z = [vec(r) for r in b["checks"]["Z"]]
    ok = all(popcount(x & z) % 2 == 0 for x in X for z in Z)
    rx, rz = rank(X), rank(Z)
    k = n - rx - rz
    wx, wz = b["witness"]["X"], b["witness"]["Z"]
    ok &= k == b["k"] and len(wx) == b["d_X"] and len(wz) == b["d_Z"]
    ok &= css_logical(wx, Span(X), Z) and css_logical(wz, Span(Z), X)
    ok &= rx == len(X) and rz == len(Z)
    return report(ok, f"base {b['slug']}: CSS, k = {k} by ranks, the {len(X)} X-checks and {len(Z)} Z-checks are "
                      f"linearly independent, witnesses of weight {len(wx)}/{len(wz)} are nontrivial logicals")


def check_amplified(slug, board):
    c = load("data", "codes", f"{slug}.json")
    con = c["construction"]
    b = load("data", "bases", f"{con['base']}.json")
    N, n0, k0 = con["N"], b["n"], b["k"]
    tag = f"{slug} (base {con['base']}, N = {N})"
    check_base(b)
    hx, hz = b["checks"]["X"], b["checks"]["Z"]
    n, HX, HZ = amplify(n0, hx, hz, N)
    report(n == c["n"] and HX == c["checks"]["X"] and HZ == c["checks"]["Z"],
           f"{tag}: rebuilt check matrices equal data/codes/{slug}.json row by row")
    orth, rx, rz, k = css_params(n, HX, HZ)
    report(orth, f"{tag}: CSS orthogonality H_X H_Z^T = 0")
    report(n == (N + 1) * n0 - k0 and k == (N - 2) * k0 == c["k"],
           f"{tag}: n = {n} = {N + 1}*{n0} - {k0}, k = {k} = {N - 2}*{k0} by ranks ({rx}, {rz})")
    w = max(len(r) for r in HX + HZ)
    wX0, wZ0 = max(len(r) for r in hx), max(len(r) for r in hz)
    qX0 = max(sum(q in set(r) for r in hx) for q in range(n0))
    qZ0 = max(sum(q in set(r) for r in hz) for q in range(n0))
    degX = [0] * n
    degZ = [0] * n
    for r in HX:
        for q in r:
            degX[q] += 1
    for r in HZ:
        for q in r:
            degZ[q] += 1
    report(w == max(wX0 + 1, wZ0 + 1, N + qZ0, N + qX0) == c["w"],
           f"{tag}: check weight {w} = max(w_X + 1, w_Z + 1, N + q_Z, N + q_X) (Lemma 7)")
    info(f"{tag}: qubit degree X {max(degX)}, Z {max(degZ)}, total {max(x + z for x, z in zip(degX, degZ))}")
    MX, MZ = [vec(r) for r in HX], [vec(r) for r in HZ]
    for side, own, other in (("X", MX, MZ), ("Z", MZ, MX)):
        wit = c["witness"][side]
        lift = sorted(N * q + a for q in b["witness"][side] for a in (0, 1))
        dP = b[f"d_{side}"]
        report(wit == lift and len(wit) == 2 * dP and css_logical(wit, Span(own), other),
               f"{tag}: {side}-witness = base witness (x) (e_0 + e_1), weight {len(wit)} = 2*{dP}, nontrivial logical")
    report(connected(n, HX + HZ), f"{tag}: Tanner graph connected")
    check_tropidist_inputs(slug, n, HX, HZ, k)
    if board:
        bp = os.path.join(board, "codes", f"{con['base']}.json")
        with open(bp, encoding="utf-8") as f:
            bd = json.load(f)
        report(bd["checks"]["X"] == b["checks"]["X"] and bd["checks"]["Z"] == b["checks"]["Z"]
               and bd["distance"]["X"]["witness"] == b["witness"]["X"]
               and bd["distance"]["Z"]["witness"] == b["witness"]["Z"],
               f"{tag}: data/bases/{con['base']}.json equals the board's codes/{con['base']}.json (checks, witnesses)")
    check_amp_instances(slug, con["base"], n0, b, n, HX, HZ)


def check_amp_instances(slug, base, n0, b, n, HX, HZ):
    d = "instances/amplified"
    bX, bZ = rowset(b["checks"]["X"]), rowset(b["checks"]["Z"])
    for tool, stem in (("A", f"A/{base}_r1"), ("B", f"B/{base}")):
        p = read_ints(f"{d}/{stem}.perm")
        _, ix = read_inst(f"{d}/{stem}.hx")
        _, iz = read_inst(f"{d}/{stem}.hz")
        ok = sorted(p) == list(range(n0)) and rowset(ix) == {frozenset(p[q] for q in r) for r in bX} \
            and rowset(iz) == {frozenset(p[q] for q in r) for r in bZ}
        gens = read_gens(f"{d}/{stem}.gens", n0)
        ok &= all(preserves(g, ix) and preserves(g, iz) for g in gens)
        report(ok, f"{slug}: {d}/{stem}.{{hx,hz}} are the base {base} relabelled by q -> p[q] ({stem}.perm); "
                   f"{len(gens)} automorphism generators")
    _, ax = read_inst(f"{d}/amp422_{base}.hx")
    _, az = read_inst(f"{d}/amp422_{base}.hz")
    ok = rowset(ax) == rowset(HX) and rowset(az) == rowset(HZ) and len(ax) == len(HX) and len(az) == len(HZ)
    for tool in ("A", "B"):
        gens = read_gens(f"{d}/{tool}/amp422_{base}.gens", n)
        grp = generated(gens, n)
        ok &= all(preserves(g, ax) and preserves(g, az) for g in gens) and len(grp) == 24
    report(ok, f"{slug}: {d}/amp422_{base}.{{hx,hz}} have the rows of the code (same labelling); the automorphism "
               "generators of both programs preserve them and generate a group of order 24 (the permutations of "
               "the four copies)")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--board", help="checkout of the qLDPC challenge board repository")
    a = ap.parse_args()
    check_480()
    for slug in AMPLIFIED:
        check_amplified(slug, a.board)
    if a.board:
        for slug in ["480-16-26"] + AMPLIFIED:
            bp = os.path.join(a.board, "codes", f"{slug}.json")
            if os.path.exists(bp):
                with open(bp, encoding="utf-8") as f:
                    bd = json.load(f)
                c = load("data", "codes", f"{slug}.json")
                report(bd["checks"]["X"] == c["checks"]["X"] and bd["checks"]["Z"] == c["checks"]["Z"]
                       and bd["distance"]["X"]["witness"] == c["witness"]["X"]
                       and bd["distance"]["Z"]["witness"] == c["witness"]["Z"],
                       f"{slug}: the board's codes/{slug}.json has the checks and witnesses of data/codes/{slug}.json")
            else:
                info(f"{slug}: no codes/{slug}.json in {a.board}")
    print("RESULT: all checks passed" if not FAILURES else f"RESULT: {len(FAILURES)} checks FAILED")
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
