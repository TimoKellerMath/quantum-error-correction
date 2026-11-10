#!/usr/bin/env python3
"""Rebuild the ten codes of this directory from their definitions and check them.

    python3 src/rebuild.py [--board DIR]

Run from the directory amplified-and-folded-codes/ (Python 3.8 or later, standard library only).

Amplified codes (CSS): the check matrices are rebuilt from the base entry in data/bases/ (a copy of the
board's codes/<base>.json), the list of dropped base rows and the amplifier length N, by Definition 1 of
README.md, and compared with data/codes/<slug>.json row by row. Checked: the base (k by ranks, witnesses),
that the retained base rows are independent and span the row spaces, CSS orthogonality, n, k by ranks,
n = (N+1) n_base - k_base, k = (N-2) k_base, the check weight w and its formula, qubit degrees, both
witnesses (in the kernel of the other side, outside the row space of their own side, weight 2 d_base, equal
to base witness (x) (e_0 + e_1)), and connectivity of the Tanner graph.

Folds (stabilizer): the parent two-block code is rebuilt from its group and supports, the ZX-involution from
its formula, and the fold from Definition 6 of README.md; compared with data/codes/<slug>.json generator by
generator. Checked: the involution (fixed-point-free, maps every X-check onto a Z-check), isotropy, k by
rank, generator weights and Y counts, the witness and every word of data/witnesses/<slug>.json (commutes
with all generators, outside the stabilizer group, Pauli weight), the generator-level test of the board's
verifier for local Clifford equivalence to a CSS code and condition (i) of Lemma 9, the automorphisms used
by the searches, and the instance files in instances/ (generators, logical representatives). The four other
folds of the gross code (instances/gross_c0.code to gross_c3.code) are rebuilt as well.

--board DIR: DIR is a checkout of github.com/unitaryfoundation/qldpc-challenge. Then data/bases/*.json are
compared with DIR/codes/<base>.json, the fold parents with DIR/codes/144-12-12.json and
DIR/codes/510-16-26.json, and, where DIR/codes/<slug>.json exists, the submitted codes with it.

Exit status 0 if every check passes. Every line of the output starts with PASS, FAIL or INFO.
"""
import argparse
import hashlib
import itertools
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
AMPLIFIED = ["880-40-22", "866-128-12", "610-40-22", "696-48-24", "606-48-22", "656-64-24", "618-48-24",
             "616-56-22"]
FOLDS = {"72-6-9-s": "gross_c4", "255-8-22": "B510_c0"}
FAILURES = []


def report(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        FAILURES.append(what)
    return ok


def info(what):
    print("INFO " + what)


def load(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return json.load(f)


def popcount(v):
    return bin(v).count("1")


def vec(support):
    v = 0
    for q in support:
        v |= 1 << q
    return v


def support(v):
    out, q = [], 0
    while v:
        if v & 1:
            out.append(q)
        v >>= 1
        q += 1
    return out


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


def checks_sha256(checks):
    """The board's binding hash (verify/check_certs.py): sorted supports, sorted rows, compact JSON."""
    if "S" in checks:
        canon = {"S": sorted([sorted(g["X"]), sorted(g["Z"])] for g in checks["S"])}
    else:
        canon = {s: sorted(sorted(r) for r in checks[s]) for s in ("X", "Z")}
    return hashlib.sha256(json.dumps(canon, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


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


# ------------------------------------------------------------------------------------------- amplified codes

def amplify(n, hx, hz, N):
    """Definition 1 of README.md: [LGCE26, Eqs. (14), (15)] with G_X = G_Z = (1, ..., 1) of length N.

    Qubits: central (q, a) -> N q + a; side t_j -> N n + j for the j-th retained Z-check (register Q_2 (x) A_0);
    side s_i -> N n + m_Z + i for the i-th retained X-check (register Q_0 (x) A_2)."""
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
    k = n - rank(X) - rank(Z)
    wx, wz = b["witness"]["X"], b["witness"]["Z"]
    ok &= k == b["k"] and len(wx) == b["d_X"] and len(wz) == b["d_Z"]
    ok &= css_logical(wx, Span(X), Z) and css_logical(wz, Span(Z), X)
    return report(ok, f"base {b['slug']}: CSS, k = {k} by ranks, witnesses of weight {len(wx)}/{len(wz)} are "
                      f"nontrivial logicals; certificate: {b['certificate']}")


def check_amplified(slug, board):
    c = load("data", "codes", f"{slug}.json")
    con = c["construction"]
    b = load("data", "bases", f"{con['base']}.json")
    N, n0, k0 = con["N"], b["n"], b["k"]
    tag = f"{slug} (base {con['base']}, N = {N})"
    check_base(b)
    keepX = [i for i in range(len(b["checks"]["X"])) if i not in set(con["drop_X"])]
    keepZ = [i for i in range(len(b["checks"]["Z"])) if i not in set(con["drop_Z"])]
    hx = [b["checks"]["X"][i] for i in keepX]
    hz = [b["checks"]["Z"][i] for i in keepZ]
    rX, rZ = rank(vec(r) for r in b["checks"]["X"]), rank(vec(r) for r in b["checks"]["Z"])
    report(rank(vec(r) for r in hx) == len(hx) == rX and rank(vec(r) for r in hz) == len(hz) == rZ,
           f"{tag}: retained base rows independent and spanning ({len(hx)} of {len(b['checks']['X'])} X-checks, "
           f"{len(hz)} of {len(b['checks']['Z'])} Z-checks)")
    n, HX, HZ = amplify(n0, hx, hz, N)
    report(n == c["n"] and HX == c["checks"]["X"] and HZ == c["checks"]["Z"],
           f"{tag}: rebuilt check matrices equal data/codes/{slug}.json row by row")
    X = [vec(r) for r in HX]
    Z = [vec(r) for r in HZ]
    report(all(popcount(x & z) % 2 == 0 for x in X for z in Z), f"{tag}: CSS condition H_X H_Z^T = 0")
    SX, SZ = Span(X), Span(Z)
    k = n - len(SX) - len(SZ)
    report(k == c["k"] == (N - 2) * k0 and n == (N + 1) * n0 - k0,
           f"{tag}: n = {n} = (N+1) n_base - k_base, k = {k} = (N-2) k_base (ranks {len(SX)}, {len(SZ)})")
    degX, degZ = [0] * n0, [0] * n0
    for r in hx:
        for q in r:
            degX[q] += 1
    for r in hz:
        for q in r:
            degZ[q] += 1
    w = max(len(r) for r in HX + HZ)
    wf = max(max(len(r) for r in hx) + 1, max(len(r) for r in hz) + 1, N + max(degZ), N + max(degX))
    qx, qz = [0] * n, [0] * n
    for r in HX:
        for q in r:
            qx[q] += 1
    for r in HZ:
        for q in r:
            qz[q] += 1
    report(w == wf, f"{tag}: check weight w = {w} = max(w_X + 1, w_Z + 1, N + q_Z, N + q_X) (Lemma 2(2)); "
                    f"largest qubit degree X/Z/total {max(qx)}/{max(qz)}/{max(a + z for a, z in zip(qx, qz))}")
    d = 2 * b["d_X"]
    lift = {s: sorted(N * q + a for q in b["witness"][s] for a in (0, 1)) for s in "XZ"}
    ok = lift["X"] == c["witness"]["X"] and lift["Z"] == c["witness"]["Z"]
    ok &= css_logical(c["witness"]["X"], SX, Z) and css_logical(c["witness"]["Z"], SZ, X)
    ok &= len(c["witness"]["X"]) == 2 * b["d_X"] and len(c["witness"]["Z"]) == 2 * b["d_Z"] == d == c["d"]
    report(ok, f"{tag}: witnesses u (x) (e_0 + e_1) of weight {2 * b['d_X']}/{2 * b['d_Z']} are nontrivial "
               f"logicals; d = 2 d_base = {d} by Proposition 4")
    report(connected(n, HX + HZ), f"{tag}: Tanner graph connected")
    info(f"{slug}: checks_sha256 {checks_sha256(c['checks'])}")
    if board:
        compare_base(b, board)
        compare_entry(c, board)


# ------------------------------------------------------------------------------------------- folds

def group(mods):
    elts = list(itertools.product(*[range(m) for m in mods]))
    idx = {g: i for i, g in enumerate(elts)}
    add = lambda g, h: tuple((a + b) % m for a, b, m in zip(g, h, mods))
    neg = lambda g: tuple((-a) % m for a, m in zip(g, mods))

    def hom(images):
        def f(g):
            out = tuple(0 for _ in mods)
            for coef, im in zip(g, images):
                out = add(out, tuple((coef * a) % m for a, m in zip(im, mods)))
            return out
        return f

    return elts, idx, add, neg, hom


def parent_checks(con):
    elts, idx, add, neg, _ = group(con["group"])
    N = len(elts)
    FA = [tuple(a) for a in con["F_A"]]
    FB = [tuple(b) for b in con["F_B"]]
    hx = [sorted([idx[add(g, a)] for a in FA] + [N + idx[add(g, b)] for b in FB]) for g in elts]
    hz = [sorted([idx[add(g, neg(b))] for b in FB] + [N + idx[add(g, neg(a))] for a in FA]) for g in elts]
    return N, hx, hz


def involution(con):
    """The qubit permutation Psi of the parent (qubit (1, h) -> idx h, (2, h) -> N + idx h)."""
    elts, idx, add, neg, hom = group(con["group"])
    N = len(elts)
    inv = con["involution"]
    psi = hom([tuple(v) for v in inv["psi"]])
    images = [psi(g) for g in elts]
    if len(set(images)) != N:
        raise ValueError("psi is not bijective")
    perm = [0] * (2 * N)
    if inv["type"] == "preserve":
        c1, c2 = tuple(inv["c1"]), tuple(inv["c2"])
        for g in elts:
            perm[idx[g]] = idx[add(c1, psi(g))]
            perm[N + idx[g]] = N + idx[add(c2, psi(g))]
    else:
        c = tuple(inv["c"])
        pinv = {psi(g): g for g in elts}
        for g in elts:
            perm[idx[g]] = N + idx[add(c, psi(g))]
            perm[N + idx[g]] = idx[pinv[add(g, neg(c))]]
    return perm


def fold_generators(hx, perm):
    """Definition 6 of README.md: Sigma = smaller element of each orbit of Psi, fold qubit = position in Sigma;
    the X-check a gives X on the sigma in a and Z on the sigma with Psi(sigma) in a."""
    sigma = sorted(q for q in range(len(perm)) if q < perm[q])
    pos = {q: i for i, q in enumerate(sigma)}
    S = []
    for a in hx:
        S.append({"X": sorted(pos[q] for q in a if q in pos), "Z": sorted(pos[perm[q]] for q in a if perm[q] in pos)})
    return sigma, S


def pvec(p, n):
    return vec(p["X"]) | (vec(p["Z"]) << n)


def sympl(u, v, n):
    m = (1 << n) - 1
    return (popcount((u & m) & (v >> n)) + popcount((u >> n) & (v & m))) % 2


def pauli_weight(v, n):
    m = (1 << n) - 1
    return popcount((v & m) | (v >> n))


def lc_css_generator_level(S, n):
    """The board's generator-level test: unknowns tau_i (generator i becomes X- or Z-type), one equation
    tau_i + tau_j = [P_iq != P_jq] per pair of generators acting on a common qubit q; a qubit carrying X, Y and
    Z is infeasible outright. Returns (feasible, reason)."""
    on = [[] for _ in range(n)]
    for i, g in enumerate(S):
        xs, zs = set(g["X"]), set(g["Z"])
        for q in xs | zs:
            on[q].append((i, "Y" if q in xs and q in zs else ("X" if q in xs else "Z")))
    parent = list(range(len(S)))
    par = [0] * len(S)

    def find(a):
        if parent[a] == a:
            return a, 0
        r, p = find(parent[a])
        parent[a], par[a] = r, par[a] ^ p
        return r, par[a]

    for q in range(n):
        if len({t for _, t in on[q]}) == 3:
            return False, f"qubit {q} carries X, Y and Z"
        if not on[q]:
            continue
        i0, t0 = on[q][0]
        for i, t in on[q][1:]:
            want = int(t != t0)
            (ri, pi), (r0, p0) = find(i), find(i0)
            if ri == r0:
                if pi ^ p0 != want:
                    return False, f"odd cycle through qubit {q} (generators {i0}, {i})"
            else:
                parent[ri], par[ri] = r0, pi ^ p0 ^ want
    return True, "the parity system has a solution"


def condition_i(S, n):
    """dim {s in <S> : supp s within supp s_i} = rank S - rank(S restricted to the complement of supp s_i)."""
    rows = [pvec(g, n) for g in S]
    r = rank(rows)
    dims = set()
    for g in S:
        keep = ((1 << n) - 1) & ~vec(set(g["X"]) | set(g["Z"]))
        mask = keep | (keep << n)
        dims.add(r - rank(v & mask for v in rows))
    return dims


def read_instance(name, n):
    with open(os.path.join(ROOT, "instances", name), encoding="ascii") as f:
        head = f.readline().split()
        nn, r, L = map(int, head)
        rows = [f.readline().strip() for _ in range(r + L)]
    if nn != n or any(len(s) != 2 * n for s in rows):
        raise ValueError(f"instances/{name}: bad header or row length")
    vs = [int(s[::-1], 2) for s in rows]
    return vs[:r], vs[r:]


def check_fold(slug, board):
    c = load("data", "codes", f"{slug}.json")
    con = c["construction"]
    N, hx, hz = parent_checks(con)
    n = c["n"]
    tag = f"{slug} (fold of {con['parent']})"
    PX = [vec(r) for r in hx]
    PZ = [vec(r) for r in hz]
    kp = 2 * N - rank(PX) - rank(PZ)
    report(all(popcount(x & z) % 2 == 0 for x in PX for z in PZ) and kp == 2 * c["k"],
           f"{tag}: parent rebuilt, CSS, n = {2 * N}, k = {kp} = 2k")
    perm = involution(con)
    zset = {tuple(r) for r in hz}
    ok = all(perm[perm[q]] == q and perm[q] != q for q in range(2 * N))
    ok &= {tuple(sorted(perm[q] for q in r)) for r in hx} == zset and len(zset) == len(hz)
    report(ok, f"{tag}: Psi is a fixed-point-free involution mapping the set of X-check supports onto the set of "
               "Z-check supports")
    sigma, S = fold_generators(hx, perm)
    report(len(sigma) == n == N and S == c["checks"]["S"],
           f"{tag}: rebuilt generators equal data/codes/{slug}.json generator by generator")
    G = [pvec(g, n) for g in S]
    report(all(sympl(u, v, n) == 0 for u, v in itertools.combinations(G, 2)), f"{tag}: isotropic")
    SS = Span(G)
    k = n - len(SS)
    report(k == c["k"], f"{tag}: rank S = {len(SS)}, k = {k}")
    wts = sorted({popcount(vec(set(g["X"]) | set(g["Z"]))) for g in S})
    ys = sum(len(set(g["X"]) & set(g["Z"])) for g in S)
    info(f"{slug}: generator weights {wts}, {ys} Y factors in the generators, "
         f"{sum(1 for g in S if set(g['X']) & set(g['Z']))} generators with a Y")
    w = c["witness"]
    v = pvec(w, n)
    report(all(sympl(v, g, n) == 0 for g in G) and v not in SS and pauli_weight(v, n) == c["d"],
           f"{tag}: witness commutes with S, lies outside S, Pauli weight {pauli_weight(v, n)} = d")
    words = load("data", "witnesses", f"{slug}.json")
    ok = all(all(sympl(pvec(p, n), g, n) == 0 for g in G) and pvec(p, n) not in SS
             and pauli_weight(pvec(p, n), n) == words["weight"] for p in words["words"])
    report(ok and words["count"] == len(words["words"]) and w in words["words"],
           f"{tag}: all {len(words['words'])} words of data/witnesses/{slug}.json are nontrivial logicals of "
           f"Pauli weight {words['weight']}")
    feas, why = lc_css_generator_level(S, n)
    dims = condition_i(S, n)
    report(not feas and dims == {1}, f"{tag}: generator-level LC-CSS test infeasible ({why}); condition (i) holds "
                                     f"for every generator (dimensions {sorted(dims)}), so by Lemma 9 no local "
                                     "Clifford makes the code CSS")
    gens, L = read_instance(FOLDS[slug] + ".code", n)
    ok = gens == G and len(L) == 2 * k
    ok &= all(sympl(l, g, n) == 0 for l in L for g in G) and rank(G + L) == len(SS) + len(L)
    inst = load("instances", FOLDS[slug] + ".json")
    ok &= [pvec(p, n) for p in inst["S"]] == G and [pvec(p, n) for p in inst["logicals"]] == L
    ok &= inst["involution"] == perm and inst["sigma"] == sigma
    report(ok, f"{tag}: instances/{FOLDS[slug]}.code and .json: generators as above, {len(L)} logical "
               "representatives commuting with S and independent modulo S")
    if slug == "255-8-22":
        gset = set(G)
        for name, mp in (("translation h -> h + 5", lambda h: (h + 5) % n), ("multiplier h -> 101 h",
                                                                                  lambda h: (101 * h) % n)):
            img = {pvec({"X": [mp(q) for q in g["X"]], "Z": [mp(q) for q in g["Z"]]}, n) for g in S}
            report(img == gset, f"{tag}: the {name} maps the generator set onto itself")
    info(f"{slug}: checks_sha256 {checks_sha256(c['checks'])}")
    if board:
        b = board_doc(board, con["parent"])
        same = b is not None and {tuple(r) for r in b["checks"]["X"]} == {tuple(r) for r in hx} \
            and {tuple(r) for r in b["checks"]["Z"]} == {tuple(r) for r in hz}
        report(same, f"{tag}: parent equals codes/{con['parent']}.json of the board (same qubit labels)")
        compare_entry(c, board)


def check_other_gross_folds():
    """The block-swapping folds of the gross code with psi = -1 and c in {(0,0), (0,1), (1,0), (1,1)}
    (instances/gross_c0.code to gross_c3.code; their distances are in the logs), and the weight-6 logicals of the
    first one (data/witnesses/gross_c0_d6.json)."""
    con = dict(load("data", "codes", "72-6-9-s.json")["construction"])
    N, hx, hz = parent_checks(con)
    zset = {tuple(r) for r in hz}
    for i, c in enumerate([[0, 0], [0, 1], [1, 0], [1, 1]]):
        con["involution"] = {"type": "swap", "psi": [[11, 0], [0, 5]], "c": c}
        perm = involution(con)
        ok = all(perm[perm[q]] == q and perm[q] != q for q in range(2 * N))
        ok &= {tuple(sorted(perm[q] for q in r)) for r in hx} == zset
        _, S = fold_generators(hx, perm)
        G = [pvec(g, N) for g in S]
        gens, L = read_instance(f"gross_c{i}.code", N)
        SS = Span(G)
        ok &= gens == G and len(SS) == 66 and len(L) == 12 and rank(G + L) == 78
        ok &= all(sympl(l, g, N) == 0 for l in L for g in G)
        msg = f"gross fold c = {tuple(c)} (Psi(1,h) = (2, c - h)): rebuilt, equal to instances/gross_c{i}.code, k = 6"
        if i == 0:
            words = load("data", "witnesses", "gross_c0_d6.json")
            ok &= all(all(sympl(pvec(p, N), g, N) == 0 for g in G) and pvec(p, N) not in SS
                      and pauli_weight(pvec(p, N), N) == 6 for p in words["words"])
            msg += f"; the {len(words['words'])} words of data/witnesses/gross_c0_d6.json are nontrivial of weight 6"
        report(ok, msg)


# ------------------------------------------------------------------------------------------- board comparison

def board_doc(board, slug):
    p = os.path.join(board, "codes", f"{slug}.json")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def compare_base(b, board):
    d = board_doc(board, b["slug"])
    ok = d is not None and d["checks"]["X"] == b["checks"]["X"] and d["checks"]["Z"] == b["checks"]["Z"] \
        and d["distance"]["X"]["witness"] == b["witness"]["X"] and d["distance"]["Z"]["witness"] == b["witness"]["Z"]
    report(ok, f"data/bases/{b['slug']}.json equals codes/{b['slug']}.json of the board (checks, witnesses)")


def compare_entry(c, board):
    d = board_doc(board, c["slug"])
    if d is None:
        info(f"codes/{c['slug']}.json is not in the board checkout")
        return
    report(checks_sha256(d["checks"]) == checks_sha256(c["checks"]) and d["n"] == c["n"] and d["k"] == c["k"],
           f"codes/{c['slug']}.json of the board has the checks of data/codes/{c['slug']}.json")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--board", help="checkout of github.com/unitaryfoundation/qldpc-challenge")
    a = ap.parse_args()
    for slug in AMPLIFIED:
        check_amplified(slug, a.board)
    for slug in FOLDS:
        check_fold(slug, a.board)
    check_other_gross_folds()
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} checks failed")
        return 1
    print("RESULT: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
