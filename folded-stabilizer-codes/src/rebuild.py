#!/usr/bin/env python3
"""Rebuild the five codes of this directory from their definitions and check them.

    python3 src/rebuild.py [--board DIR]

Run from the directory folded-stabilizer-codes/ (Python 3.8 or later, standard library only).

For each code: the parent two-block code is rebuilt from its group and supports (Definition 1 of README.md) and
compared, after the relabelling of the group elements stored in data/parents/<parent>.json, with the parent's
check matrices copied from the board; the ZX-duality Psi_c is built from its formula (Definition 2) and checked
(fixed-point-free involution mapping the X-check supports onto the Z-check supports); the fold (Definition 3) is
compared with data/codes/<slug>.json generator by generator and with the Khesin-Lu mirror form (Remark 5).
Checked on the fold: isotropy, k by rank (k = k_parent / 2), generator weights, Y factors and qubit degrees, the
witness and every word of data/witnesses/<slug>.json (commutes with every generator, outside the stabilizer group,
Pauli weight d; its image under lambda is a nontrivial Z-logical of the parent of weight 2y + t, Lemma 4), the
generator-level test of the board's verifier for local Clifford equivalence to a CSS code, the group-level
criterion of Lemma 6, the board's two connectivity rules, and that instances/<slug>.txt and magma/<slug>.m encode
these generators.

--board DIR: DIR is a checkout of github.com/unitaryfoundation/qldpc-challenge. Then data/parents/*.json are compared
with DIR/codes/<parent>.json and, where DIR/codes/<slug>.json exists, the codes with their board entries.

Exit status 0 if every check passes. Every line of the output starts with PASS, FAIL or INFO.
"""
import argparse
import hashlib
import itertools
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CODES = ["60-9-8", "56-6-8-s", "32-9-6", "48-6-8-s", "33-2-7"]
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


# ------------------------------------------------------------------------------------------- GF(2) helpers

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


def rref_rows(rows):
    """The rows of the reduced row echelon form, columns in increasing bit order (bit j = column j)."""
    piv = {}
    for r in rows:
        v = r
        for c, p in piv.items():
            if v >> c & 1:
                v ^= p
        if v:
            c = (v & -v).bit_length() - 1
            for c2 in list(piv):
                if piv[c2] >> c & 1:
                    piv[c2] ^= v
            piv[c] = v
    return [piv[c] for c in sorted(piv)]


def pvec(p, n):
    """A Pauli operator {"X": [...], "Z": [...]} as the integer x | z << n."""
    return vec(p["X"]) | (vec(p["Z"]) << n)


def sympl(u, v, n):
    m = (1 << n) - 1
    return (popcount((u & m) & (v >> n)) + popcount((u >> n) & (v & m))) % 2


def pauli_weight(v, n):
    m = (1 << n) - 1
    return popcount((v & m) | (v >> n))


def checks_sha256(checks):
    """The board's binding hash (verify/check_certs.py): sorted supports, sorted generators, compact JSON."""
    canon = {"S": sorted([sorted(g["X"]), sorted(g["Z"])] for g in checks["S"])}
    return hashlib.sha256(json.dumps(canon, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


# ------------------------------------------------------------------------------------------- groups, parents, folds

def group(orders):
    """G = Z/o_1 x ... x Z/o_r; elements in lexicographic order of the coordinates (iota = position)."""
    elts = list(itertools.product(*[range(o) for o in orders]))
    idx = {g: i for i, g in enumerate(elts)}

    def add(g, h):
        return tuple((a + b) % o for a, b, o in zip(g, h, orders))

    def neg(g):
        return tuple((-a) % o for a, o in zip(g, orders))

    return elts, idx, add, neg


def parent_checks(con):
    """Definition 1: X-check g on (1, g + F_A) and (2, g + F_B), Z-check g on (1, g - F_B) and (2, g - F_A);
    qubit (1, h) -> iota(h), (2, h) -> |G| + iota(h). Rows in the order of g."""
    elts, idx, add, neg = group(con["group"])
    N = len(elts)
    FA = [tuple(a) for a in con["F_A"]]
    FB = [tuple(b) for b in con["F_B"]]
    hx = [sorted([idx[add(g, a)] for a in FA] + [N + idx[add(g, b)] for b in FB]) for g in elts]
    hz = [sorted([idx[add(g, neg(b))] for b in FB] + [N + idx[add(g, neg(a))] for a in FA]) for g in elts]
    return N, hx, hz


def board_relabelling(par, orders):
    """The map from the board's qubit labels of the parent to ours: board qubit i < N is (1, E_i), board qubit
    N + i is (2, E_i), where E_i is the group element with board index i (data/parents/<parent>.json)."""
    _, idx, _, _ = group(orders)
    E = [tuple(e) for e in par["labels"]["elements"]]
    N = len(E)
    if sorted(E) != sorted(idx) or len(set(E)) != N:
        raise ValueError(f"{par['slug']}: the label list is not a bijection onto G")
    return [idx[e] for e in E] + [N + idx[e] for e in E]


def involution(con):
    """Definition 2: Psi_c(1, h) = (2, c - h), Psi_c(2, h) = (1, c - h), as a permutation of the parent's qubits."""
    elts, idx, add, neg = group(con["group"])
    N = len(elts)
    c = tuple(con["involution"]["c"])
    perm = [0] * (2 * N)
    for h in elts:
        perm[idx[h]] = N + idx[add(c, neg(h))]
        perm[N + idx[h]] = idx[add(c, neg(h))]
    return perm


def fold_generators(hx, perm):
    """Definition 3: Sigma = the qubits q with q < Psi(q), in increasing order, fold qubit i = sigma_i; the X-check a
    gives X on {i : sigma_i in a} and Z on {i : Psi(sigma_i) in a}, in the order of the X-checks."""
    sigma = sorted(q for q in range(len(perm)) if q < perm[q])
    pos = {q: i for i, q in enumerate(sigma)}
    S = [{"X": sorted(pos[q] for q in a if q in pos), "Z": sorted(pos[perm[q]] for q in a if perm[q] in pos)}
         for a in hx]
    return sigma, S


def mirror_generators(con):
    """Remark 5: the mirror code of Khesin and Lu with A = c - F_B and B = F_A, generator Z(A + g) X(B - g) for g in
    G, in the order of -g (so that it lines up with the X-check g)."""
    elts, idx, add, neg = group(con["group"])
    c = tuple(con["involution"]["c"])
    A = [add(c, neg(tuple(b))) for b in con["F_B"]]
    B = [tuple(a) for a in con["F_A"]]
    out = []
    for g in elts:
        h = neg(g)
        out.append({"X": sorted(idx[add(b, neg(h))] for b in B), "Z": sorted(idx[add(a, h)] for a in A)})
    return out


# ------------------------------------------------------------------------------------------- local Clifford tests

def letters_on_qubits(S, n):
    on = [[] for _ in range(n)]
    for i, g in enumerate(S):
        xs, zs = set(g["X"]), set(g["Z"])
        for q in xs | zs:
            on[q].append((i, "Y" if q in xs and q in zs else ("X" if q in xs else "Z")))
    return on


def parity_system(S, n, use):
    """The board's generator-level test, restricted to the generators in `use`: unknowns t_i (generator i becomes
    X- or Z-type), t_i = t_j for two generators with the same letter on a common qubit, t_i != t_j for different
    letters; a qubit on which the generators of `use` carry X, Y and Z is infeasible outright. Returns
    (feasible, reason)."""
    on = letters_on_qubits(S, n)
    parent = {i: i for i in use}
    par = {i: 0 for i in use}

    def find(a):
        path = []
        while parent[a] != a:
            path.append(a)
            a = parent[a]
        acc = 0
        for v in reversed(path):
            acc ^= par[v]
            par[v] = acc
            parent[v] = a
        return a, (par[path[0]] if path else 0)

    for q in range(n):
        here = [(i, t) for i, t in on[q] if i in parent]
        if len({t for _, t in here}) == 3:
            return False, f"qubit {q} carries X, Y and Z"
        if not here:
            continue
        i0, t0 = here[0]
        for i, t in here[1:]:
            want = int(t != t0)
            (ri, pi), (r0, p0) = find(i), find(i0)
            if ri == r0:
                if pi ^ p0 != want:
                    return False, f"an odd cycle closes at qubit {q} with generators {i0} and {i}"
            else:
                parent[ri], par[ri] = r0, pi ^ p0 ^ want
    return True, "the parity system has a solution"


def isolated_generators(S, n):
    """The generators s with {t in <S> : supp t within supp s} = {1, s} (condition (i) of Lemma 6), computed as
    dim = rank S - rank(S restricted to the complement of supp s) = 1."""
    rows = [pvec(g, n) for g in S]
    r = rank(rows)
    out = []
    for i, g in enumerate(S):
        keep = ((1 << n) - 1) & ~vec(set(g["X"]) | set(g["Z"]))
        mask = keep | (keep << n)
        if r - rank(v & mask for v in rows) == 1:
            out.append(i)
    return out


def tanner_components(S, n):
    """Components of the generator/qubit incidence graph, isolated vertices included (board: tanner_connected)."""
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    empty = 0
    for g in S:
        sup = sorted(set(g["X"]) | set(g["Z"]))
        if not sup:
            empty += 1
        for q in sup[1:]:
            parent[find(q)] = find(sup[0])
    return len({find(q) for q in range(n)}) + empty


def symplectic_blocks(S, n):
    """Disjoint qubit blocks of the stabilizer group: union of the supports of the rows of the reduced row echelon
    form of (A | B), columns X_0..X_{n-1}, Z_0..Z_{n-1} (board: stabilizer_group_connected)."""
    m = (1 << n) - 1
    rows = [support((v & m) | (v >> n)) for v in rref_rows(pvec(g, n) for g in S)]
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for sup in rows:
        for q in sup[1:]:
            parent[find(q)] = find(sup[0])
    return len({find(q) for q in range(n)})


# ------------------------------------------------------------------------------------------- input files

def read_mitm_instance(path):
    """instances/<slug>.txt: "n m", then m lines "nx x_1 .. x_nx nz z_1 .. z_nz"."""
    with open(path, encoding="ascii") as f:
        tok = f.read().split()
    n, m = int(tok[0]), int(tok[1])
    pos, out = 2, []
    for _ in range(m):
        nx = int(tok[pos])
        X = [int(t) for t in tok[pos + 1:pos + 1 + nx]]
        pos += 1 + nx
        nz = int(tok[pos])
        Z = [int(t) for t in tok[pos + 1:pos + 1 + nz]]
        pos += 1 + nz
        out.append({"X": X, "Z": Z})
    if pos != len(tok):
        raise ValueError(f"{path}: trailing data")
    return n, out


MAGMA_ROW = re.compile(r"for q in \[([0-9, ]*)\] do v\[q\+1\] := 1; end for; for q in \[([0-9, ]*)\] do "
                       r"v\[n\+q\+1\] := 1; end for; Append\(~rows, v\);")


def read_magma_script(path):
    with open(path, encoding="ascii") as f:
        text = f.read()
    n = int(re.search(r"^n := (\d+);", text, re.M).group(1))
    gens = [{"X": [int(t) for t in a.replace(",", " ").split()], "Z": [int(t) for t in b.replace(",", " ").split()]}
            for a, b in MAGMA_ROW.findall(text)]
    return n, gens


# ------------------------------------------------------------------------------------------- the checks

def check_parent(par, con):
    N, hx, hz = parent_checks(con)
    X = [vec(r) for r in hx]
    Z = [vec(r) for r in hz]
    tag = f"parent {par['slug']}"
    k = 2 * N - rank(X) - rank(Z)
    report(all(popcount(x & z) % 2 == 0 for x in X for z in Z) and k == par["k"] and 2 * N == par["n"],
           f"{tag}: rebuilt from G = {' x '.join(f'Z/{o}' for o in con['group'])}, F_A, F_B; CSS, n = {2 * N}, "
           f"k = {k} by ranks")
    pi = board_relabelling(par, con["group"])
    bx = {tuple(sorted(pi[q] for q in r)) for r in par["checks"]["X"]}
    bz = {tuple(sorted(pi[q] for q in r)) for r in par["checks"]["Z"]}
    report(bx == {tuple(r) for r in hx} and bz == {tuple(r) for r in hz} and len(bx) == len(hx) == N
           and len(par["checks"]["X"]) == len(par["checks"]["Z"]) == N,
           f"{tag}: the board's checks (data/parents/{par['slug']}.json), relabelled by {par['labels']['rule']}, "
           "equal the rebuilt X- and Z-checks as sets of rows")
    SX, SZ = Span(X), Span(Z)
    wx = vec(pi[q] for q in par["witness"]["X"])
    wz = vec(pi[q] for q in par["witness"]["Z"])
    ok = all(popcount(wx & z) % 2 == 0 for z in Z) and wx not in SX and popcount(wx) == par["d_X"]
    ok &= all(popcount(wz & x) % 2 == 0 for x in X) and wz not in SZ and popcount(wz) == par["d_Z"]
    report(ok, f"{tag}: the board's witnesses, relabelled, are nontrivial logicals of weight {par['d_X']}/"
               f"{par['d_Z']}; certificate: {par['certificate']}")
    return N, hx, hz, SZ


def check_code(slug, board):
    c = load("data", "codes", f"{slug}.json")
    con = c["construction"]
    par = load("data", "parents", f"{con['parent']}.json")
    n, d = c["n"], c["d"]
    tag = f"{slug} (fold of {con['parent']}, c = {tuple(con['involution']['c'])})"
    N, hx, hz, SZp = check_parent(par, con)
    perm = involution(con)
    zset = {tuple(r) for r in hz}
    ok = all(perm[perm[q]] == q and perm[q] != q for q in range(2 * N))
    ok &= {tuple(sorted(perm[q] for q in r)) for r in hx} == zset and len(zset) == N
    report(ok, f"{tag}: Psi_c is a fixed-point-free involution mapping the X-check supports onto the Z-check "
               "supports")
    sigma, S = fold_generators(hx, perm)
    report(sigma == list(range(N)) and n == N and S == c["checks"]["S"],
           f"{tag}: Sigma is the first block; the rebuilt generators equal data/codes/{slug}.json generator by "
           "generator")
    report(mirror_generators(con) == S, f"{tag}: the generators are the mirror code Z(A + g) X(B - g) with "
                                        "A = c - F_B, B = F_A (Remark 5)")
    G = [pvec(g, n) for g in S]
    report(all(sympl(u, v, n) == 0 for u, v in itertools.combinations(G, 2)), f"{tag}: isotropic")
    SS = Span(G)
    k = n - len(SS)
    report(k == c["k"] == par["k"] // 2, f"{tag}: rank {len(SS)}, k = {k} = k(parent) / 2")
    wts = [len(set(g["X"]) | set(g["Z"])) for g in S]
    ys = [len(set(g["X"]) & set(g["Z"])) for g in S]
    deg = [0] * n
    for g in S:
        for q in set(g["X"]) | set(g["Z"]):
            deg[q] += 1
    hist = {w: wts.count(w) for w in sorted(set(wts))}
    report(max(wts) == c["w"] and sum(ys) == c["y_factors"] and max(deg) == c["qubit_degree"],
           f"{tag}: generator weights {hist} (w = {max(wts)}), {sum(ys)} Y factors in "
           f"{sum(1 for y in ys if y)} generators, qubit degree at most {max(deg)}")
    w = c["witness"]
    v = pvec(w, n)
    report(all(sympl(v, g, n) == 0 for g in G) and v not in SS and pauli_weight(v, n) == d,
           f"{tag}: the witness commutes with S, lies outside S, Pauli weight {pauli_weight(v, n)} = d")
    words = load("data", "witnesses", f"{slug}.json")
    ok = words["count"] == len(words["words"]) and w in words["words"] and words["weight"] == d
    lam_ok = True
    for p in words["words"]:
        u = pvec(p, n)
        ok &= all(sympl(u, g, n) == 0 for g in G) and u not in SS and pauli_weight(u, n) == d
        # Lemma 4: lambda(p) has z_i on sigma_i and x_i on Psi(sigma_i); a nontrivial Z-logical of the parent
        lam = vec(sigma[i] for i in p["Z"]) | vec(perm[sigma[i]] for i in p["X"])
        y, t = len(set(p["X"]) & set(p["Z"])), len(set(p["X"]) ^ set(p["Z"]))
        lam_ok &= all(popcount(lam & r) % 2 == 0 for r in (vec(a) for a in hx)) and lam not in SZp
        lam_ok &= popcount(lam) == 2 * y + t
    report(ok, f"{tag}: all {len(words['words'])} words of data/witnesses/{slug}.json are nontrivial logicals of "
               f"Pauli weight {d}")
    report(lam_ok, f"{tag}: lambda maps each of them to a nontrivial Z-logical of the parent of weight 2y + t "
                   "(Lemma 4)")
    feas, why = parity_system(S, n, range(len(S)))
    report(not feas, f"{tag}: the board's generator-level test finds no single-qubit Cliffords making every "
                     f"generator pure ({why})")
    iso = isolated_generators(S, n)
    feas2, why2 = parity_system(S, n, iso)
    report(not feas2, f"{tag}: {len(iso)} of {len(S)} generators satisfy condition (i) of Lemma 6, and already "
                      f"for these the parity system has no solution ({why2}); so no product of single-qubit "
                      "Cliffords maps the stabilizer group onto a CSS group")
    report(tanner_components(S, n) == 1 and symplectic_blocks(S, n) == 1,
           f"{tag}: the board's connectivity rules: Tanner graph connected, stabilizer group in one block")
    nn, inst = read_mitm_instance(os.path.join(ROOT, "instances", f"{slug}.txt"))
    report(nn == n and inst == S, f"{tag}: instances/{slug}.txt holds these generators")
    nm, mag = read_magma_script(os.path.join(ROOT, "magma", f"{slug}.m"))
    report(nm == n and mag == S, f"{tag}: magma/{slug}.m holds these generators")
    info(f"{slug}: checks_sha256 {checks_sha256(c['checks'])}")
    if board:
        b = board_doc(board, con["parent"])
        same = b is not None and b["checks"]["X"] == par["checks"]["X"] and b["checks"]["Z"] == par["checks"]["Z"] \
            and b["distance"]["X"]["witness"] == par["witness"]["X"] \
            and b["distance"]["Z"]["witness"] == par["witness"]["Z"]
        report(same, f"data/parents/{con['parent']}.json equals codes/{con['parent']}.json of the board (checks, "
                     "witnesses)")
        e = board_doc(board, slug)
        if e is None:
            info(f"codes/{slug}.json is not in the board checkout")
        else:
            report(checks_sha256(e["checks"]) == checks_sha256(c["checks"]) and e["n"] == n and e["k"] == k,
                   f"codes/{slug}.json of the board has the generators of data/codes/{slug}.json")


def board_doc(board, slug):
    p = os.path.join(board, "codes", f"{slug}.json")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--board", help="checkout of github.com/unitaryfoundation/qldpc-challenge")
    a = ap.parse_args()
    for slug in CODES:
        check_code(slug, a.board)
    if FAILURES:
        print(f"RESULT: {len(FAILURES)} checks failed")
        return 1
    print("RESULT: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
