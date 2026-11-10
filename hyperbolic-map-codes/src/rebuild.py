#!/usr/bin/env python3
"""Rebuild the nine codes from their group presentations and from the stored darts, and check them.

    python3 src/rebuild.py [--data DIR] [SLUG ...]

Python 3.8 or later, standard library only. For each code (all codes of data/parameters.json by default):

 (a) Stored darts. data/maps.json holds the dart permutations R, S (right action, 0-based; dart 0 is the identity
     of the group) of the map whose code is in data/codes/; for 720-74-8 those of the base map and the voltage
     of the cover on the base edges. The code of these darts must equal data/codes/<slug>.hx, .hz row for row:
     qubit e = the e-th edge {x, x^(RS)} in the order of its smallest dart, X-check i = the i-th cycle of S
     (vertex), Z-check j = the j-th cycle of R (face), cycles in the order of their smallest dart, incidences
     counted mod 2.
 (b) Presentation. Coset enumeration (Todd-Coxeter, HLT with coincidences) of <r,s | relators> of
     data/presentations.json over the trivial subgroup gives the regular representation of the group; it must
     have order 2E, and the map of (r, s) must be isomorphic to the stored map (a dart bijection commuting with R
     and S); the induced qubit permutation must map the rebuilt check rows onto the stored ones. For 720-74-8 the
     base map is rebuilt this way and the double cover is built from the word list of data/cover_voltage.txt;
     the isomorphism with the stored cover is searched over all roots (a cover is not a regular map).
 (c) Conder's relators. Every relator of Conder's entry that involves only R and S must hold for the stored
     darts in the stated orientation (R -> r, S -> s for "map"; R -> r^-1, S -> s^-1 for "mirror"); for an
     entry of the list of reflexible maps also its mirror image w(R^-1, S^-1).
 (d) H_X H_Z^T = 0, ranks over GF(2), k = n - rank H_X - rank H_Z = 2 - chi, maximal row weights, column weights.
 (e) Witnesses of data/parameters.json: an X-witness has weight d_X, lies in ker H_Z and outside rowsp H_X;
     a Z-witness mirrors it.
Prints one JSON line per code and exits with status 1 if any check fails.
"""
import argparse
import json
import os
import sys
import time
from collections import deque

HERE = os.path.dirname(os.path.abspath(__file__))
INV = (1, 0, 3, 2)                                   # letters: 0 = r, 1 = r^-1, 2 = s, 3 = s^-1


# ---------------------------------------------------------------------------------------------- words
def parse_word(s):
    """'s^-1*r^3*(r*s)^2' -> list of letters 0..3; the letters may be r, s or R, S."""
    s = s.replace(" ", "")
    pos = 0

    def word():
        nonlocal pos
        out = factor()
        while pos < len(s) and s[pos] == "*":
            pos += 1
            out += factor()
        return out

    def factor():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            base = word()
            assert s[pos] == ")", s
            pos += 1
        else:
            c = s[pos].lower()
            assert c in "rs", s
            base = [0 if c == "r" else 2]
            pos += 1
        e = 1
        if pos < len(s) and s[pos] == "^":
            pos += 1
            j = pos
            if s[pos] == "-":
                pos += 1
            while pos < len(s) and s[pos].isdigit():
                pos += 1
            e = int(s[j:pos])
        if e < 0:
            base = [INV[x] for x in reversed(base)]
            e = -e
        return base * e

    w = word()
    assert pos == len(s), s
    return w


def apply_word(w, gens, x):
    """dart reached from x by the word w (right action: letters applied left to right)."""
    for a in w:
        x = gens[a][x]
    return x


def inverse(P):
    Q = [0] * len(P)
    for i, j in enumerate(P):
        Q[j] = i
    return Q


# ---------------------------------------------------------------------------------- coset enumeration
def enumerate_cosets(rels, limit=3000000):
    """HLT coset enumeration over the trivial subgroup with coincidence processing (Holt, Handbook of
    Computational Group Theory, Sec. 5.1). Returns the permutations of r and s on the cosets."""
    T = [-1, -1, -1, -1]
    p = [0]

    def new(c, x):
        d = len(p)
        if d >= limit:
            raise OverflowError("coset limit")
        p.append(d)
        T.extend((-1, -1, -1, -1))
        T[4 * c + x] = d
        T[4 * d + INV[x]] = c
        return d

    def rep(c):
        r = c
        while p[r] != r:
            r = p[r]
        while p[c] != r:
            p[c], c = r, p[c]
        return r

    def merge(a, b, q):
        a, b = rep(a), rep(b)
        if a != b:
            if a > b:
                a, b = b, a
            p[b] = a
            q.append(b)

    def coincidence(a, b):
        q = []
        merge(a, b, q)
        i = 0
        while i < len(q):
            e = q[i]
            i += 1
            for x in range(4):
                f = T[4 * e + x]
                if f != -1:
                    ix = INV[x]
                    T[4 * f + ix] = -1
                    e1, f1 = rep(e), rep(f)
                    if T[4 * e1 + x] != -1:
                        merge(f1, T[4 * e1 + x], q)
                    elif T[4 * f1 + ix] != -1:
                        merge(e1, T[4 * f1 + ix], q)
                    else:
                        T[4 * e1 + x] = f1
                        T[4 * f1 + ix] = e1

    def scan_and_fill(c, w):
        f, b, i, j = c, c, 0, len(w) - 1
        while True:
            while i <= j and T[4 * f + w[i]] != -1:
                f = T[4 * f + w[i]]
                i += 1
            if i > j:
                if f != b:
                    coincidence(f, b)
                return
            while j >= i and T[4 * b + INV[w[j]]] != -1:
                b = T[4 * b + INV[w[j]]]
                j -= 1
            if j < i:
                coincidence(f, b)
                return
            if i == j:
                T[4 * f + w[i]] = b
                T[4 * b + INV[w[i]]] = f
                return
            new(f, w[i])

    c = 0
    while c < len(p):
        if p[c] == c:
            for w in rels:
                scan_and_fill(c, w)
                if p[c] != c:
                    break
            if p[c] == c:
                for x in range(4):
                    if T[4 * c + x] == -1:
                        new(c, x)
        c += 1
    live = [c for c in range(len(p)) if p[c] == c]
    num = {c: i for i, c in enumerate(live)}
    R = [num[T[4 * c + 0]] for c in live]
    S = [num[T[4 * c + 2]] for c in live]
    return R, S, len(p)


# ------------------------------------------------------------------------------------- maps and codes
def cycles(P):
    seen = [False] * len(P)
    out = []
    for x in range(len(P)):
        if not seen[x]:
            c = []
            y = x
            while not seen[y]:
                seen[y] = True
                c.append(y)
                y = P[y]
            out.append(c)
    return out


def edges(R, S):
    """edge label of every dart and the dart pairs; edge partner of x is x^(RS) = S[R[x]]."""
    lab = [-1] * len(R)
    pairs = []
    for x in range(len(R)):
        if lab[x] < 0:
            y = S[R[x]]
            if y == x or S[R[y]] != x:
                raise ValueError("RS is not a fixed-point-free involution")
            lab[x] = lab[y] = len(pairs)
            pairs.append((x, y))
    return lab, pairs


def code(R, S):
    """(hx rows, hz rows, V, E, F, lab, pairs); rows are sorted lists of edge labels (incidence mod 2)."""
    lab, pairs = edges(R, S)

    def rows(cyc):
        out = []
        for c in cyc:
            odd = set()
            for x in c:
                odd ^= {lab[x]}
            out.append(sorted(odd))
        return out
    vc, fc = cycles(S), cycles(R)
    return rows(vc), rows(fc), len(vc), len(pairs), len(fc), lab, pairs


def cover_perms(R, S, dart_voltage, m):
    """darts (x, i) -> x + D i; R~(x, i) = (xR, i + omega(x)), S~(x, i) = (xS, i)."""
    D = len(R)
    Rc = [R[x % D] + D * (((x // D) + dart_voltage[x % D]) % m) for x in range(D * m)]
    Sc = [S[x % D] + D * (x // D) for x in range(D * m)]
    return Rc, Sc


def bfs_key(gens, root):
    """relabel darts in BFS order from root (generators in the given order); returns (key, order)."""
    D = len(gens[0])
    pos = [-1] * D
    pos[root] = 0
    order = [root]
    dq = deque([root])
    while dq:
        x = dq.popleft()
        for g in gens:
            y = g[x]
            if pos[y] < 0:
                pos[y] = len(order)
                order.append(y)
                dq.append(y)
    if len(order) != D:
        return None, None
    return tuple(tuple(pos[g[x]] for x in order) for g in gens), order


def isomorphism(gens1, gens2, roots2):
    """dart bijection phi with phi(g1(x)) = g2(phi(x)) for both generators, or None (root 0 of map 1
    against the given roots of map 2)."""
    k1, o1 = bfs_key(gens1, 0)
    for t in roots2:
        k2, o2 = bfs_key(gens2, t)
        if k2 is not None and k2 == k1:
            phi = [0] * len(o1)
            for a, b in zip(o1, o2):
                phi[a] = b
            return phi
    return None


def relabelled_rows(rows, pi):
    return sorted(sorted(pi[e] for e in r) for r in rows)


# ---------------------------------------------------------------------------------------------- GF(2)
def to_int(row):
    v = 0
    for j in row:
        v ^= 1 << j
    return v


def echelon(rows):
    piv = {}
    for r in rows:
        v = to_int(r)
        while v:
            h = v.bit_length() - 1
            if h not in piv:
                piv[h] = v
                break
            v ^= piv[h]
    return piv


def reduce(piv, v):
    while v:
        h = v.bit_length() - 1
        if h not in piv:
            return v
        v ^= piv[h]
    return 0


def read_rows(path):
    with open(path) as f:
        n, m = map(int, f.readline().split())
        rows = [sorted(map(int, f.readline().split())) for _ in range(m)]
    return n, rows


# ---------------------------------------------------------------------------------------------- checks
def check(slug, data, maps, pres, params, cover_words):
    t0 = time.process_time()
    P, M, C = pres[slug], maps[slug], params[slug]
    n, HX = read_rows(os.path.join(data, "codes", slug + ".hx"))
    n2, HZ = read_rows(os.path.join(data, "codes", slug + ".hz"))
    out = {"slug": slug}
    R, S = M["R"], M["S"]
    D = len(R)
    is_cover = "voltage" in M
    # (a) stored darts -> stored matrices, row for row
    if is_cover:
        lab, _ = edges(R, S)
        R_st, S_st = cover_perms(R, S, [M["voltage"][lab[x]] for x in range(D)], M["m"])
    else:
        R_st, S_st = R, S
    hx, hz, V, E, F, lab_st, pairs_st = code(R_st, S_st)
    out["a_darts_give_matrices"] = (n == n2 == E and hx == HX and hz == HZ)
    out.update(V=V, E=E, F=F, chi=V - E + F)
    # (b) presentation -> regular representation -> map -> code
    rels = [parse_word(w) for w in P["relators"]]
    try:
        Rr, Sr, defined = enumerate_cosets(sorted(rels, key=len))
    except OverflowError:                            # the enumeration did not close within the coset limit
        Rr, Sr, defined = [], [], None
    out["group_order"] = len(Rr) if defined else None
    out["cosets_defined"] = defined
    out["order_is_2E"] = len(Rr) == 2 * P["E"] == D
    if len(Rr) == D:
        phi_base = isomorphism([R, S], [Rr, Sr], [0])
    else:
        phi_base = None
    out["b_base_isomorphic"] = phi_base is not None
    if phi_base is not None:
        if is_cover:
            # voltage from the word list on the rebuilt base, then the cover
            Ri, Si = inverse(Rr), inverse(Sr)
            gens = [Rr, Ri, Sr, Si]
            labr, _ = edges(Rr, Sr)
            wedge = [0] * (D // 2)
            for w in cover_words:
                wedge[labr[apply_word(parse_word(w), gens, 0)]] ^= 1
            R2, S2 = cover_perms(Rr, Sr, [wedge[labr[x]] for x in range(D)], M["m"])
            phi = isomorphism([R_st, S_st], [R2, S2], range(len(R2)))
            out["b_voltage_edges"] = sum(wedge)
        else:
            R2, S2, phi = Rr, Sr, phi_base
        out["b_map_isomorphic"] = phi is not None
        if phi is not None:
            hx2, hz2, _, _, _, lab2, _ = code(R2, S2)
            pi = [lab2[phi[x]] for (x, y) in pairs_st]
            out["b_rows_equal"] = (sorted(pi) == list(range(E)) and relabelled_rows(HX, pi) == sorted(hx2)
                                   and relabelled_rows(HZ, pi) == sorted(hz2))
        else:
            out["b_rows_equal"] = False
    else:
        out["b_map_isomorphic"] = out["b_rows_equal"] = False
    # (c) Conder's relators on the stored (base) darts
    Ri, Si = inverse(R), inverse(S)
    if P["variant"] == "mirror":
        gens = [Ri, R, Si, S]                       # R -> r^-1, S -> s^-1
    else:
        gens = [R, Ri, S, Si]
    words = [w for w in P["conder_relators"] if "T" not in w]
    ok = True
    for w in words:
        lw = parse_word(w)
        variants = [lw]
        if P["kind"] == "R":                         # rotation group of a reflexible map: w and w(R^-1, S^-1)
            variants.append([INV[a] for a in lw])
        for v in variants:
            ok &= all(apply_word(v, gens, x) == x for x in range(D))
    out["c_conder_relators_hold"] = ok
    out["c_conder_words_checked"] = len(words)
    # (d) code parameters
    rx, rz = len(echelon(HX)), len(echelon(HZ))
    k = n - rx - rz
    out["commute"] = all(bin(to_int(a) & to_int(b)).count("1") % 2 == 0 for a in HX for b in HZ)
    colx, colz = [0] * n, [0] * n
    for r in HX:
        for j in r:
            colx[j] += 1
    for r in HZ:
        for j in r:
            colz[j] += 1
    out.update(n=n, k=k, rankX=rx, rankZ=rz, k_is_2_minus_chi=(k == 2 - out["chi"]),
               wX=max(map(len, HX)), wZ=max(map(len, HZ)), colX=sorted(set(colx)), colZ=sorted(set(colz)))
    out["d_parameters_match"] = (n == C["n"] and k == C["k"])
    # (e) witnesses
    pivX, pivZ = echelon(HX), echelon(HZ)
    for side, ker_rows, piv in (("X", HZ, pivX), ("Z", HX, pivZ)):
        w = C["witness_" + side]
        v = to_int(w)
        okw = (len(set(w)) == len(w) == C["d_" + side]
               and all(bin(v & to_int(r)).count("1") % 2 == 0 for r in ker_rows) and reduce(piv, v) != 0)
        out["e_witness_" + side] = okw
    out["ok"] = all(out[x] for x in ("a_darts_give_matrices", "order_is_2E", "b_base_isomorphic",
                                     "b_map_isomorphic", "b_rows_equal", "c_conder_relators_hold", "commute",
                                     "k_is_2_minus_chi", "d_parameters_match", "e_witness_X", "e_witness_Z"))
    out["cpu_secs"] = round(time.process_time() - t0, 2)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("slugs", nargs="*")
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    a = ap.parse_args(argv)

    def load(name):
        with open(os.path.join(a.data, name)) as f:
            return json.load(f)
    maps, pres, params = load("maps.json"), load("presentations.json"), load("parameters.json")
    with open(os.path.join(a.data, "cover_voltage.txt")) as f:
        cover_words = [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]
    bad = 0
    for slug in a.slugs or list(params):
        try:
            r = check(slug, a.data, maps, pres, params, cover_words if "voltage" in maps[slug] else [])
        except (ValueError, KeyError, IndexError, AssertionError) as e:
            r = {"slug": slug, "error": f"{type(e).__name__}: {e}", "ok": False}
        bad += not r["ok"]
        print(json.dumps(r), flush=True)
    print(json.dumps({"codes": len(a.slugs or params), "failed": bad}))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
