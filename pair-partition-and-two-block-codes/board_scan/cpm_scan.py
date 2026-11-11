#!/usr/bin/env python3
"""Which CSS board entries of check weight <= 6 are quasi-cyclic CPM codes in a natural qubit layout?

    cpm_scan.py CODESDIR [PMIN]

For every factorization n = H * M * P (P >= PMIN, default 3) the qubit q = h*M*P + s*M + c is position s in Z_P of
block column (h, c) (contiguous blocks: M = 1; interleaved blocks: H = 1; per-half interleaving: H = 2). The code is a
CPM code over Z_P in that layout when (1) no check has two qubits in one block column and (2) the shift
s -> s + 1 (mod P) maps the X-checks onto X-checks and the Z-checks onto Z-checks. Then every block of H_X and H_Z is
0 or one circulant permutation matrix, H_X H_Z^T = 0 is the paired-difference condition (every exponent difference
on an X-row/Z-row overlap occurs an even number of times), and the template has J_X, J_Z rows (check orbits) on
L = H*M columns. Reported per entry: the layouts found (largest P first), the template shape, whether it is all-one,
the column degrees, and whether two checks of one type share two qubits (the girth rule fails then).
"""
import json
import os
import sys
from collections import Counter

CODES, PMIN = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 3


def factorizations(n):
    for P in range(PMIN, n // 2 + 1):
        if n % P:
            continue
        L = n // P
        for H in range(1, L + 1):
            if L % H == 0:
                yield H, L // H, P


def test(rows_x, rows_z, n, H, M, P):
    def split(q):
        h, r = divmod(q, M * P)
        s, c = divmod(r, M)
        return h * M + c, s
    for r in rows_x + rows_z:
        bl = [split(q)[0] for q in r]
        if len(set(bl)) != len(bl):
            return None
    def shift(q):
        h, r = divmod(q, M * P)
        s, c = divmod(r, M)
        return h * M * P + ((s + 1) % P) * M + c
    out = {}
    for T, rows in (("X", rows_x), ("Z", rows_z)):
        S = {frozenset(r) for r in rows}
        if any(frozenset(shift(q) for q in r) not in S for r in S):
            return None
        # orbits of checks under the shift
        seen, orbits = set(), []
        for r in S:
            if r in seen:
                continue
            orb, x = [], r
            while x not in seen:
                seen.add(x)
                orb.append(x)
                x = frozenset(shift(q) for q in x)
            orbits.append(orb)
        tmpl = [sorted({split(q)[0] for q in orb[0]}) for orb in orbits]
        out[T] = {"J": len(orbits), "orbit_sizes": sorted(Counter(len(o) for o in orbits).items()), "rows": tmpl}
    L = H * M
    for T in "XZ":
        deg = Counter(u for row in out[T]["rows"] for u in row)
        out[T]["coldeg"] = sorted(Counter(deg[u] for u in range(L)).items())
        out[T]["all_one"] = all(len(row) == L for row in out[T]["rows"])
    return out


def shares_two(rows):
    cnt = Counter()
    for r in rows:
        rs = sorted(r)
        for i in range(len(rs)):
            for j in range(i + 1, len(rs)):
                cnt[(rs[i], rs[j])] += 1
    return any(v > 1 for v in cnt.values())


res = []
for fn in sorted(os.listdir(CODES)):
    if not fn.endswith(".json"):
        continue
    d = json.load(open(os.path.join(CODES, fn)))
    if d.get("code_type") != "CSS":
        continue
    X, Z = d["checks"]["X"], d["checks"]["Z"]
    w = max(len(r) for r in X + Z)
    if w > 6:
        continue
    n = d["n"]
    found = []
    for H, M, P in factorizations(n):
        t = test(X, Z, n, H, M, P)
        if t:
            found.append((P, H, M, t))
    if not found:
        continue
    found.sort(key=lambda f: -f[0])
    P, H, M, t = found[0]
    rec = {"slug": fn[:-5], "n": n, "k": d["k"], "d": d["distance"]["d"], "w": w, "family": d.get("family"),
           "authors": d["provenance"].get("authors"), "P": P, "H": H, "M": M, "L": H * M,
           "JX": t["X"]["J"], "JZ": t["Z"]["J"], "all_one": t["X"]["all_one"] and t["Z"]["all_one"],
           "coldegX": t["X"]["coldeg"], "coldegZ": t["Z"]["coldeg"],
           "orbit_sizes": [t["X"]["orbit_sizes"], t["Z"]["orbit_sizes"]],
           "girth_rule": not shares_two(X) and not shares_two(Z),
           "layouts": [(f[0], f[1], f[2]) for f in found][:8],
           "name": d["name"][:80]}
    res.append(rec)
    print(json.dumps(rec))
print("entries with a CPM layout:", len(res), file=sys.stderr)
