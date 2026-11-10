#!/usr/bin/env python3
"""Build the five codes of exact_distances.tex and the two validation codes, check them, and write the
search instances.

    python3 src/codes.py [--board DIR]

Writes into instances/:
  vinst_<name>_<T>.txt  input of vsearch_orb (role-ordered weight-6 checks, majority family, dual basis),
  cinst_<name>_<T>.txt  input of cluster2 (header n mX mZ k N a b c, then H_X, H_Z and the dual basis),
  jobs_C600.txt         anchors with exclusions for vsearch_orb -L (orbits of the automorphism group of C600),
where T = Z is the search for Z-logical operators (checks H_X, dual basis = X-logical operators) and T = X the
mirror image. Checks, for every code: H_X H_Z^T = 0, n, k, check weights, odd column weights (Lemma 2.2), the
X/Z exchange of Lemma 2.3 (cyclic codes), the automorphism generators of C600 (data/C600_automorphisms.json),
and the witnesses in data/. With --board DIR (the codes/ directory of the qLDPC challenge repository) the
matrices of the cyclic board codes are compared with the board's files.

Requires Python 3 only (standard library).
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent.parent
DATA = HERE / "data"
INST = HERE / "instances"


# ---------------------------------------------------------------- GF(2) with Python ints
def mask(supp):
    m = 0
    for q in supp:
        m |= 1 << q
    return m


def bits(v):
    out, i = [], 0
    while v:
        if v & 1:
            out.append(i)
        v >>= 1
        i += 1
    return out


def rref(rows):
    piv = {}
    for r in rows:
        for p, pr in piv.items():
            if r >> p & 1:
                r ^= pr
        if r:
            p = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> p & 1:
                    piv[q] ^= r
            piv[p] = r
    return piv


def reduce(v, piv):
    for p, pr in piv.items():
        if v >> p & 1:
            v ^= pr
    return v


def kernel(rows, n):
    piv = rref(rows)
    pivcols = set(piv)
    basis = []
    for f in range(n):
        if f in pivcols:
            continue
        v = 1 << f
        for p, pr in piv.items():
            if pr >> f & 1:
                v |= 1 << p
        basis.append(v)
    return basis


def dual_basis(comm_rows, stab_rows, n):
    """Basis of ker(comm_rows) modulo rowsp(stab_rows)."""
    piv = rref(stab_rows)
    out = []
    for v in kernel(comm_rows, n):
        r = reduce(v, piv)
        if r:
            out.append(v)
            p = r.bit_length() - 1
            for q in list(piv):
                if piv[q] >> p & 1:
                    piv[q] ^= r
            piv[p] = r
    return out


def parity(x):
    return bin(x).count("1") & 1


# ---------------------------------------------------------------- constructions (Section 1 of the note)
def twobga(m, q, r, A, B):
    """C(G;a,b) for G = <x,y | x^m, y^q, y x y^-1 = x^r>, x^i y^j at index q*i + j."""
    els = [(i, j) for i in range(m) for j in range(q)]
    idx = {e: q * e[0] + e[1] for e in els}

    def mul(g, h):
        (i1, j1), (i2, j2) = g, h
        return ((i1 + pow(r, j1, m) * i2) % m, (j1 + j2) % q)

    inv = {}
    for g in els:
        for h in els:
            if mul(g, h) == (0, 0):
                inv[g] = h
                break
    N = len(els)
    A = [(i % m, j % q) for i, j in A]
    B = [(i % m, j % q) for i, j in B]
    HX = [[idx[mul(inv[s], c)] for s in A] + [N + idx[mul(c, inv[t])] for t in B] for c in els]
    HZ = [[idx[mul(c, t)] for t in B] + [N + idx[mul(s, c)] for s in A] for c in els]
    return HX, HZ, None, (m, q, r)


def cyclic(m, a, b):
    """H_X = [circ(a) | circ(b)], H_Z = [circ(b)^T | circ(a)^T], circ(v)[i, (i + j) mod m] = v[j]."""
    HX = [[(i + s) % m for s in a] + [m + (i + t) % m for t in b] for i in range(m)]
    HZ = [[(i - t) % m for t in b] + [m + (i - s) % m for s in a] for i in range(m)]
    perm = [m + (-i) % m for i in range(m)] + [(-i) % m for i in range(m)]  # Lemma 2.3
    return HX, HZ, perm, (m, 0, 1)


def abelian(ell, m, a, b):
    """Bivariate version of cyclic() on Z/ell x Z/m, (i, j) at index m*i + j; a, b lists of exponents (u, v)."""
    N = ell * m
    idx = lambda i, j: (i % ell) * m + (j % m)
    HX = [[idx(i + u, j + v) for u, v in a] + [N + idx(i + u, j + v) for u, v in b] for i in range(ell) for j in range(m)]
    HZ = [[idx(i - u, j - v) for u, v in b] + [N + idx(i - u, j - v) for u, v in a] for i in range(ell) for j in range(m)]
    perm = [N + idx(-(g // m), -(g % m)) for g in range(N)] + [idx(-(g // m), -(g % m)) for g in range(N)]
    return HX, HZ, perm, (ell, 0, m)


# name: (construction, parameters, board slug, claimed distance); V144 and V288 are the validation codes
# (gross and two-gross code of Bravyi et al., A = x^3 + y + y^2, B = y^3 + x + x^2 resp. A = x^3 + y^2 + y^7)
CODES = {
    "V144": ("abelian", (12, 6, [(3, 0), (0, 1), (0, 2)], [(0, 3), (1, 0), (2, 0)]), None, 12),
    "V288": ("abelian", (12, 12, [(3, 0), (0, 2), (0, 7)], [(0, 3), (1, 0), (2, 0)]), None, 18),
    "C600": ("twobga", (30, 10, 19, [(0, 0), (24, 8), (13, 9)], [(0, 0), (22, 4), (24, 7)]), None, 26),
    "C504": ("cyclic", (252, [0, 1, 58], [0, 48, 166]), None, 28),
    "B510": ("cyclic", (255, [0, 26, 76], [0, 30, 66]), "510-16-26", 26),
    "B558": ("cyclic", (279, [0, 37, 140], [0, 61, 197]), "558-14-28", 28),
    "B434": ("cyclic", (217, [0, 29, 117], [0, 19, 111]), "434-10-28", 28),
}


def majority(stabs):
    sm = [mask(r) for r in stabs]
    fam = [sorted(r) for r in stabs]
    seen = {tuple(f) for f in fam}
    for i in range(len(sm)):
        for j in range(i + 1, len(sm)):
            if sm[i] & sm[j]:
                s = tuple(bits(sm[i] ^ sm[j]))
                if s and s not in seen:
                    seen.add(s)
                    fam.append(list(s))
    return fam


def write_instances(name, T, checks, stabs, dual, n, N, k, lat):
    fam = majority(stabs)
    with open(INST / f"vinst_{name}_{T}.txt", "w") as fh:
        fh.write(f"{n} {N} {k} {len(fam)}\n")
        for r in checks:
            fh.write(" ".join(map(str, r)) + "\n")
        for s in fam:
            fh.write(f"{len(s)} " + " ".join(map(str, s)) + "\n")
        for v in dual:
            b = bits(v)
            fh.write(f"{len(b)} " + " ".join(map(str, b)) + "\n")
    with open(INST / f"cinst_{name}_{T}.txt", "w") as fh:
        fh.write(f"{n} {len(checks)} {len(stabs)} {k} {N} {lat[0]} {lat[1]} {lat[2]}\n")
        for M in (checks, stabs):
            for r in M:
                fh.write(f"{len(r)} " + " ".join(map(str, sorted(r))) + "\n")
        for v in dual:
            b = bits(v)
            fh.write(f"{len(b)} " + " ".join(map(str, b)) + "\n")


def orbits(gens, n):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for g in gens:
        for x in range(n):
            a, b = find(x), find(g[x])
            if a != b:
                parent[max(a, b)] = min(a, b)
    orb = {}
    for x in range(n):
        orb.setdefault(find(x), []).append(x)
    return sorted(orb.values(), key=lambda o: o[0])


def main():
    board = None
    if "--board" in sys.argv:
        board = Path(sys.argv[sys.argv.index("--board") + 1])
    INST.mkdir(exist_ok=True)
    ok = True
    for name, (kind, args, slug, dclaim) in CODES.items():
        HX, HZ, perm, lat = {"twobga": twobga, "cyclic": cyclic, "abelian": abelian}[kind](*args)
        N = len(HX)
        n = 2 * N
        mX = [mask(r) for r in HX]
        mZ = [mask(r) for r in HZ]
        css = all(parity(x & z) == 0 for x in mX for z in mZ)
        rX, rZ = len(rref(mX)), len(rref(mZ))
        k = n - rX - rZ
        LX = dual_basis(mZ, mX, n)
        LZ = dual_basis(mX, mZ, n)
        weights = sorted({len(set(r)) for r in HX + HZ})
        colX = [0] * n
        colZ = [0] * n
        for r in HX:
            for q in r:
                colX[q] += 1
        for r in HZ:
            for q in r:
                colZ[q] += 1
        odd = all(c % 2 for c in colX) and all(c % 2 for c in colZ)
        rowsX = {tuple(sorted(r)) for r in HX}
        rowsZ = {tuple(sorted(r)) for r in HZ}
        line = f"{name}: n={n} k={k} CSS={css} check weights={weights} odd column weights={odd}"
        if perm is not None:
            xz = {tuple(sorted(perm[q] for q in r)) for r in HX} == rowsZ and \
                 {tuple(sorted(perm[q] for q in r)) for r in HZ} == rowsX
            line += f" X/Z exchange (Lemma 2.3)={xz}"
            ok &= xz
        kclaim = {"V144": 12, "V288": 12, "C600": 24, "C504": 12, "B510": 16, "B558": 14, "B434": 10}[name]
        ok &= css and odd and weights == [6] and len(LX) == k == len(LZ) == kclaim
        json.dump({"name": name, "construction": kind, "parameters": args, "n": n, "k": k,
                   "HX": [sorted(r) for r in HX], "HZ": [sorted(r) for r in HZ]},
                  open(DATA / f"{name}.json", "w"))
        write_instances(name, "Z", HX, HZ, LX, n, N, k, lat)
        if perm is None:
            write_instances(name, "X", HZ, HX, LZ, n, N, k, lat)
        if name == "C600":
            aut = json.load(open(DATA / "C600_automorphisms.json"))
            good = all({tuple(sorted(g[q] for q in r)) for r in HX} == rowsX and
                       {tuple(sorted(g[q] for q in r)) for r in HZ} == rowsZ and
                       sorted(g) == list(range(n)) for g in aut["generators"])
            orb = orbits(aut["generators"], n)
            blocks_ok = all(all(q < N for q in o) or all(q >= N for q in o) for o in orb)
            with open(INST / "jobs_C600.txt", "w") as fh:
                excl = []
                for o in orb:
                    fh.write(f"{o[0]} {len(excl)} " + " ".join(map(str, excl)) + "\n")
                    excl += o
            line += f" automorphisms ok={good} orbits={[len(o) for o in orb]} block-preserving={blocks_ok}"
            ok &= good and blocks_ok
        for T, comm, dual in (("Z", mX, LX), ("X", mZ, LZ)):
            wf = DATA / f"{name}_witnesses_{T}.json"
            if not wf.exists():
                continue
            words = json.load(open(wf))
            valid = [z for z in words if all(parity(mask(z) & c) == 0 for c in comm)
                     and any(parity(mask(z) & l) for l in dual) and len(set(z)) == dclaim]
            line += f" {T}-witnesses {len(valid)}/{len(words)} valid of weight {dclaim}"
            ok &= len(valid) == len(words) > 0
        if board is not None and slug is not None:
            d = json.load(open(board / f"{slug}.json"))
            same = {tuple(sorted(r)) for r in d["checks"]["X"]} == rowsX and \
                   {tuple(sorted(r)) for r in d["checks"]["Z"]} == rowsZ
            line += f" equal to board {slug}: {same}"
            ok &= same
        print(line, flush=True)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
