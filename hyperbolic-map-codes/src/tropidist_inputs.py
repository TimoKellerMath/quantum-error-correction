#!/usr/bin/env python3
"""Write the tropidist problems of the nine codes in the convention of the qLDPC challenge board.

    python3 src/tropidist_inputs.py OUTDIR [SLUG ...]

For each code and side, OUTDIR/<slug>.<side>.H.mtx and OUTDIR/<slug>.<side>.L.mtx (MatrixMarket coordinate
pattern, 1-based, entries in row-major order). Side X asks for the least weight of an X-type logical operator:
H = H_Z and L = representatives of ker H_X modulo rowsp H_Z (Z-type logical operators), so that e in ker H_Z
is outside rowsp H_X exactly when L e != 0. Side Z mirrors it. This is the "css-sides" encoding of the board's
tropidist backend (research/certs/tropidist_certify.py in github.com/unitaryfoundation/qldpc-challenge).
Any basis of the logical operators gives the same problem. Standard library only.

tropidist (github.com/Fanerst/tropidist, MIT license, commit b03ec2d4077f) is then run from its checkout as

    python -m cdist --matrices OUTDIR/<slug>.<side>.H.mtx OUTDIR/<slug>.<side>.L.mtx \\
        --threads 1 --time-limit 600 --out <slug>.<side>.json

Our records are in logs/tropidist/.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def to_int(row):
    v = 0
    for j in row:
        v ^= 1 << j
    return v


def to_list(v):
    out, j = [], 0
    while v:
        if v & 1:
            out.append(j)
        v >>= 1
        j += 1
    return out


def kernel_basis(rows, n):
    piv = {}
    for r in rows:
        v = to_int(r)
        for c, p in piv.items():
            if (v >> c) & 1:
                v ^= p
        if not v:
            continue
        c = (v & -v).bit_length() - 1
        for c2 in list(piv):
            if (piv[c2] >> c) & 1:
                piv[c2] ^= v
        piv[c] = v
    basis = []
    for f in range(n):
        if f in piv:
            continue
        x = 1 << f
        for c, p in piv.items():
            if (p >> f) & 1:
                x |= 1 << c
        basis.append(x)
    return basis


def logical_basis(H_same, H_opp, n):
    """ker(H_same) modulo rowsp(H_opp)."""
    piv = {}

    def add(v):
        while v:
            h = v.bit_length() - 1
            if h not in piv:
                piv[h] = v
                return True
            v ^= piv[h]
        return False
    for r in H_opp:
        add(to_int(r))
    return [x for x in kernel_basis(H_same, n) if add(x)]


def read_rows(path):
    with open(path) as f:
        n, m = map(int, f.readline().split())
        return n, [sorted(map(int, f.readline().split())) for _ in range(m)]


def write_mtx(path, rows, ncols):
    ent = [(i, j) for i, r in enumerate(rows) for j in sorted(r)]
    with open(path, "w") as f:
        f.write("%%MatrixMarket matrix coordinate pattern general\n")
        f.write(f"{len(rows)} {ncols} {len(ent)}\n")
        for i, j in ent:
            f.write(f"{i + 1} {j + 1}\n")


def main(argv):
    if not argv:
        sys.exit(__doc__)
    out = argv[0]
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(ROOT, "data", "parameters.json")) as f:
        params = json.load(f)
    for slug in argv[1:] or list(params):
        n, HX = read_rows(os.path.join(ROOT, "data", "codes", slug + ".hx"))
        _, HZ = read_rows(os.path.join(ROOT, "data", "codes", slug + ".hz"))
        for side, H_same, H_opp in (("X", HX, HZ), ("Z", HZ, HX)):
            L = logical_basis(H_same, H_opp, n)
            assert len(L) == params[slug]["k"], (slug, side, len(L))
            write_mtx(os.path.join(out, f"{slug}.{side}.H.mtx"), H_opp, n)
            write_mtx(os.path.join(out, f"{slug}.{side}.L.mtx"), [to_list(v) for v in L], n)
        print(slug, "k =", params[slug]["k"], "written")


if __name__ == "__main__":
    main(sys.argv[1:])
