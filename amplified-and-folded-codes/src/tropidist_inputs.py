#!/usr/bin/env python3
"""Inputs for tropidist (github.com/Fanerst/tropidist) in the encodings of the board's certifier, and the recheck
of a tropidist record.

    python3 src/tropidist_inputs.py write OUTDIR SLUG.SIDE ...   -> OUTDIR/SLUG.SIDE.H.mtx, OUTDIR/SLUG.SIDE.L.mtx
    python3 src/tropidist_inputs.py check SLUG.SIDE RECORD       -> recheck the vector of a tropidist record

tropidist computes the least weight of a vector e with H e = 0 and L e != 0. Stabilizer code, side P: three
columns X_i, Z_i, Y_i per qubit (blocks of n, in this order); a row of H is the symplectic product of a generator
(x | z) with every column, i.e. (z | x | x + z), and L is the same for a basis of N(S)/S; a minimum never uses two
columns of one qubit, so the least weight is the Pauli distance. CSS code, side X: H = H_Z and L a basis of
ker H_X modulo rowsp H_Z; side Z mirrored. Run from the directory amplified-and-folded-codes/.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rebuild as R  # noqa: E402


def nullspace(rows, ncols):
    """Basis of {v : r.v = 0 for every row r}, vectors as integers with ncols bits."""
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
    basis = []
    for f in range(ncols):
        if f in piv:
            continue
        v = 1 << f
        for c, p in piv.items():
            if p >> f & 1:
                v |= 1 << c
        basis.append(v)
    return basis


def problem(slug, side):
    c = R.load("data", "codes", f"{slug}.json")
    n = c["n"]
    if c["code_type"] == "stabilizer":
        if side != "P":
            raise ValueError("a stabilizer code has the single side P")
        G = [R.pvec(g, n) for g in c["checks"]["S"]]
        mask = (1 << n) - 1
        swapped = [((g >> n) & mask) | ((g & mask) << n) for g in G]
        span = R.Span(G)
        L = [v for v in nullspace(swapped, 2 * n) if span.add(v)]
        assert len(L) == 2 * c["k"], len(L)

        def mech(v):
            x, z = v & mask, v >> n
            return z | (x << n) | ((x ^ z) << (2 * n))
        return 3 * n, [mech(g) for g in G], [mech(v) for v in L]
    X = [R.vec(r) for r in c["checks"]["X"]]
    Z = [R.vec(r) for r in c["checks"]["Z"]]
    same, opp = (X, Z) if side == "X" else (Z, X)
    span = R.Span(opp)
    L = [v for v in nullspace(same, n) if span.add(v)]
    assert len(L) == c["k"], len(L)
    return n, opp, L


def write_mtx(path, rows, ncols):
    ent = [(i + 1, j + 1) for i, r in enumerate(rows) for j in R.support(r)]
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("%%MatrixMarket matrix coordinate pattern general\n")
        f.write(f"{len(rows)} {ncols} {len(ent)}\n")
        for i, j in ent:
            f.write(f"{i} {j}\n")


def check(job, record):
    slug, side = job.rsplit(".", 1)
    c = R.load("data", "codes", f"{slug}.json")
    n = c["n"]
    with open(record, encoding="utf-8") as f:
        rec = json.load(f)
    cols = rec.get("witness") or []
    if not cols:
        print(f"{job}: status {rec.get('status')}, lower bound {rec.get('lower_bound')}, no vector returned, "
              f"{rec.get('seconds')} s")
        return 0 if rec.get("status") != "exact" else 1
    if c["code_type"] == "stabilizer":
        x = z = 0
        for col in cols:
            kind, q = divmod(int(col), n)
            if kind in (0, 2):
                x ^= 1 << q
            if kind in (1, 2):
                z ^= 1 << q
        v = x | (z << n)
        G = [R.pvec(g, n) for g in c["checks"]["S"]]
        ok = all(R.sympl(v, g, n) == 0 for g in G) and v not in R.Span(G)
        weight = R.pauli_weight(v, n)
    else:
        v = 0
        for col in cols:
            v ^= 1 << int(col)
        X = [R.vec(r) for r in c["checks"]["X"]]
        Z = [R.vec(r) for r in c["checks"]["Z"]]
        same, opp = (X, Z) if side == "X" else (Z, X)
        ok = all(R.popcount(v & r) % 2 == 0 for r in opp) and v not in R.Span(same)
        weight = R.popcount(v)
    good = ok and weight == rec.get("upper_bound")
    print(f"{job}: status {rec.get('status')}, lower bound {rec.get('lower_bound')}, upper bound "
          f"{rec.get('upper_bound')}, {rec.get('seconds')} s; returned vector is a nontrivial logical of weight "
          f"{weight}: {good}")
    return 0 if good else 1


def main():
    if len(sys.argv) >= 4 and sys.argv[1] == "write":
        os.makedirs(sys.argv[2], exist_ok=True)
        for job in sys.argv[3:]:
            slug, side = job.rsplit(".", 1)
            ncols, H, L = problem(slug, side)
            write_mtx(os.path.join(sys.argv[2], f"{job}.H.mtx"), H, ncols)
            write_mtx(os.path.join(sys.argv[2], f"{job}.L.mtx"), L, ncols)
            print(f"{job}: H {len(H)} x {ncols}, L {len(L)} x {ncols}")
        return 0
    if len(sys.argv) == 4 and sys.argv[1] == "check":
        return check(sys.argv[2], sys.argv[3])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
