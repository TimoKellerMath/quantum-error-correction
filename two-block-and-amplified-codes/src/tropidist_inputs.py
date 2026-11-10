#!/usr/bin/env python3
"""Inputs for tropidist (github.com/Fanerst/tropidist) in the encoding of the board's certifier, and the recheck
of a tropidist record.

    python3 src/tropidist_inputs.py write OUTDIR SLUG.SIDE ...   -> OUTDIR/SLUG.SIDE.H.mtx, OUTDIR/SLUG.SIDE.L.mtx
    python3 src/tropidist_inputs.py check SLUG.SIDE RECORD       -> recheck the vector of a tropidist record

tropidist computes the least weight of a vector e with H e = 0 and L e != 0. Side Z (d_Z): H = H_X and L a basis
of the X-logicals (ker H_Z modulo rowsp H_X); side X mirrored. The basis is the one of the echelon procedure below
(the procedure of our first runs on 480-16-26, so that the inputs of all runs are reproduced byte for byte).
Run from the directory two-block-and-amplified-codes/.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return json.load(f)


def mask(row):
    m = 0
    for q in row:
        m |= 1 << q
    return m


def echelon(rows):
    piv = {}
    for v in rows:
        for p in sorted(piv, reverse=True):
            if v >> p & 1:
                v ^= piv[p]
        if v:
            p = v.bit_length() - 1
            for p2 in list(piv):
                if piv[p2] >> p & 1:
                    piv[p2] ^= v
            piv[p] = v
    return piv


def reduce(v, piv):
    for p in sorted(piv, reverse=True):
        if v >> p & 1:
            v ^= piv[p]
    return v


def kernel(rows, n):
    piv = echelon(rows)
    out = []
    for f in (c for c in range(n) if c not in piv):
        v = 1 << f
        for p, rw in piv.items():
            if rw >> f & 1:
                v |= 1 << p
        out.append(v)
    return out


def logicals(comm, stab, n):
    """Representatives of ker(comm) / rowsp(stab)."""
    piv = echelon(stab)
    out = []
    for v in kernel(comm, n):
        w = reduce(v, piv)
        if w:
            out.append(v)
            p = w.bit_length() - 1
            for p2 in list(piv):
                if piv[p2] >> p & 1:
                    piv[p2] ^= w
            piv[p] = w
    return out


def problem(slug, side):
    c = load("data", "codes", f"{slug}.json")
    n = c["n"]
    MX = [mask(r) for r in c["checks"]["X"]]
    MZ = [mask(r) for r in c["checks"]["Z"]]
    if side == "Z":
        H, L = MX, logicals(MZ, MX, n)
    elif side == "X":
        H, L = MZ, logicals(MX, MZ, n)
    else:
        raise ValueError("side must be X or Z")
    assert len(L) == c["k"], (slug, side, len(L))
    return n, H, L


def write_mtx(path, rows, n):
    ent = [(i + 1, q + 1) for i, v in enumerate(rows) for q in range(n) if v >> q & 1]
    with open(path, "w", encoding="ascii", newline="\n") as fh:
        fh.write("%%MatrixMarket matrix coordinate pattern general\n")
        fh.write(f"{len(rows)} {n} {len(ent)}\n")
        for i, q in ent:
            fh.write(f"{i} {q}\n")


def check(job, record):
    slug, side = job.rsplit(".", 1)
    c = load("data", "codes", f"{slug}.json")
    with open(record, encoding="utf-8") as f:
        rec = json.load(f)
    cols = rec.get("witness") or []
    if not cols:
        print(f"{job}: status {rec.get('status')}, lower bound {rec.get('lower_bound')}, no vector returned, "
              f"{rec.get('seconds')} s")
        return 0 if rec.get("status") != "exact" else 1
    v = 0
    for col in cols:
        v ^= 1 << int(col)
    MX = [mask(r) for r in c["checks"]["X"]]
    MZ = [mask(r) for r in c["checks"]["Z"]]
    comm, stab = (MX, MZ) if side == "Z" else (MZ, MX)
    ok = all(bin(v & r).count("1") % 2 == 0 for r in comm) and reduce(v, echelon(stab)) != 0
    weight = bin(v).count("1")
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
            n, H, L = problem(slug, side)
            write_mtx(os.path.join(sys.argv[2], f"{job}.H.mtx"), H, n)
            write_mtx(os.path.join(sys.argv[2], f"{job}.L.mtx"), L, n)
            print(f"{job}: H {len(H)} x {n}, L {len(L)} x {n}")
        return 0
    if len(sys.argv) == 4 and sys.argv[1] == "check":
        return check(sys.argv[2], sys.argv[3])
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
