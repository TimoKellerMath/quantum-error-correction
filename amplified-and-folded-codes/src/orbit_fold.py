#!/usr/bin/env python3
"""Orbit fold of a stabilizer code along a cyclic group H = <pi> of qubit automorphisms (emulates the
board gate's orbit-fold pass, exactly).  A Pauli p constant on the H-orbits is the lift of a Pauli q on
the orbit set O; for every Pauli g, <p, g> = <q, fold(g)>, fold(g)_o = sum of g over the orbit o.
Hence p in N(S) iff q is orthogonal to fold(S), and p in S (given p in N(S)) iff <q, fold(L_l)> = 0
for all logical representatives L_l.  Output: a .code file for mitm.c / fewnony2.c-style tools with
n' = |O|, an independent subset of the folded generators, and the folded logicals (labels only).
All orbits must have the same size m (true for fixed-point-free pi); Pauli weight = m * orbit weight.

Usage: python3 -I orbit_fold.py CODE.json SHIFT OUT.code   (pi: h -> h + SHIFT mod n)"""
import json, sys


def load(path):
    d = json.load(open(path))
    n = d["n"]
    def vec(g):
        return sum(1 << q for q in g["X"]) | sum(1 << (n + q) for q in g["Z"])
    return n, [vec(g) for g in d["S"]], [vec(g) for g in d["logicals"]]


def main():
    n, S, Lg = load(sys.argv[1])
    shift = int(sys.argv[2])
    perm = [(h + shift) % n for h in range(n)]
    # orbits
    seen, orbits = set(), []
    for h in range(n):
        if h in seen:
            continue
        o, x = [], h
        while x not in o:
            o.append(x)
            x = perm[x]
        seen |= set(o)
        orbits.append(o)
    sizes = {len(o) for o in orbits}
    assert len(sizes) == 1, sizes
    m = sizes.pop()
    no = len(orbits)
    orb_of = {}
    for i, o in enumerate(orbits):
        for q in o:
            orb_of[q] = i
    mask = (1 << n) - 1

    def permv(v):
        x, z = v & mask, v >> n
        w = 0
        for q in range(n):
            if x >> q & 1:
                w |= 1 << perm[q]
            if z >> q & 1:
                w |= 1 << (n + perm[q])
        return w
    # pi must preserve the stabilizer group (row space)
    def rank(rows):
        piv, r = {}, 0
        for v in rows:
            while v:
                b = v.bit_length() - 1
                if b in piv:
                    v ^= piv[b]
                else:
                    piv[b] = v; r += 1; break
        return r
    rS = rank(S)
    assert rank(S + [permv(s) for s in S]) == rS, "shift is not an automorphism of S"

    def fold(v):
        x, z = v & mask, v >> n
        fx = fz = 0
        for q in range(n):
            if x >> q & 1:
                fx ^= 1 << orb_of[q]
            if z >> q & 1:
                fz ^= 1 << orb_of[q]
        return fx | (fz << no)
    FS = [fold(s) for s in S]
    # independent subset of the folded generators (same span)
    piv, indep = {}, []
    for v in FS:
        w = v
        while w:
            b = w.bit_length() - 1
            if b in piv:
                w ^= piv[b]
            else:
                piv[b] = w; indep.append(v); break
    FL = [fold(l) for l in Lg]
    rho = len(indep)
    print(f"n {n}, shift {shift}: {no} orbits of size {m}; rank S {rS}; rank of folded generators {rho}; "
          f"dim of H-invariant normalizer = 2*{no} - {rho} = {2 * no - rho}")
    with open(sys.argv[3], "w") as f:
        f.write(f"{no} {rho} {len(FL)}\n")
        for v in indep + FL:
            f.write("".join("1" if v >> i & 1 else "0" for i in range(2 * no)) + "\n")


if __name__ == "__main__":
    main()
