#!/usr/bin/env python3
"""Independent exact cross-check by exhaustive enumeration of circuits (simple cycles) of bounded length.

usage: circuits_check.py NAME.hx NAME.hz SIDE LMAX
SIDE is Z (graph of H_X, nontrivial = not in rowsp H_Z) or X (graph of H_Z, nontrivial = not in rowsp H_X).
Enumerates every circuit of length <= LMAX of the check graph (DFS over simple paths from each start vertex
through larger vertices only, so each circuit is met from its smallest vertex, once per direction, and is
counted for one direction only) and tests each for membership in the row space of the other side.
Prints one JSON line with the number of circuits and of nontrivial circuits per length, and the least
nontrivial length (null if there is none of length <= LMAX).

Correctness (NOTES.md, Lemma 12): a minimum-weight element of ker A outside rowsp B is a circuit, so
d_side = least length of a nontrivial circuit; "no nontrivial circuit of length <= LMAX" proves d_side > LMAX,
and the number of minimum-weight logical operators of that type is the number of nontrivial circuits of
length d_side. Shares no code with graphlike_B.c (no covers, no logical basis).
"""
import json
import sys
import time


def read_checks(path):
    with open(path) as f:
        lines = f.read().split("\n")
    n, m = map(int, lines[0].split())
    return n, [[int(t) for t in lines[1 + i].split()] for i in range(m)]


def rref(rows):
    """Fully reduced echelon basis {pivot: vector}: vector has bit pivot and no other pivot bit."""
    piv = {}
    for v in rows:
        for p in sorted(piv, reverse=True):
            if (v >> p) & 1:
                v ^= piv[p]
        if v:
            h = v.bit_length() - 1
            for p in piv:
                if (piv[p] >> h) & 1:
                    piv[p] ^= v
            piv[h] = v
    return piv


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    t0 = time.time()
    n, HX = read_checks(sys.argv[1])
    _, HZ = read_checks(sys.argv[2])
    side, lmax = sys.argv[3], int(sys.argv[4])
    A, B = (HX, HZ) if side == "Z" else (HZ, HX)
    cols = [[] for _ in range(n)]
    for i, r in enumerate(A):
        for c in r:
            cols[c].append(i)
    if any(len(c) > 2 for c in cols):
        print(json.dumps({"side": side, "error": "not graph-like"}))
        return
    inf = len(A)
    adj = [[] for _ in range(inf + 1)]
    loops = []
    for e, c in enumerate(cols):
        u, w = (c + [inf, inf])[:2]
        if u == w:
            loops.append(e)
        else:
            adj[u].append((w, e))
            adj[w].append((u, e))
    piv = rref([sum(1 << c for c in r) for r in B])

    def trivial(edges, mask):
        # in the span iff mask equals the sum of the basis vectors at its pivot positions
        w = 0
        for e in edges:
            b = piv.get(e)
            if b is not None:
                w ^= b
        return w == mask

    count = {}
    nontriv = {}
    for e in loops:
        count[1] = count.get(1, 0) + 1
        if not trivial([e], 1 << e):
            nontriv[1] = nontriv.get(1, 0) + 1

    path = []

    def dfs(s, u, mask, first_edge):
        depth = len(path)
        for w, e in adj[u]:
            if (mask >> e) & 1:
                continue
            if w == s:
                if e > first_edge:  # each circuit is met in two directions; count one
                    L = depth + 1
                    count[L] = count.get(L, 0) + 1
                    path.append(e)
                    if not trivial(path, mask | (1 << e)):
                        nontriv[L] = nontriv.get(L, 0) + 1
                    path.pop()
                continue
            if w < s or depth + 1 >= lmax or w in onpath:
                continue
            onpath.add(w)
            path.append(e)
            dfs(s, w, mask | (1 << e), first_edge)
            path.pop()
            onpath.discard(w)

    if lmax >= 2:
        for s in range(inf + 1):
            for w, e in adj[s]:
                if w > s:
                    onpath = {s, w}
                    path.append(e)
                    dfs(s, w, 1 << e, e)
                    path.pop()
    lengths = sorted(count)
    out = {"side": side, "lmax": lmax,
           "circuits": {L: count[L] for L in lengths},
           "nontrivial": {L: nontriv.get(L, 0) for L in lengths},
           "d": min(nontriv) if nontriv else None,
           "secs": round(time.time() - t0, 2)}
    print(json.dumps(out))


if __name__ == "__main__":
    main()
