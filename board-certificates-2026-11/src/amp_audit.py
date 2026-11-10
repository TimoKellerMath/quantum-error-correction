#!/usr/bin/env python3
"""Is a board entry E the amplified code Amp_N(S_X, S_Z) of a presentation (S_X, S_Z) of a board base B?

    python3 src/amp_audit.py BOARD_DIR ENTRY BASE N

Definition 1 of the README (after Liang, Gu, Chen, Eisert, Wang, arXiv:2609.37231, Eqs. (14), (15), with the
[[N,N-2,2]] amplifier): central qubits (q, a), q < n_B, a < N, numbered q*N + a ("qN+a") or a*n_B + q ("a*n+q"),
followed by side qubits. The audit
  1. reads, from the X-rows of E whose central part is {(q, a) : q in s} for a single a and which contain exactly
     one side qubit, the supports s of S_X with their side qubits s_c (and S_Z with t_z from the Z-rows);
  2. checks rowsp S_X = rowsp H_X(B) and rowsp S_Z = rowsp H_Z(B) (ranks of S, H and [S; H]);
  3. rebuilds Amp_N(S_X, S_Z) from Definition 1 with the side qubits found in E and compares its X- and Z-rows
     with those of E as sets of supports (identical rows), or else their row spaces.
If it passes, E and Amp_N(S_X, S_Z) have the same row spaces, so d_X(E) = 2 d_X(B) and d_Z(E) = 2 d_Z(B) by the
amplification proposition of the README (Proposition 3 there). Standard library only. Adapted from the audit
written for the track that found these certificates (own code, not the construction program of the codes).
"""
import json
import os
import sys


def mask(s):
    m = 0
    for q in s:
        m |= 1 << q
    return m


def rank(rows):
    piv = {}
    for r in rows:
        while r:
            p = r.bit_length() - 1
            if p in piv:
                r ^= piv[p]
            else:
                piv[p] = r
                break
    return len(piv)


def try_conv(E, B, N, conv):
    nE, XE, ZE = E["n"], E["checks"]["X"], E["checks"]["Z"]
    nB, HX, HZ = B["n"], B["checks"]["X"], B["checks"]["Z"]
    ncen = N * nB
    if conv == "qN+a":
        coord = lambda i: (i // N, i % N)  # noqa: E731
        idx = lambda q, a: q * N + a  # noqa: E731
    else:
        coord = lambda i: (i % nB, i // nB)  # noqa: E731
        idx = lambda q, a: a * nB + q  # noqa: E731

    def copies(rows):
        S = {}
        for r in rows:
            cen = [coord(q) for q in r if q < ncen]
            side = [q for q in r if q >= ncen]
            cols = {a for _, a in cen}
            if len(cols) == 1 and len(side) == 1:
                s = frozenset(q for q, _ in cen)
                S.setdefault(s, {}).setdefault(next(iter(cols)), []).append(side[0])
        return S

    SXd, SZd = copies(XE), copies(ZE)
    sX, sZ, sideX, sideZ = [], [], [], []
    for S, out, side in ((SXd, sX, sideX), (SZd, sZ, sideZ)):
        for s, per_a in S.items():
            if sorted(per_a) != list(range(N)):
                return None, f"support of size {len(s)} not copied to all {N} columns"
            labels = [tuple(per_a[a]) for a in range(N)]
            if len(set(labels)) != 1 or len(labels[0]) != 1:
                return None, "copies of one base row with different or repeated side qubits"
            out.append(sorted(s))
            side.append(labels[0][0])
    if len(set(sideX + sideZ)) != len(sideX) + len(sideZ):
        return None, "side qubit used twice"
    if set(sideX + sideZ) != set(range(ncen, nE)):
        return None, "side qubits do not cover the non-central qubits"
    rep = {}
    for name, S, H in (("X", sX, HX), ("Z", sZ, HZ)):
        mS, mH = [mask(r) for r in S], [mask(r) for r in H]
        rS, rH, rSH = rank(mS), rank(mH), rank(mS + mH)
        rep[f"rank_S{name}"], rep[f"rank_H{name}"], rep[f"rank_S{name}+H{name}"] = rS, rH, rSH
        rep[f"m_S{name}"] = len(S)
        if not rS == rH == rSH:
            return None, f"rowsp S_{name} != rowsp H_{name}(B)"
    zcols = [[] for _ in range(nB)]
    for z, r in enumerate(sZ):
        for q in r:
            zcols[q].append(z)
    xcols = [[] for _ in range(nB)]
    for c, r in enumerate(sX):
        for q in r:
            xcols[q].append(c)
    RX = [frozenset([idx(q, a) for q in sX[c]] + [sideX[c]]) for c in range(len(sX)) for a in range(N)]
    RX += [frozenset([idx(q, a) for a in range(N)] + [sideZ[z] for z in zcols[q]]) for q in range(nB)]
    RZ = [frozenset([idx(q, a) for q in sZ[z]] + [sideZ[z]]) for z in range(len(sZ)) for a in range(N)]
    RZ += [frozenset([idx(q, a) for a in range(N)] + [sideX[c] for c in xcols[q]]) for q in range(nB)]
    EX, EZ = [frozenset(r) for r in XE], [frozenset(r) for r in ZE]
    rep.update(n_E=nE, n_rebuilt=ncen + len(sX) + len(sZ), mX_E=len(EX), mX_rebuilt=len(RX),
               mZ_E=len(EZ), mZ_rebuilt=len(RZ),
               X_rows_equal=sorted(map(sorted, EX)) == sorted(map(sorted, RX)),
               Z_rows_equal=sorted(map(sorted, EZ)) == sorted(map(sorted, RZ)))
    rx, rz = rank([mask(r) for r in EX]), rank([mask(r) for r in EZ])
    rep.update(rowsp_X_equal=rx == rank([mask(r) for r in RX]) == rank([mask(r) for r in EX + RX]),
               rowsp_Z_equal=rz == rank([mask(r) for r in RZ]) == rank([mask(r) for r in EZ + RZ]))
    return rep, "ok"


def audit(E, B, N):
    """Return a dict with 'match' (the convention, or None), 'identical_rows', 'rowsp_equal' and the details."""
    out = {"N": N, "n_E": E["n"], "n_B": B["n"], "match": None}
    for conv in ("qN+a", "a*n+q"):
        rep, msg = try_conv(E, B, N, conv)
        out[conv] = rep if rep else msg
        if rep and rep["rowsp_X_equal"] and rep["rowsp_Z_equal"]:
            out["match"] = conv
            out["identical_rows"] = rep["X_rows_equal"] and rep["Z_rows_equal"]
            out["rowsp_equal"] = True
            break
    return out


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    board, e, b, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
    E = json.load(open(os.path.join(board, "codes", e + ".json")))
    B = json.load(open(os.path.join(board, "codes", b + ".json")))
    r = audit(E, B, N)
    r.update(entry=e, base=b)
    print(json.dumps(r))
    sys.exit(0 if r["match"] else 1)


if __name__ == "__main__":
    main()
