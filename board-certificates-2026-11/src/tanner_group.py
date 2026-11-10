#!/usr/bin/env python3
"""Generators of the automorphism group of the Tanner graph of a board entry (needs nauty's dreadnaut).

    python3 src/tanner_group.py BOARD_DIR SLUG OUT.gens

Vertex classes qubits | X-checks | Z-checks are fixed; the generators are restricted to the qubits and checked
to map the set of X-rows and the set of Z-rows onto themselves. Output in the format of exactA (first line the
number of generators, then one permutation per line, image of qubit q in position q); exactB reads the same file
with the first line commented out (src/certify.py writes that copy). Both search programs verify every generator
again before using it, so a wrong file stops the run rather than weakening it. Standard library and dreadnaut.
"""
import json
import os
import subprocess
import sys


def generators(n, hx, hz):
    N = n + len(hx) + len(hz)
    lines = ["l=0", f"n={N} g"]
    for i, r in enumerate(hx):
        lines.append(f"{n + i}: " + " ".join(map(str, r)) + ";")
    for i, r in enumerate(hz):
        lines.append(f"{n + len(hx) + i}: " + " ".join(map(str, r)) + ";")
    lines.append(".")
    cells = [f"0:{n - 1}"]
    if hx:
        cells.append(f"{n}:{n + len(hx) - 1}")
    if hz:
        cells.append(f"{n + len(hx)}:{N - 1}")
    lines.append("f=[" + "|".join(cells) + "]")
    lines += ["x", "q"]
    out = subprocess.run(["dreadnaut"], input="\n".join(lines) + "\n", capture_output=True, text=True,
                         check=True).stdout
    gens, cur, grpsize = [], None, None
    for line in out.splitlines():
        s = line.strip()
        if line.startswith("("):
            if cur is not None:
                gens.append(cur)
            cur = s
        elif cur is not None and line[:1].isspace() and s[:1] in "(0123456789":
            cur += " " + s
        else:
            if cur is not None:
                gens.append(cur)
                cur = None
            if "grpsize=" in s:
                grpsize = s.split("grpsize=")[1].split(";")[0]
    if cur is not None:
        gens.append(cur)
    perms = []
    for g in gens:
        p = list(range(N))
        for cyc in g.replace(")", ")\n").split("\n"):
            cyc = cyc.strip().strip("()")
            if not cyc:
                continue
            pts = list(map(int, cyc.split()))
            for a, b in zip(pts, pts[1:] + pts[:1]):
                p[a] = b
        perms.append(p)
    setx = {tuple(sorted(r)) for r in hx}
    setz = {tuple(sorted(r)) for r in hz}
    qperms = []
    for p in perms:
        assert all(p[q] < n for q in range(n))
        assert {tuple(sorted(p[q] for q in r)) for r in hx} == setx
        assert {tuple(sorted(p[q] for q in r)) for r in hz} == setz
        if p[:n] != list(range(n)):
            qperms.append(p[:n])
    return qperms, grpsize


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    board, slug, out = sys.argv[1:]
    doc = json.load(open(os.path.join(board, "codes", slug + ".json")))
    qperms, grpsize = generators(doc["n"], doc["checks"]["X"], doc["checks"]["Z"])
    with open(out, "w") as fh:
        fh.write(f"{len(qperms)}\n")
        for qp in qperms:
            fh.write(" ".join(map(str, qp)) + "\n")
    print(json.dumps({"slug": slug, "generators": len(qperms), "grpsize": grpsize}))


if __name__ == "__main__":
    main()
