#!/usr/bin/env python3
"""Exact d_X, d_Z of the nine codes by three independent programs, compared with data/parameters.json.

    python3 src/run_distances.py [--bin DIR] [--out FILE] [SLUG ...]

Runs, on data/codes/<slug>.hx, .hz:
  bin/graphlike_A NAME.hx NAME.hz --count   (src/graphlike_A: shortest cycles with F_2^k labels; also counts
                                             the minimum-weight logical operators)
  bin/graphlike_B --ref NAME.hx NAME.hz     (src/graphlike_B: double covers; --ref adds the unpruned
                                             reference search)
  python3 src/circuits_check.py NAME.hx NAME.hz SIDE d_SIDE   (all circuits of length <= d of the check
                                             graph, membership test in the other row space)
and requires: n, k, d_X, d_Z as claimed; graphlike_A's witnesses verified; graphlike_B's reference search equal;
the circuit enumeration finds no nontrivial circuit shorter than d and as many of length d as graphlike_A counts
and data/parameters.json states. Exactness of each program: src/graphlike_A/NOTES.md, src/graphlike_B/NOTES.md
(Lemma 12 for the circuit enumeration). One JSON line per code; exit status 1 on any disagreement.
"""
import argparse
import json
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def run(cmd):
    t = time.perf_counter()
    p = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t
    try:
        js = json.loads(p.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError):
        js = None
    return js, p.returncode, round(wall, 3), p.stderr.strip()[-300:]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("slugs", nargs="*")
    ap.add_argument("--bin", default=os.path.join(ROOT, "bin"))
    ap.add_argument("--data", default=os.path.join(ROOT, "data"))
    ap.add_argument("--out", default=None, help="also append the JSON lines to this file")
    a = ap.parse_args(argv)
    with open(os.path.join(a.data, "parameters.json")) as f:
        params = json.load(f)
    fo = open(a.out, "a") if a.out else None
    bad = 0
    for slug in a.slugs or list(params):
        C = params[slug]
        hx = os.path.join(a.data, "codes", slug + ".hx")
        hz = os.path.join(a.data, "codes", slug + ".hz")
        A, rcA, wA, eA = run([os.path.join(a.bin, "graphlike_A"), hx, hz, "--count"])
        B, rcB, wB, eB = run([os.path.join(a.bin, "graphlike_B"), "--ref", hx, hz])
        res = {"slug": slug, "graphlike_A": A, "exit_A": rcA, "wall_A": wA,
               "graphlike_B": B, "exit_B": rcB, "wall_B": wB}
        ok = rcA == 0 and rcB == 0 and A is not None and B is not None
        if ok:
            ok = (A["n"] == B["n"] == C["n"] and A["k"] == B["k"] == C["k"]
                  and A["dX"] == B["dX"] == B["refX"] == C["d_X"] and A["dZ"] == B["dZ"] == B["refZ"] == C["d_Z"]
                  and A["verX"] and A["verZ"] and A["lcheckX"] and A["lcheckZ"]
                  and A["countX"] == C["min_weight_logicals_X"] and A["countZ"] == C["min_weight_logicals_Z"])
        for side in ("X", "Z"):
            d = C["d_" + side]
            cc, rc, w, e = run([sys.executable, os.path.join(HERE, "circuits_check.py"), hx, hz, side, str(d)])
            res["circuits_" + side] = cc
            res["wall_circuits_" + side] = w
            good = (rc == 0 and cc is not None and cc.get("d") == d
                    and all(int(v) == 0 for L, v in cc["nontrivial"].items() if int(L) < d)
                    and int(cc["nontrivial"][str(d)]) == C["min_weight_logicals_" + side])
            res["circuits_ok_" + side] = good
            ok = ok and good
        res["ok"] = ok
        bad += not ok
        line = json.dumps(res)
        print(line, flush=True)
        if fo:
            fo.write(line + "\n")
    print(json.dumps({"codes": len(a.slugs or params), "failed": bad}))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
