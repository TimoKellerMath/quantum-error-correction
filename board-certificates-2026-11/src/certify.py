#!/usr/bin/env python3
"""Rebuild every certified instance from a board checkout and rerun the checks behind its certificate.

    make                                   # bin/graphlike_A, bin/graphlike_B, bin/exactA, bin/exactB
    python3 src/certify.py --board DIR [--slug S ...] [--category C ...] [--long] [--out FILE]

DIR is a checkout of github.com/unitaryfoundation/qldpc-challenge (only DIR/codes/<slug>.json is read). For every
entry of data/entries.json (or only the named ones) the driver
  * reads codes/<slug>.json, checks that its matrices are the certified ones (the board's binding hash
    checks_sha256 recorded in data/entries.json; a distance correction does not change it), recomputes k and
    the CSS condition;
  * writes the instance (common format: first line "n m", then one line of 0-based qubits per check) into a
    scratch directory;
  * runs the programs of each of the entry's methods and checks their verdicts:
      graphlike  graphlike_A and graphlike_B (both sides): the exact side distances d_X, d_Z;
      exact      exactA and exactB on each side: nothing at W = d_side - 1 (with the parity reduction the
                 programs search W - 1 when every logical of the side has even weight) and a word of weight
                 d_side at W = d_side; Tanner-graph group from data/groups/<slug>.gens (both programs verify it);
      partial    exactA and exactB on each side at W = d_lower - 1: nothing found (long runs: only with --long,
                 otherwise the recorded runs of data/results/<slug>.json are reported);
      amplified  src/amp_audit.py: the entry is Amp_N of a presentation of its base (same row spaces); the base is
                 an entry of data/entries.json certified by the methods above, and d_P(entry) = 2 d_P(base) by
                 Proposition 3 of the README;
    an entry with several methods (185-64-5, 185-65-5: graphlike and exact) must get the same values from each;
  * rechecks every word a program returns, in Python, on the board matrices (weight, kernel, outside the other
    row space);
  * compares the outcome with data/results/<slug>.json (verdicts and values, not node counts or times).
With --record the outcome is written to data/results/<slug>.json instead (that is how the package was made).
Output: one JSON line per entry (to stdout and, with --out, to FILE), then a summary; exit status 1 if any entry
fails. Standard library only; the four C programs are built by the Makefile.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
BIN = os.path.join(PKG, "bin")
DATA = os.path.join(PKG, "data")
sys.path.insert(0, HERE)
import amp_audit  # noqa: E402


# --------------------------------------------------------------------------------------------- GF(2), hash
def mask(r):
    v = 0
    for q in r:
        v ^= 1 << q
    return v


def reduce(v, piv):
    while v:
        p = piv.get(v.bit_length() - 1)
        if p is None:
            return v
        v ^= p
    return 0


def echelon(rows):
    piv = {}
    for r in rows:
        v = reduce(mask(r), piv)
        if v:
            piv[v.bit_length() - 1] = v
    return piv


def recheck(doc, side, word):
    """Board convention: an X-type word lies in ker H_Z and outside rowsp H_X; a Z-type word dually."""
    HX, HZ = doc["checks"]["X"], doc["checks"]["Z"]
    comm, stab = (HZ, HX) if side == "X" else (HX, HZ)
    v = mask(word)
    in_ker = all(bin(v & mask(r)).count("1") % 2 == 0 for r in comm)
    return {"weight": len(set(word)), "in_ker": in_ker, "nontrivial": reduce(v, echelon(stab)) != 0}


def checks_sha256(doc):
    """The board's binding hash of (H_X, H_Z) (verify/check_certs.py of the board, reimplemented)."""
    canon = {s: sorted(sorted(int(q) for q in row) for row in doc["checks"][s]) for s in ("X", "Z")}
    blob = json.dumps(canon, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(blob.encode("ascii")).hexdigest()


def basic(doc):
    HX, HZ = doc["checks"]["X"], doc["checks"]["Z"]
    zm = [mask(r) for r in HZ]
    css = all(bin(mask(x) & z).count("1") % 2 == 0 for x in HX for z in zm)
    k = doc["n"] - len(echelon(HX)) - len(echelon(HZ))
    return css, k


# --------------------------------------------------------------------------------------------- programs
def write_inst(prefix, doc):
    for side, ext in (("X", ".hx"), ("Z", ".hz")):
        rows = doc["checks"][side]
        with open(prefix + ext, "w") as f:
            f.write(f"{doc['n']} {len(rows)}\n")
            for r in rows:
                f.write(" ".join(str(q) for q in r) + "\n")


def run(cmd):
    t = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr, round(time.time() - t, 3)


def rel(cmd, work):
    return " ".join(os.path.relpath(c, PKG) if os.path.isabs(c) and not c.startswith(work) else
                    ("WORK/" + os.path.relpath(c, work) if os.path.isabs(c) else c) for c in cmd)


def graphlike(doc, prefix, work):
    out = {}
    for tool in ("graphlike_A", "graphlike_B"):
        cmd = [os.path.join(BIN, tool), prefix + ".hx", prefix + ".hz"]
        rc, so, se, secs = run(cmd)
        line = next((l for l in so.splitlines() if l.startswith("{")), None)
        rec = {"cmd": rel(cmd, work), "rc": rc, "wall": secs, "output": json.loads(line) if line else None}
        if rc or line is None:
            rec["stderr"] = se.strip()[-500:]
        else:
            o = rec["output"]
            rec["dX"], rec["dZ"] = o.get("dX"), o.get("dZ")
            rec["recheck"] = {s: recheck(doc, s, o["wit" + s]) for s in ("X", "Z") if o.get("wit" + s) is not None}
        out[tool] = rec
    return out


def parse_search(tool, so):
    r = {"verdict": None}
    for l in so.splitlines():
        if l.startswith("RESULT"):
            kv = dict(t.split("=", 1) for t in l.split()[1:] if "=" in t)
            r["W"] = int(kv["W"])
            r["searched"] = int(kv.get("Weff", kv.get("searched")))
            r["verdict"] = "found" if kv.get("verdict", kv.get("found")) in ("found", "yes") else \
                ("none" if kv.get("verdict", kv.get("found")) == "none" else kv.get("verdict", kv.get("found")))
            r["nodes"] = int(kv["nodes"])
            r["cpu"] = float(kv["cpu"])
            for key in ("jobs", "mode", "recheck_ker", "recheck_nontrivial", "generated"):
                if key in kv:
                    r[key] = kv[key]
        elif l.startswith("WORD"):
            r["word"] = sorted(int(q) for q in l.split(":", 1)[1].split())
        elif "parity:" in l and "even weight" in l:
            r["parity"] = True
    r.setdefault("parity", False)
    return r


def search(tool, prefix, side, W, gens, work):
    cmd = [os.path.join(BIN, tool), prefix + ".hx", prefix + ".hz", "--side", side, "--W", str(W)]
    if gens:
        cmd += ["--group", gens[tool]]
    if tool == "exactA":
        cmd += ["--mode", "anchor"]
    rc, so, se, secs = run(cmd)
    rec = {"tool": tool, "side": side, "cmd": rel(cmd, work), "rc": rc, "wall": secs}
    rec.update(parse_search(tool, so))
    if rc:
        rec["stderr"] = se.strip()[-500:]
    return rec


def group_files(slug, prefix):
    src = os.path.join(DATA, "groups", slug + ".gens")
    if not os.path.exists(src):
        return None
    with open(src) as f:
        lines = f.read().splitlines()
    if not lines or int(lines[0]) == 0:
        return None
    a, b = prefix + ".gens", prefix + ".gensB"
    shutil.copy(src, a)
    with open(b, "w") as f:
        f.write("# " + lines[0] + " generators\n" + "\n".join(lines[1:]) + "\n")
    return {"exactA": a, "exactB": b}


# --------------------------------------------------------------------------------------------- methods
def m_graphlike(e, doc, prefix, work, bad):
    exp = e["sides"]
    g = graphlike(doc, prefix, work)
    for tool, r in g.items():
        if r.get("output") is None:
            bad(f"{tool} failed (exit {r['rc']})")
            continue
        for s in ("X", "Z"):
            if r[f"d{s}"] != exp[s]:
                bad(f"{tool}: d_{s} = {r[f'd{s}']}, expected {exp[s]}")
            rc = r["recheck"].get(s)
            if not rc or not (rc["in_ker"] and rc["nontrivial"] and rc["weight"] == exp[s]):
                bad(f"{tool}: the {s} word fails the recheck: {rc}")
    return g, {s: exp[s] for s in ("X", "Z")}


def m_search(e, doc, prefix, work, bad, partial):
    exp = e.get("sides") or {}
    gens = group_files(e["slug"], prefix)
    runs = []
    for s in ("X", "Z"):
        Ws = [e["d_lower"] - 1] if partial else [exp[s] - 1, exp[s]]
        for W in Ws:
            for tool in ("exactA", "exactB"):
                r = search(tool, prefix, s, W, gens, work)
                last = not partial and W == exp[s]
                if r["rc"] or r["verdict"] is None:
                    bad(f"{tool} side {s} W={W} failed (exit {r['rc']})")
                elif not last and r["verdict"] != "none":
                    bad(f"{tool} side {s} W={W}: {r['verdict']}, expected none")
                elif last:
                    if r["verdict"] != "found":
                        bad(f"{tool} side {s} W={W}: no word, expected one of weight {W}")
                    else:
                        rc = r["recheck"] = recheck(doc, s, r["word"])
                        if not (rc["in_ker"] and rc["nontrivial"] and rc["weight"] == W):
                            bad(f"{tool} side {s}: word fails the recheck {rc}")
                runs.append(r)
    return runs, ({"d_lower": e["d_lower"]} if partial else {s: exp[s] for s in ("X", "Z")})


def m_amplified(e, doc, board, done, bad):
    b = e["base"]
    B = json.load(open(os.path.join(board, "codes", b + ".json")))
    a = amp_audit.audit(doc, B, e["N"])
    aud = {k2: a[k2] for k2 in ("N", "match", "identical_rows", "rowsp_equal") if k2 in a}
    aud["base"] = b
    aud["detail"] = a.get(a["match"]) if a["match"] else {c: a.get(c) for c in ("qN+a", "a*n+q")}
    if not a["match"]:
        bad(f"amp_audit: {e['slug']} is not Amp_{e['N']} of a presentation of {b}")
    base = done.get(b)
    if base is None or not base["ok"] or base.get("values") is None:
        bad(f"base {b} not certified in this run")
        return aud, None
    bv = base["values"]
    if "d_lower" in bv:
        v = {"d_lower": 2 * bv["d_lower"]}
        if v["d_lower"] != e["d_lower"]:
            bad(f"2 d_lower(base) = {v['d_lower']}, expected {e['d_lower']}")
    else:
        v = {s: 2 * bv[s] for s in ("X", "Z")}
        if any(v[s] != e["sides"][s] for s in ("X", "Z")):
            bad(f"2 d(base) = {v}, expected {e['sides']}")
    return aud, v


def do_entry(e, board, work, long_runs, done):
    slug = e["slug"]
    path = os.path.join(board, "codes", slug + ".json")
    res = {"slug": slug, "category": e["category"], "methods": e["methods"], "ok": True, "problems": []}

    def bad(msg):
        res["ok"] = False
        res["problems"].append(msg)

    if not os.path.exists(path):
        bad(f"codes/{slug}.json not in the board checkout")
        done[slug] = res
        return res
    doc = json.load(open(path))
    res["checks_sha256"] = h = checks_sha256(doc)
    if h != e["checks_sha256"]:
        bad("the matrices differ from the certified ones (checks_sha256)")
        done[slug] = res
        return res
    css, k = basic(doc)
    res["n"], res["k"] = doc["n"], k
    if not css or k != e["k"] or doc["n"] != e["n"]:
        bad(f"CSS {css}, n {doc['n']}, k {k}: expected n {e['n']}, k {e['k']}")
    res["claim"] = {s: doc["distance"][s]["value"] for s in ("X", "Z")}
    prefix = os.path.join(work, slug)
    write_inst(prefix, doc)
    res["runs"], vals = {}, []
    for m in e["methods"]:
        if m == "graphlike":
            res["runs"][m], v = m_graphlike(e, doc, prefix, work, bad)
        elif m == "exact":
            res["runs"][m], v = m_search(e, doc, prefix, work, bad, False)
        elif m == "partial":
            if not long_runs:
                res["runs"][m] = "recorded in data/results (rerun with --long)"
                res["skipped_long"] = True
                v = {"d_lower": e["d_lower"]}
            else:
                res["runs"][m], v = m_search(e, doc, prefix, work, bad, True)
        elif m == "amplified":
            res["audit"], v = m_amplified(e, doc, board, done, bad)
        else:
            bad(f"unknown method {m}")
            v = None
        vals.append(v)
    if any(v != vals[0] for v in vals):
        bad(f"the methods disagree: {vals}")
    res["values"] = vals[0] if res["ok"] else None
    cor = e.get("correction")
    if cor:
        rc = recheck(doc, cor["side"], cor["witness"])
        res["correction_witness"] = rc
        if not (rc["in_ker"] and rc["nontrivial"] and rc["weight"] == cor["new"] == e["sides"][cor["side"]]):
            bad(f"correction witness of side {cor['side']} fails: {rc}")
    done[slug] = res
    return res


def strip(res):
    """The part of a result that must agree between runs (no node counts, times, words)."""
    out = {k: res.get(k) for k in ("slug", "category", "methods", "ok", "checks_sha256", "n", "k", "values")}
    if res.get("audit"):
        out["audit"] = {k: res["audit"].get(k) for k in ("base", "match", "identical_rows", "rowsp_equal")}
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--board", required=True, help="checkout of the board (directory with codes/)")
    ap.add_argument("--slug", nargs="*", help="only these entries (bases of amplified entries are added)")
    ap.add_argument("--category", nargs="*", help="only these categories of data/entries.json")
    ap.add_argument("--long", action="store_true", help="also rerun the long partial searches (hours)")
    ap.add_argument("--record", action="store_true", help="write data/results/<slug>.json")
    ap.add_argument("--out", help="also write the JSON lines to this file")
    ap.add_argument("--keep", help="keep the instances in this directory")
    a = ap.parse_args()
    for b in ("graphlike_A", "graphlike_B", "exactA", "exactB"):
        if not os.access(os.path.join(BIN, b), os.X_OK):
            sys.exit(f"bin/{b} missing: run make first")
    man = json.load(open(os.path.join(DATA, "entries.json")))
    ents = man["entries"]
    want = list(ents)
    if a.category:
        want = [s for s in want if ents[s]["category"] in a.category]
    if a.slug:
        want = [s for s in want if s in a.slug]
        missing = [s for s in a.slug if s not in ents]
        if missing:
            sys.exit(f"not in data/entries.json: {missing}")
    for s in list(want):                       # bases first
        b = ents[s].get("base")
        if b and b not in want:
            want.append(b)
    want.sort(key=lambda s: ("amplified" in ents[s]["methods"], s))
    work = a.keep or tempfile.mkdtemp(prefix="certify-")
    os.makedirs(work, exist_ok=True)
    outf = open(a.out, "w") if a.out else None
    done, nbad, nlong = {}, 0, 0
    t0 = time.time()
    for s in want:
        res = do_entry(ents[s], os.path.abspath(a.board), work, a.long, done)
        recf = os.path.join(DATA, "results", s + ".json")
        if a.record and not res.get("skipped_long"):
            os.makedirs(os.path.dirname(recf), exist_ok=True)
            with open(recf, "w") as f:
                json.dump(res, f, indent=1)
                f.write("\n")
        elif os.path.exists(recf):
            with open(recf) as f:
                old = json.load(f)
            res["agrees_with_record"] = strip(old) == strip(res)
            if not res["agrees_with_record"]:
                res["ok"] = False
                res["problems"].append("differs from data/results/" + s + ".json")
        nlong += bool(res.get("skipped_long"))
        nbad += not res["ok"]
        line = json.dumps(res, separators=(",", ":"))
        print(line)
        if outf:
            outf.write(line + "\n")
    if outf:
        outf.close()
    if not a.keep:
        shutil.rmtree(work, ignore_errors=True)
    print(f"# {len(want)} entries, {len(want) - nbad} ok, {nbad} failed, {nlong} long partial runs not rerun "
          f"(--long); {time.time() - t0:.1f} s", file=sys.stderr)
    sys.exit(1 if nbad else 0)


if __name__ == "__main__":
    main()
