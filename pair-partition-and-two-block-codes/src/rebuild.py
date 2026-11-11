#!/usr/bin/env python3
"""Rebuild the six codes of this directory from their definitions and check them, their witnesses and the records.

    python3 src/rebuild.py [SLUG ...] [--entry FILE ...]

Run from anywhere (Python 3.8 or later, standard library only); the package is the parent directory of src/.
Every output line starts with PASS, FAIL or INFO; the last line is "RESULT: all checks passed" or
"RESULT: <m> checks failed", and the exit status is 0 only in the first case.

For each code <slug>/ (all six if no SLUG is given):

1. Definition. Pair-partition codes (Part A of README.md): the template (rows of weight 6, column degrees 2 or 3 on
   each side, six columns of degree 2 per side forming a connected unicyclic graph on the rows, T_X T_Z^T = 0), the
   pairing of every overlap X_i & Z_j and the paired-difference equations, the girth rule of the pairing array, and
   the exponents against the job's exponent string. Two-block code (Part B): the group law of C12 x|_5 C28 (closure,
   identity, inverses, associativity by Light's test on the generators x, y), the relations of the presentation,
   the element numbering, and a, b against their indices.
2. Matrices. H_X and H_Z built from the definition must equal hx.txt and hz.txt row by row (same rows, same order).
   Checked on them: n, the ranks and k, CSS orthogonality, check weights, qubit degrees, connectivity of the Tanner
   graph, whether two rows of H_X (of H_Z) share two qubits, and whether the all-one vector lies in the row spaces.
3. Witnesses. witness_X.txt lies in ker H_Z outside rowsp H_X, witness_Z.txt in ker H_X outside rowsp H_Z, and their
   weights are d_X and d_Z of definition.json.
4. Records. The instance files of the exact searches equal hx.txt, hz.txt and the witnesses; every permutation in
   the generator files maps the rows of H_X onto rows of H_X and the rows of H_Z onto rows of H_Z (the closure of
   the generators is computed: order and orbits); the two programs' logs of both sides report no word at
   W = d - 1 with the parity decision, n, k and ranks found here; the JSON records of the job agree with the logs;
   records/searches.json agrees with the logs. Together with item 3 this gives d_X = d_Z = d.
5. Third-party tools. The MatrixMarket inputs equal H_X, H_Z, the logical bases L_X, L_Z complete the row spaces to
   the kernels, every run listed in third_tools/runs.json reports what runs.json states, and every word the tools
   returned is a nontrivial logical of the stated side and weight d.

--entry FILE: FILE is a board entry codes/<slug>.json. Its checks must equal hx.txt, hz.txt of <slug>/ row by row,
its n, k, d and side values must be the ones established here, and its witnesses must be nontrivial logicals of
the stated weights.
"""
import argparse
import hashlib
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SLUGS = ["294-42-11", "322-46-12", "434-62-13", "518-74-14", "742-106-15", "672-28-24"]
FAILURES = []


def report(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        FAILURES.append(what)
    return ok


def info(what):
    print("INFO " + what)


def path(*parts):
    return os.path.join(ROOT, *parts)


def read_text(*parts):
    with open(path(*parts), encoding="utf-8") as f:
        return f.read()


def load(*parts):
    return json.loads(read_text(*parts))


def sha256(*parts):
    with open(path(*parts), "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


# ------------------------------------------------------------------------------------------------ GF(2)

def popcount(v):
    return bin(v).count("1")


def vec(support):
    v = 0
    for q in support:
        v |= 1 << q
    return v


def support(v):
    out, i = [], 0
    while v:
        if v & 1:
            out.append(i)
        v >>= 1
        i += 1
    return out


class Span:
    """Row space over GF(2); rows are Python integers with distinct leading bits."""

    def __init__(self, rows=()):
        self.piv = {}
        for r in rows:
            self.add(r)

    def reduce(self, v):
        while v:
            r = self.piv.get(v.bit_length() - 1)
            if r is None:
                return v
            v ^= r
        return 0

    def add(self, v):
        v = self.reduce(v)
        if v:
            self.piv[v.bit_length() - 1] = v
            return True
        return False

    def contains(self, v):
        return self.reduce(v) == 0

    def rank(self):
        return len(self.piv)


def in_kernel(rows, v):
    return all(popcount(r & v) % 2 == 0 for r in rows)


def read_rows(*parts):
    rows = []
    for ln in read_text(*parts).splitlines():
        if ln.strip():
            rows.append([int(t) for t in ln.split()])
    return rows


def read_word(*parts):
    return [int(t) for t in read_text(*parts).split()]


def rows_text(rows):
    return "".join(" ".join(str(q) for q in r) + "\n" for r in rows)


class Code:
    def __init__(self, n, hx, hz):
        self.n, self.hx, self.hz = n, hx, hz
        self.HX = [vec(r) for r in hx]
        self.HZ = [vec(r) for r in hz]
        self.SX, self.SZ = Span(self.HX), Span(self.HZ)
        self.k = n - self.SX.rank() - self.SZ.rank()

    def logical(self, side, word):
        """(in the right kernel, outside the same-side row space) for an X- or a Z-logical."""
        v = vec(word)
        if side == "X":
            return in_kernel(self.HZ, v), not self.SX.contains(v)
        return in_kernel(self.HX, v), not self.SZ.contains(v)

    def check_word(self, side, word, weight, label):
        ok_set = len(set(word)) == len(word) and all(0 <= q < self.n for q in word)
        ker, nontriv = self.logical(side, word)
        other = "H_Z" if side == "X" else "H_X"
        return report(ok_set and ker and nontriv and len(word) == weight,
                      f"{label}: {len(word)} distinct qubits, in ker {other}: {ker}, outside rowsp H_{side}: "
                      f"{nontriv}, weight {len(word)} (expected {weight})")


# ------------------------------------------------------------------------------------------------ Part A

def template_checks(slug, D):
    J, L, P = D["J"], D["L"], D["P"]
    X, Z = D["template"]["X"], D["template"]["Z"]
    ok = len(X) == len(Z) == J and all(len(r) == 6 and r == sorted(set(r)) and all(0 <= u < L for u in r)
                                        for r in X + Z)
    report(ok, f"{slug}: template of {J} + {J} rows of weight 6 on {L} block columns")
    for side, T in (("X", X), ("Z", Z)):
        deg = [sum(u in r for r in T) for u in range(L)]
        two = [u for u in range(L) if deg[u] == 2]
        report(set(deg) <= {2, 3} and len(two) == 6,
               f"{slug}: T_{side} column degrees {deg}, six of degree 2: {two}")
        # the degree-2 columns as edges between the two rows containing them: connected and unicyclic
        edges = [tuple(i for i in range(J) if u in T[i]) for u in two]
        parent = list(range(J))

        def find(a):
            while parent[a] != a:
                parent[a] = parent[parent[a]]
                a = parent[a]
            return a
        cycles = 0
        for a, b in edges:
            ra, rb = find(a), find(b)
            if ra == rb:
                cycles += 1
            else:
                parent[ra] = rb
        comps = len({find(i) for i in range(J)})
        report(comps == 1 and cycles == 1 and len(edges) == J,
               f"{slug}: degree-2 graph of T_{side} ({len(edges)} edges on {J} rows) connected and unicyclic")
    report(all(len(set(x) & set(z)) % 2 == 0 for x in X for z in Z), f"{slug}: T_X T_Z^T = 0 (even overlaps)")

    # exponents: rows aligned with the template rows, and the job's exponent string
    eX, eZ = D["exponents"]["X"], D["exponents"]["Z"]
    ok = (len(eX) == len(eZ) == J and all(len(e) == 6 for e in eX + eZ)
          and all(0 <= t < P for e in eX + eZ for t in e))
    flat = [t for e in eX for t in e] + [t for e in eZ for t in e]
    report(ok and flat == [int(t) for t in D["job"]["exps"].split(",")],
           f"{slug}: exponents ({len(flat)} edges, values mod {P}) equal the job's exponent string")
    ex = {(i, u): e for i in range(J) for u, e in zip(X[i], eX[i])}
    ez = {(j, u): e for j in range(J) for u, e in zip(Z[j], eZ[j])}

    # pair partition (paired differences) and the girth rule of the pairing array
    cells = {}
    for c in D["pairing"]:
        cells[(c["X"], c["Z"])] = [tuple(p) for p in c["pairs"]]
    want = {(i, j) for i in range(J) for j in range(J) if set(X[i]) & set(Z[j])}
    report(set(cells) == want, f"{slug}: pairings given for exactly the {len(want)} nonempty overlaps X_i & Z_j")
    bad = []
    for (i, j), pairs in cells.items():
        flatp = [u for p in pairs for u in p]
        if sorted(flatp) != sorted(set(X[i]) & set(Z[j])) or len(set(flatp)) != len(flatp):
            bad.append((i, j, "not a partition into pairs"))
            continue
        for u, v in pairs:
            if (ex[(i, u)] - ez[(j, u)] - ex[(i, v)] + ez[(j, v)]) % P:
                bad.append((i, j, (u, v)))
    report(not bad, f"{slug}: every overlap is partitioned into pairs (u, v) with "
                    f"e_iu - d_ju = e_iv - d_jv mod {P} {bad[:3] if bad else ''}")
    dup = []
    for fix in (0, 1):
        for a in range(J):
            seen = {}
            for key, pairs in cells.items():
                if key[fix] != a:
                    continue
                for p in pairs:
                    p = tuple(sorted(p))
                    if p in seen and seen[p] != key:
                        dup.append((key, seen[p], p))
                    seen[p] = key
    report(not dup, f"{slug}: girth rule of the pairing array (no pair twice in an array row or column) "
                    f"{dup[:3] if dup else ''}")
    return ex, ez


def build_l3(D):
    P = D["P"]
    X, Z = D["template"]["X"], D["template"]["Z"]
    eX, eZ = D["exponents"]["X"], D["exponents"]["Z"]

    def rows(T, E):
        return [sorted(u * P + (r + e) % P for u, e in zip(T[i], E[i])) for i in range(len(T)) for r in range(P)]
    return D["L"] * P, rows(X, eX), rows(Z, eZ)


# ------------------------------------------------------------------------------------------------ Part B

def build_2bga(slug, D):
    g = D["group"]
    m, q, r = g["m"], g["q"], g["r"]
    order = m * q
    els = [(i, j) for i in range(m) for j in range(q)]

    def mul(a, b):
        return ((a[0] + pow(r, a[1], m) * b[0]) % m, (a[1] + b[1]) % q)

    def index(e):
        i, j = e
        return (j % 2) * (order // 2) + i * (q // 2) + j // 2
    e0, x, y = (0, 0), (1, 0), (0, 1)
    idx = [index(e) for e in els]
    report(sorted(idx) == list(range(order)), f"{slug}: element index {g['index']} is a bijection onto 0..{order - 1}")
    report(all(mul(e0, e) == e == mul(e, e0) for e in els)
           and all(sum(mul(e, f) == e0 for f in els) == 1 for e in els),
           f"{slug}: identity x^0 y^0 and two-sided inverses (every element has exactly one right inverse)")
    report(all(mul(mul(a, s), b) == mul(a, mul(s, b)) for s in (x, y) for a in els for b in els),
           f"{slug}: associativity by Light's test on the generators x, y ({2 * order * order} products)")

    def power(e, k):
        out = e0
        for _ in range(k):
            out = mul(out, e)
        return out

    def ordr(e):
        k, t = 1, e
        while t != e0:
            t, k = mul(t, e), k + 1
        return k
    yinv = power(y, q - 1)
    report(ordr(x) == m and ordr(y) == q and mul(mul(y, x), yinv) == power(x, r) and pow(r, 2, m) == 1,
           f"{slug}: |x| = {m}, |y| = {q}, y x y^-1 = x^{r} (and {r}^2 = 1 mod {m}, so the rule is a semidirect "
           f"product); every element is x^i y^j, so x and y generate")
    centre = [e for e in els if all(mul(e, f) == mul(f, e) for f in (x, y))]
    info(f"{slug}: |G| = {order}, centre of order {len(centre)}")

    def word(e):
        return power(x, e[0]) if e[1] == 0 else mul(power(x, e[0]), power(y, e[1]))
    inv = {e: next(f for f in els if mul(e, f) == e0) for e in els}
    a = [tuple(t) for t in D["a"]["exponents"]]
    b = [tuple(t) for t in D["b"]["exponents"]]
    report([index(word(e)) for e in a] == D["a"]["indices"] and [index(word(e)) for e in b] == D["b"]["indices"],
           f"{slug}: a = {D['a']['word']} has indices {D['a']['indices']}, b = {D['b']['word']} has indices "
           f"{D['b']['indices']}")
    a, b = [word(e) for e in a], [word(e) for e in b]
    hx, hz = [None] * order, [None] * order
    for c in els:
        hx[index(c)] = sorted([index(mul(inv[s], c)) for s in a] + [order + index(mul(c, inv[t])) for t in b])
        hz[index(c)] = sorted([index(mul(c, t)) for t in b] + [order + index(mul(s, c)) for s in a])
    return 2 * order, hx, hz


# ------------------------------------------------------------------------------------------------ the code

def code_checks(slug, D, code):
    n, k = code.n, code.k
    P = D["parameters"]
    report(n == P["n"] and k == P["k"],
           f"{slug}: n = {n}, rank H_X = {code.SX.rank()}, rank H_Z = {code.SZ.rank()}, k = {k} "
           f"(definition: n = {P['n']}, k = {P['k']})")
    report(all(popcount(a & b) % 2 == 0 for a in code.HX for b in code.HZ), f"{slug}: H_X H_Z^T = 0")
    wts = sorted({len(r) for r in code.hx + code.hz})
    report(wts == [P["w"]] and all(len(set(r)) == len(r) for r in code.hx + code.hz),
           f"{slug}: every check has weight {P['w']} ({len(code.hx)} X-checks, {len(code.hz)} Z-checks)")
    for side, rows in (("X", code.hx), ("Z", code.hz)):
        deg = [0] * n
        for r in rows:
            for qb in r:
                deg[qb] += 1
        hist = {t: deg.count(t) for t in sorted(set(deg))}
        report(sorted(hist) == P["qubit_degrees"][side],
               f"{slug}: qubit {side}-degrees {hist} (expected values {P['qubit_degrees'][side]})")
    # connectivity of the Tanner graph (qubits, X-checks, Z-checks)
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for r in code.hx + code.hz:
        for qb in r[1:]:
            ra, rb = find(r[0]), find(qb)
            if ra != rb:
                parent[ra] = rb
    report(len({find(i) for i in range(n)}) == 1, f"{slug}: Tanner graph connected")
    for side, rows in (("X", code.HX), ("Z", code.HZ)):
        share2 = sum(1 for i in range(len(rows)) for j in range(i) if popcount(rows[i] & rows[j]) >= 2)
        if D["kind"] == "pair-partition":
            report(share2 == 0, f"{slug}: no two rows of H_{side} share two qubits (no 4-cycle within H_{side})")
        else:
            info(f"{slug}: pairs of rows of H_{side} sharing at least two qubits: {share2}")
    one = (1 << n) - 1
    par = {"X": code.SX.contains(one), "Z": code.SZ.contains(one)}
    report(par == {s: P["all_one_in_rowspace"][s] for s in "XZ"},
           f"{slug}: all-one vector in rowsp H_X: {par['X']}, in rowsp H_Z: {par['Z']} (so every logical of side Z, "
           f"resp. X, has even weight iff True)")
    return par


def witness_checks(slug, D, code):
    P = D["parameters"]
    for side in "XZ":
        code.check_word(side, read_word(slug, f"witness_{side}.txt"), P[f"d_{side}"],
                        f"{slug}: witness_{side}.txt is a nontrivial {side}-logical")


# ------------------------------------------------------------------------------------------------ records

def read_perms(*parts):
    perms = []
    for ln in read_text(*parts).splitlines():
        t = ln.split()
        if not t or t[0].startswith("#"):
            continue
        perms.append([int(v) for v in t])
    return perms


def closure(gens, n):
    ident = tuple(range(n))
    seen, todo = {ident}, [ident]
    while todo:
        g = todo.pop()
        for s in gens:
            h = tuple(s[g[i]] for i in range(n))
            if h not in seen:
                seen.add(h)
                todo.append(h)
                if len(seen) > 100000:
                    return None
    return seen


def orbits(gens, n):
    parent = list(range(n))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    for s in gens:
        for i in range(n):
            ra, rb = find(i), find(s[i])
            if ra != rb:
                parent[ra] = rb
    out = {}
    for i in range(n):
        out.setdefault(find(i), []).append(i)
    return sorted(out.values())


RES_A = re.compile(r"^RESULT side=(\w) W=(\d+) Weff=(\d+) verdict=(\w+) nodes=(\d+) .*cpu=([\d.]+)")
RES_B = re.compile(r"^RESULT side=(\w) W=(\d+) searched=(\d+) found=(\w+) nodes=(\d+) .*cpu=([\d.]+)")
CON_A = re.compile(r"^CONCLUSION none of weight <= (\d+); d_(\w) >= (\d+)")
CON_B = re.compile(r"^CONCLUSION: no nontrivial (\w)-logical of weight <= (\d+), d_(\w) >= (\d+)")


def parse_log(text, impl):
    """The fields of an exactA or exactB log that the checks below use."""
    out = {"impl": impl}
    for ln in text.splitlines():
        s = ln.strip()
        if impl == "A":
            m = re.match(r"exactA: side (\w), n=(\d+)", s)
            if m:
                out["side_hdr"], out["n"] = m.group(1), int(m.group(2))
            m = re.match(r"rank H_(\w)=(\d+), rank H_(\w)=(\d+), k=(\d+)", s)
            if m:
                out["rank"] = {m.group(1): int(m.group(2)), m.group(3): int(m.group(4))}
                out["k"] = int(m.group(5))
            m = re.match(r"group: (\d+) generators from (\S+), preserving", s)
            if m:
                out["group_file"] = m.group(2).rstrip(",")
            m = re.match(r"mode (\w+)", s)
            if m:
                out["mode"] = m.group(1)
            m = RES_A.match(s)
            if m:
                out.update(side=m.group(1), W=int(m.group(2)), searched=int(m.group(3)), verdict=m.group(4),
                           nodes=int(m.group(5)), cpu=float(m.group(6)))
            m = CON_A.match(s)
            if m:
                out["conclusion"] = (m.group(2), int(m.group(3)))
        else:
            m = re.match(r"exactB: (\S+) (\S+) side (\w): n=(\d+) k=(\d+) \| A \(H_(\w)\): m=\d+ rank=(\d+) .*"
                         r"\| B \(H_(\w)\): m=\d+ rank=(\d+)", s)
            if m:
                out["files"] = (m.group(1), m.group(2))
                out["side_hdr"], out["n"], out["k"] = m.group(3), int(m.group(4)), int(m.group(5))
                out["rank"] = {m.group(6): int(m.group(7)), m.group(8): int(m.group(9))}
            m = re.match(r"group: (\d+) generators from (\S+) verified .*; (\d+) orbits", s)
            if m:
                out["group_file"], out["orbits"] = m.group(2), int(m.group(3))
            m = re.match(r"run \d+: anchor \d+, excluded \d+, nodes (\d+)", s)
            if m:
                out.setdefault("run_nodes", []).append(int(m.group(1)))
            m = RES_B.match(s)
            if m:
                out.update(side=m.group(1), W=int(m.group(2)), searched=int(m.group(3)), verdict=m.group(4),
                           nodes=int(m.group(5)), cpu=float(m.group(6)))
            m = CON_B.match(s)
            if m:
                out["conclusion"] = (m.group(3), int(m.group(4)))
        if s.startswith("parity:"):
            out["parity_even"] = (" in rowsp" in s) and ("not in rowsp" not in s)
        if s.startswith("exit="):
            out["exit"] = int(s.split("=")[1])
    return out


def record_checks(slug, D, code, par):
    R = D["records"]
    P = D["parameters"]
    inst = R["instance"]
    hdr_x = f"{code.n} {len(code.hx)}\n"
    report(read_text(slug, "records", inst + ".hx") == hdr_x + rows_text(code.hx)
           and read_text(slug, "records", inst + ".hz") == f"{code.n} {len(code.hz)}\n" + rows_text(code.hz),
           f"{slug}: records/{inst}.hx, .hz are hx.txt, hz.txt with the header 'n m' (the searched instance)")
    for side in "XZ":
        w = read_word(slug, "records", f"{inst}.w{side}")
        report(sorted(w) == sorted(read_word(slug, f"witness_{side}.txt")),
               f"{slug}: records/{inst}.w{side} is witness_{side}.txt")
    # generators
    gens_all = []
    for fn in R["generator_files"]:
        perms = read_perms(slug, "records", fn)
        if fn.endswith(".gensA"):
            cnt = int(perms[0][0])
            perms = perms[1:]
            report(cnt == len(perms), f"{slug}: {fn} announces {cnt} permutations and has {len(perms)}")
        rowsx, rowsz = set(map(tuple, code.hx)), set(map(tuple, code.hz))
        ok = all(sorted(p) == list(range(code.n)) for p in perms)
        ok = ok and all(tuple(sorted(p[qb] for qb in r)) in rowsx for p in perms for r in code.hx)
        ok = ok and all(tuple(sorted(p[qb] for qb in r)) in rowsz for p in perms for r in code.hz)
        report(ok, f"{slug}: the {len(perms)} permutations of records/{fn} map the rows of H_X onto rows of H_X "
                   f"and the rows of H_Z onto rows of H_Z")
        gens_all.append({tuple(p) for p in perms})
    report(all(g == gens_all[0] for g in gens_all),
           f"{slug}: the generator files {', '.join(R['generator_files'])} hold the same permutations")
    gens = [list(p) for p in sorted(gens_all[0])]
    grp = closure(gens, code.n)
    orb = orbits(gens, code.n)
    report(grp is not None and len(grp) == R["group_order"] and sorted(len(o) for o in orb) == R["orbit_sizes"],
           f"{slug}: the generated group has order {len(grp) if grp else '> 100000'} "
           f"(expected {R['group_order']}) and {len(orb)} qubit orbits of sizes {sorted(set(len(o) for o in orb))}")
    b = R["orbit_block"]
    report(orb == [list(range(i * b, (i + 1) * b)) for i in range(code.n // b)],
           f"{slug}: the orbits are the {code.n // b} ranges of {b} consecutive qubits")

    # the logs of the four searches
    summary = []
    for run in R["searches"]:
        log = read_text(slug, "records", run["log"])
        L = parse_log(log, run["impl"])
        side, W = run["side"], P["d"] - 1
        rel = "Z" if side == "X" else "X"
        even = par[rel]
        want_searched = W - 1 if (even and W % 2 == 1) else W
        ok = (L.get("side") == side == L.get("side_hdr") and L.get("W") == W and L.get("verdict") == "none"
              and L.get("searched") == want_searched and L.get("parity_even") == even
              and L.get("conclusion") == (side, P["d"]) and L.get("n") == code.n and L.get("k") == code.k
              and L.get("rank") == {"X": code.SX.rank(), "Z": code.SZ.rank()}
              and L.get("group_file") == run["args"][run["args"].index("--group") + 1]
              and L.get("exit", 0) == 0)
        if run["impl"] == "B":
            ok = ok and L.get("files") == (inst + ".hx", inst + ".hz") and sum(L.get("run_nodes", [])) == L["nodes"]
        else:
            want_mode = run["args"][run["args"].index("--mode") + 1] if "--mode" in run["args"] else "anchor"
            ok = ok and L.get("mode") == want_mode
        report(ok, f"{slug}: records/{run['log']}: exact{run['impl']}, side {side}, W = {W}, searched "
                   f"{L.get('searched')} (parity: all logicals even: {even}), verdict none, d_{side} >= {P['d']}, "
                   f"{L.get('nodes')} nodes; n, k, ranks and generator file as computed here")
        summary.append({"impl": run["impl"], "side": side, "W": L.get("W"), "searched": L.get("searched"),
                        "verdict": L.get("verdict"), "nodes": L.get("nodes"), "cpu_s": L.get("cpu"),
                        "lower_bound": L.get("conclusion", (None, None))[1], "log": run["log"],
                        "args": run["args"]})
        if run.get("json"):
            J = load(slug, "records", run["json"])
            rd = J["rounds"]
            ok = (J["verdict"] == "none" and J["lower"] == J["upper"] == P["d"] and J["exact"] is True
                  and J["dcode"] == P["d"] and J["side"] == side and J["impl"] == run["impl"] and len(rd) == 1
                  and rd[0]["W"] == W and rd[0]["verdict"] == "none" and rd[0]["nodes"] == L.get("nodes")
                  and rd[0]["log"] == run["log"] and rd[0]["returncode"] == 0)
            report(ok, f"{slug}: records/{run['json']} (the job's record) agrees with the log: one round at W = {W}, "
                       f"none, lower = upper = {P['d']}")
    for side in "XZ":
        lows = sorted(s["impl"] for s in summary if s["side"] == side and s["lower_bound"] == P["d"])
        report(lows == ["A", "B"] and P[f"d_{side}"] == P["d"],
               f"{slug}: d_{side} >= {P['d']} by exactA and by exactB, d_{side} <= {P['d']} by witness_{side}.txt, "
               f"so d_{side} = {P['d']}")
    S = load(slug, "records", "searches.json")
    report(S["runs"] == summary, f"{slug}: records/searches.json agrees with the four logs")


# ------------------------------------------------------------------------------------------------ third tools

MM = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d+)\s*$")
CCW = re.compile(r"CC codeword of weight (\d+) \(1-based columns\):((?: \d+)+)")


def read_mtx(*parts):
    lines = [ln for ln in read_text(*parts).splitlines() if ln and not ln.startswith("%")]
    m, n, nnz = map(int, lines[0].split())
    rows = [[] for _ in range(m)]
    for ln in lines[1:1 + nnz]:
        i, j, v = map(int, ln.split())
        if v % 2:
            rows[i - 1].append(j - 1)
    return m, n, [sorted(r) for r in rows]


def read_nz(*parts):
    words = []
    for ln in read_text(*parts).splitlines():
        if not ln or ln.startswith("%"):
            continue
        t = [int(v) for v in ln.split()]
        words.append((t[0], [v - 1 for v in t[1:]]))
    return words


def third_checks(slug, D, code):
    T = load(slug, "third_tools", "runs.json")
    d = D["parameters"]["d"]
    for name, fn in T["inputs"].items():
        m, n, rows = read_mtx(slug, "third_tools", fn)
        if name in ("HX", "HZ"):
            report(n == code.n and rows == (code.hx if name == "HX" else code.hz),
                   f"{slug}: third_tools/{fn} is H_{name[1]} ({m} x {n}, rows in the order of hx.txt/hz.txt)")
        else:
            side = name[1]
            vecs = [vec(r) for r in rows]
            S = Span(code.HX if side == "X" else code.HZ)
            other = code.HZ if side == "X" else code.HX
            ok = n == code.n and m == code.k and all(in_kernel(other, v) for v in vecs)
            ok = ok and all(S.add(v) for v in vecs)
            report(ok, f"{slug}: third_tools/{fn} is a basis of {code.k} {side}-logical classes (rows in ker "
                       f"H_{'Z' if side == 'X' else 'X'}, independent modulo rowsp H_{side})")
    for run in T["runs"]:
        side, tool = run["side"], run["tool"]
        lab = f"{slug}: third_tools/{run['files'][0]} ({tool}, side {side})"
        text = "\n".join(read_text(slug, "third_tools", fn) for fn in run["files"])
        if tool == "dist-m4ri":
            res = [tuple(map(int, m.groups())) for m in map(MM.match, text.splitlines()) if m]
            ok = len(res) == 1 and res[0][0] == run["lower_bound"] and res[0][1] == run["found_weight"]
            hdr = re.search(r"finH='([^']+)'", text)
            want_h = "HX.mtx" if side == "Z" else "HZ.mtx"
            ok = ok and hdr is not None and hdr.group(1).endswith(want_h)
            sym = bool(re.search(r"start=\d|start list", text))
            sw = re.search(r"CC from w=(\d+)", text)
            ok = ok and sym == run["symmetry"] and (int(sw.group(1)) if sw else 1) == run["start_weight"]
            report(ok, f"{lab}: result line {res[0] if res else None} = (lower bound {run['lower_bound']}, found "
                       f"weight {run['found_weight']}), finH = {want_h}, "
                       f"{'anchored at orbit representatives' if sym else 'no symmetry'}, rounds from w = "
                       f"{run['start_weight']}")
            if run["lower_bound"] > d or run["found_weight"] not in (0, d):
                report(False, f"{lab}: lower bound {run['lower_bound']} or found weight {run['found_weight']} "
                              f"contradicts d = {d}")
            words = [w for _, w in read_nz(slug, "third_tools", run["words"])] if run.get("words") else []
            words += [[int(v) - 1 for v in m.group(2).split()] for m in CCW.finditer(text)]
            if run.get("n_words") is not None:
                report(len(words) == run["n_words"] and len({tuple(sorted(w)) for w in words}) == len(words),
                       f"{lab}: {len(words)} distinct words listed (runs.json: {run['n_words']})")
            bad = [w for w in words if len(w) != d or not all(code.logical(side, w))]
            if words:
                report(not bad, f"{lab}: all {len(words)} words are nontrivial {side}-logicals of weight {d}")
        else:
            rec = load(slug, "third_tools", run["record"])
            lab = f"{slug}: third_tools/{run['record']} ({tool}, side {side})"
            inst = {"Z": "HX", "X": "HZ"}[side]
            ok = (rec["status"] == "exact" and rec["lower_bound"] == rec["upper_bound"] == d
                  and rec["instance"] == inst and rec["threads"] == 1 and rec["cols"] == code.n
                  and rec["obs"] == code.k)
            report(ok, f"{lab}: status {rec['status']}, lower bound {rec['lower_bound']} = upper bound "
                       f"{rec['upper_bound']} (d = {d}), instance {rec['instance']}, one thread")
            code.check_word(side, rec["witness"], d, f"{lab}: its witness is a nontrivial {side}-logical")
    for side in "XZ":
        free = [r for r in T["runs"] if r["tool"] == "dist-m4ri" and r["side"] == side and not r["symmetry"]]
        base = max([r["lower_bound"] for r in free if r["start_weight"] == 1] + [0])
        best = max([base] + [r["lower_bound"] for r in free if 1 < r["start_weight"] <= base])
        report(best == d, f"{slug}: dist-m4ri without any symmetry: d_{side} >= {best} (a run whose rounds start at "
                          f"w = s counts only after a run from w = 1 proved d_{side} >= s); d = {d}")
        trop = [r for r in T["runs"] if r["tool"] == "tropidist" and r["side"] == side]
        if trop:
            report(all(r["status"] == "exact" and r["lower_bound"] == r["upper_bound"] == d for r in trop),
                   f"{slug}: tropidist, without symmetry: d_{side} = {d} exactly")


# ------------------------------------------------------------------------------------------------ entries

def entry_checks(fn):
    with open(fn, encoding="utf-8") as f:
        doc = json.load(f)
    slug = os.path.basename(fn)[:-5]
    if slug not in SLUGS:
        report(False, f"entry {fn}: {slug} is not a code of this package")
        return
    D = load(slug, "definition.json")
    P = D["parameters"]
    hx, hz = read_rows(slug, "hx.txt"), read_rows(slug, "hz.txt")
    code = Code(P["n"], hx, hz)
    ch = doc.get("checks", {})
    report(ch.get("X") == hx and ch.get("Z") == hz,
           f"entry {slug}: checks.X = hx.txt and checks.Z = hz.txt row by row")
    report(doc.get("code_type") == "CSS" and doc.get("n") == code.n and doc.get("k") == code.k,
           f"entry {slug}: CSS, n = {doc.get('n')}, k = {doc.get('k')} (here n = {code.n}, k = {code.k})")
    dist = doc.get("distance", {})
    report(dist.get("d") == P["d"], f"entry {slug}: distance.d = {dist.get('d')} (here {P['d']})")
    for side in "XZ":
        s = dist.get(side, {})
        report(s.get("value") == P[f"d_{side}"], f"entry {slug}: distance.{side}.value = {s.get('value')}")
        code.check_word(side, s.get("witness", []), P[f"d_{side}"],
                        f"entry {slug}: distance.{side}.witness is a nontrivial {side}-logical")
        if sorted(s.get("witness", [])) != sorted(read_word(slug, f"witness_{side}.txt")):
            info(f"entry {slug}: distance.{side}.witness differs from witness_{side}.txt (both checked)")


# ------------------------------------------------------------------------------------------------ main

def check_code(slug):
    D = load(slug, "definition.json")
    report(D["slug"] == slug, f"{slug}: definition.json describes {D['name']}")
    if D["kind"] == "pair-partition":
        template_checks(slug, D)
        n, hx, hz = build_l3(D)
    else:
        n, hx, hz = build_2bga(slug, D)
    report(hx == read_rows(slug, "hx.txt") and hz == read_rows(slug, "hz.txt")
           and read_text(slug, "hx.txt") == rows_text(hx) and read_text(slug, "hz.txt") == rows_text(hz),
           f"{slug}: H_X, H_Z built from definition.json equal hx.txt, hz.txt row by row ({len(hx)} + {len(hz)} rows)")
    code = Code(n, hx, hz)
    par = code_checks(slug, D, code)
    witness_checks(slug, D, code)
    record_checks(slug, D, code, par)
    third_checks(slug, D, code)
    P = D["parameters"]
    info(f"{slug}: [[{code.n},{code.k},{P['d']}]], w = {P['w']}, d_X = {P['d_X']}, d_Z = {P['d_Z']}, "
         f"kd^2/n = {code.k * P['d'] ** 2 / code.n:.2f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("slugs", nargs="*")
    ap.add_argument("--entry", action="append", default=[])
    a = ap.parse_args()
    slugs = a.slugs or ([] if a.entry else SLUGS)
    for s in slugs:
        if s not in SLUGS:
            report(False, f"unknown slug {s}")
            continue
        check_code(s)
    for fn in a.entry:
        entry_checks(fn)
    print("RESULT: all checks passed" if not FAILURES else f"RESULT: {len(FAILURES)} checks failed")
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
