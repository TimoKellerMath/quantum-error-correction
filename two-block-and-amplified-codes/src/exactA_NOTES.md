# exactA: exact search for CSS codes of any check weight (track T-E, implementation A)

`src/exactA.c` generalizes `vsearch_orb.c` (implementation A of
`quantum-error-correction/weight6-exact-distances/exact_distances.tex`, cited below as [ED]) to

1. arbitrary and mixed check weights $\le16$ (CSR instances; role masks of 16 bits);
2. any number $k$ of logical qubits (classes in $\lceil k/64\rceil$ words);
3. any number of qubits (no block structure);
4. anchors per qubit (trivial group) or per orbit of a group given by generators ([ED], Lemma 2.4), and
   pattern classes for such a group (prompt Section 8, C1 item 2; Lemma A below);
5. the parity reduction by one rank computation ($\mathbf 1\in\operatorname{rowsp}H_X$).

Written from vsearch_orb.c, [ED] and prompt C1 items 1 to 7 only (independence rule of track T-E: cluster2.c,
cluster3.c, vs2.c, vsearch.c and `new_codes/tools/exactB/` were not read).

## Usage

```
make -C new_codes/tools/exactA            # bin/exactA (gcc -O3 -march=native)
bin/exactA NAME.hx NAME.hz --side Z|X --W W [--group GENS] [--mode anchor|pattern] [--threads t]
    [--order index|bfs] [--maj 0|1|2] [--pack 0|1] [--tie 0|1] [--noparity] [--split 0|1]
    [--jobs a:b] [--tabmax w] [--enumerate FILE] [--keeporder] [-v]
```

* `NAME.hx`, `NAME.hz`: common format (first line `n m`, then `m` lines of 0-based qubits). Empty and repeated
  rows are dropped (they do not change $\ker$ or row spaces); the CSS condition is checked.
* Side `Z`: least weight of $z\in\ker H_X\setminus\operatorname{rowsp}H_Z$; checks are the rows of $H_X$, the
  majority family and the nontriviality test come from $H_Z$. Side `X` is the mirror image.
* `--group GENS`: first the number $g$ of generators, then $g$ permutations of $\{0,\dots,n-1\}$ ($n$ images
  each, image of qubit $q$ in position $q$). The program checks that every generator maps the set of rows of
  $H_X$ onto itself and the set of rows of $H_Z$ onto itself (by hashing the sorted supports) and derives the
  action on checks from the action on qubits; otherwise it stops. `src/inst.py aut NAME` writes `NAME.gens`
  for the automorphism group of the Tanner graph with the three vertex classes fixed (nauty's dreadnaut).
* `--mode anchor` (default): one job per qubit orbit ([ED], Lemma 2.4); without `--group`, one per qubit.
  `--mode pattern`: one job per pair (check orbit, pattern class) (Lemma A, Corollary B).
* `--W W`: searches weights $\le W$, or $\le W-1$ if $W$ is odd and all logicals of the side are even (parity,
  reported in the header); `--noparity` switches the reduction off.
* `--threads t`: work items (jobs, or with `--split 1` the iterations of the root call of Grow of each job)
  are distributed dynamically; every item runs from the empty state, so its search tree does not depend on
  the thread. Default `--split 1` if $t>1$ and there are fewer than $4t$ jobs.
* `--jobs a:b`: only jobs $a,\dots,b-1$ (a partial run proves nothing by itself; the verdict says so).
* `--enumerate FILE`: do not stop at the first word; write every reported word.
* `--order bfs`: anchors (and check orbits) by breadth-first order from qubit 0 instead of the index order;
  better for the trivial group and small groups (230-16-12: 80k against 170k nodes).
* `--keeporder`: the input order of the qubits of a check is its role order (default: sorted); used to
  reproduce vsearch_orb.c exactly.

Output: a header (weights, ranks, $k$, parity, group, jobs), then
`RESULT side= W= Weff= verdict=none|found weight= nodes= leaves= wall= cpu= ...`, for a found word
`WORD w : q_1 ... q_w` (rechecked: in $\ker$ of the checks by the check lists, outside the row space of the
dual matrix by its echelon form, which is independent of the class bits used in the search), and
`CONCLUSION d_Z >= ...` or `CONCLUSION d_Z <= w`. Nodes are calls of Grow, as in vsearch_orb.c.

## Correspondence with [ED]

Numbers are those of the compiled [ED]: Lemma 2.1 (minimal words), Lemma 2.2 (parity), Lemma 2.4 (orbit
anchors), Lemma 2.5 (patterns), Algorithm 1, Proposition 3.1 (exactness), Corollary 3.2 (runs). The brief
refers to "Lemma 2bgasym(b)"; that label is not in [ED] (presumably `gross_framework.tex`, not read here), and
the statement used for orbit anchors is [ED], Lemma 2.4, which is stated for any group of qubit permutations
preserving $\operatorname{rowsp}H_X$ and $\operatorname{rowsp}H_Z$.

| code (src/exactA.c) | statement | change for general weights |
|---|---|---|
| `leaf`: XOR of the class bits of $T$ | [ED], Section 2, the nontriviality test before Lemma 2.1 | none; $k$ bits in $\lceil k/64\rceil$ words |
| `propagate`, plain checks: fail iff no free qubit and $\lvert T\cap c\rvert$ odd; force the only free qubit | rule of Propagate (Algorithm 1) for $\mathcal P_c$ = all even subsets | Lemma D below (the table of vsearch_orb.c for $w=6$ and no forbidden pattern computes the same) |
| `propagate`, restricted checks (`rrule`): $3^w$ table ($w\le12$) or enumeration with the counting test; forced qubits applied in role order | rule of Propagate for $\mathcal P_c=\{\emptyset\}\cup\{P: \text{class}\ge\Pi\}$ | Lemmas C and E below |
| `propagate`, majority family | Algorithm 1, majority rules; Lemma 2.1(b) | none: $\mathcal S$ = rows of the dual matrix and sums of two rows sharing a qubit (pairs $(r,s)$, $r<s$, in the order of codes.py), threshold $\lfloor\lvert s\rvert/2\rfloor$, as in [ED], Section 3 |
| `lbound` | $\beta(T,X)$ of [ED], Section 3: $\max(\lceil u/h\rceil,\ \text{greedy packing})$ | none (no weight enters); $h$ from per-qubit counts of unsatisfied checks kept by `toggle_unsat`; both bounds stop once the test $\lvert T\rvert+\beta\le W$ is decided (Remark H) |
| `grow` | Grow of Algorithm 1; check with fewest free qubits, ties by smallest index, candidates in role order | none (the proof of Prop. 3.1 allows any unsatisfied check and any order) |
| `root_branch` (`--split 1`) | the iteration $i$ of the loop of the root call of Grow | Remark F below |
| `setup_job`, anchor mode: $T_0=\{q_i\}$, $X_0=O_1\cup\dots\cup O_{i-1}$ | Lemma 2.4 with Corollary 3.2 | none; orbits of the generated group, ordered by least index (or least BFS rank with `--order bfs`); trivial group = one anchor per qubit |
| `setup_job`, pattern mode | Lemma A and Corollary B below | replaces Lemma 2.5 (which needs a simply transitive action and one weight) |
| parity test (`ech_reduce` of $\mathbf 1$ by the echelon form of $H_X$) | Lemma 2.2 | hypothesis "odd column weights" replaced by $\mathbf 1\in\operatorname{rowsp}H_X$, same proof (Lemma G) |
| weight-one test before the jobs | Lemma 2.1(a) | qubits in no check (Lemma A, hypothesis) |

### What changes in the proof

Nothing in Lemma 2.1, Lemma 2.4, Proposition 3.1 or Corollary 3.2 of [ED] uses the check weight or a block
structure: Proposition 3.1 holds for any families $\mathcal P_c$ of even sets containing $\emptyset$, any
$\mathcal S\subseteq\operatorname{rowsp}H_Z$ and any lower bound $\beta$ for $\lvert z\setminus T\rvert$. The weight
enters only (i) in Lemma 2.5, replaced by Lemma A, (ii) in the hypothesis of Lemma 2.2, replaced by Lemma G,
and (iii) in the implementation of the rule of Propagate, justified by Lemmas C, D, E.

Throughout: $H_X\in\mathbb F_2^{m_X\times n}$, $H_Z\in\mathbb F_2^{m_Z\times n}$, $H_XH_Z^{\mathsf T}=0$; side $Z$
(side $X$ is symmetric); $\mathcal C$ is the set of supports of the rows of $H_X$ (the checks).

**Lemma A (pattern classes; prompt C1 item 2).** Let $\Gamma$ be a finite group acting on $\{0,\dots,n-1\}$
such that every $\gamma\in\Gamma$ maps $\mathcal C$ onto $\mathcal C$ and $\operatorname{rowsp}H_Z$ onto
itself; let $\Gamma$ act on $\mathcal C$ by $c\mapsto\gamma(c)$. Let $c_1,\dots,c_r$ represent the
$\Gamma$-orbits on $\mathcal C$, $S_i=\operatorname{Stab}_\Gamma(c_i)$, and for $c\in\Gamma c_i$ fix
$\gamma_c\in\Gamma$ with $\gamma_cc_i=c$. For $z\subseteq\{0,\dots,n-1\}$ and $c\in\Gamma c_i$ let
$\operatorname{cl}_z(c)$ be the $S_i$-orbit of $\gamma_c^{-1}(z\cap c)$, a set of subsets of $c_i$.
1. $\operatorname{cl}_z(c)$ does not depend on the choice of $\gamma_c$.
2. $\operatorname{cl}_{\gamma z}(\gamma c)=\operatorname{cl}_z(c)$ for all $\gamma\in\Gamma$, $c\in\mathcal C$.
3. Fix a total order on the pairs $(i,\Pi)$, $\Pi$ an $S_i$-orbit of nonempty even subsets of $c_i$, in which
   $(i,\Pi)<(i',\Pi')$ whenever $i<i'$, and fix $\pi_\Pi\in\Pi$. Assume $d_Z<\infty$ and that no $\{q\}$ with
   $q$ in no check is a nontrivial $Z$-logical operator. Then there are a pair $(i,\Pi)$ and a nontrivial
   $Z$-logical operator $z$ with $\lvert z\rvert=d_Z$, $z\cap c_i=\pi_\Pi$, $z\cap c=\emptyset$ for every check $c$
   in an orbit $\Gamma c_{i'}$ with $i'<i$, and $z\cap c=\emptyset$ or $\operatorname{cl}_z(c)\ge\Pi$ (in the order of
   the pairs with first entry $i$) for every $c\in\Gamma c_i$.

*Proof.* (1) Another choice is $\gamma_cs$ with $s\in S_i$, and $(\gamma_cs)^{-1}(z\cap c)=s^{-1}\gamma_c^{-1}(z\cap c)$
lies in the same $S_i$-orbit. (2) $\gamma\gamma_c$ maps $c_i$ to $\gamma c$, and $\gamma z\cap\gamma c=\gamma(z\cap c)$,
so by (1) $\operatorname{cl}_{\gamma z}(\gamma c)$ is the $S_i$-orbit of $(\gamma\gamma_c)^{-1}\gamma(z\cap c)=\gamma_c^{-1}(z\cap c)$.
(3) Let $z_0$ be a nontrivial $Z$-logical operator of weight $d_Z$. If $z_0$ contains a qubit $q$ in no check,
then $\{q\}\in\ker H_X$, so $z_0=\{q\}$ by [ED], Lemma 2.1(a), which the hypothesis excludes. Hence every qubit
of $z_0\ne\emptyset$ lies in a check, and the set $P(z_0)$ of pairs $(i(c),\operatorname{cl}_{z_0}(c))$ over the
checks $c$ with $z_0\cap c\ne\emptyset$ is nonempty; here $i(c)$ is the index of the orbit of $c$, and the class
consists of nonempty even sets because $z_0\in\ker H_X$. Let $(i,\Pi)$ be its least element, attained at $c$,
and $s\in S_i$ with $s\gamma_c^{-1}(z_0\cap c)=\pi_\Pi$. Put $\gamma=s\gamma_c^{-1}$ and $z=\gamma z_0$. As $\gamma$
permutes $\mathcal C$, it preserves $\ker H_X$; it preserves $\operatorname{rowsp}H_Z$ and weights, so $z$ is a
nontrivial $Z$-logical operator of weight $d_Z$. Since $\gamma c=c_i$, $z\cap c_i=\gamma(z_0\cap c)=\pi_\Pi$. By (2),
and because $\gamma$ permutes $\mathcal C$ preserving orbits with $z\cap\gamma c'=\gamma(z_0\cap c')$, $P(z)=P(z_0)$,
whose least element is $(i,\Pi)$; this gives the last two assertions. $\square$

For $r=1$ and $S_1=1$ this is [ED], Lemma 2.5. For the trivial group, the orbits are the single checks and
the classes the single masks (a "first check met" decomposition).

**Corollary B.** In the situation of Lemma A, for each pair $(i,\Pi)$ let $R(i,\Pi)$ be Algorithm 1 of [ED] with
$T_0=\pi_\Pi$, $X_0=(c_i\setminus\pi_\Pi)\cup\bigcup_{i'<i}\bigcup_{c\in\Gamma c_{i'}}c$, $\mathcal P_c=\{\emptyset\}\cup
\{P\subseteq c:\lvert P\rvert\text{ even},\ P\ne\emptyset,\ \operatorname{cl}_P(c)\ge\Pi\}$ for $c\in\Gamma c_i$ and $\mathcal P_c$
the even subsets of $c$ otherwise. If no qubit in no check is a nontrivial $Z$-logical operator and no run
$R(i,\Pi)$ reports a set, then $d_Z\ge W+1$.

*Proof.* Otherwise $d_Z\le W$, and Lemma A(3) gives $(i,\Pi)$ and $z$. Then $T_0\subseteq z$; $z\cap X_0=\emptyset$
because $z\cap c_i=\pi_\Pi$ and $z$ meets no check of an earlier orbit; and $z\cap c\in\mathcal P_c$ for every check
$c$ (even because $z\in\ker H_X$; for $c\in\Gamma c_i$ empty or of class $\ge\Pi$). By [ED], Proposition 3.1,
$R(i,\Pi)$ reports a set. $\square$

In the program, $\gamma_c$ is the element of the Schreier tree of the orbit (breadth-first search from $c_i$
with the generators), the qubits of $c$ are stored in the role order $\rho_c=\gamma_c\circ\rho_{c_i}$ ($\rho_{c_i}$
the sorted order of $c_i$), so the role mask of $z$ at $c$ is $\rho_{c_i}^{-1}\gamma_c^{-1}(z\cap c)$. The action of
$S_i$ on roles is generated by the Schreier generators $\gamma_{gc}^{-1}g\gamma_c$ ($c\in\Gamma c_i$, $g$ a generator),
which generate $S_i$ (Schreier's lemma; e.g. Seress, *Permutation group algorithms*, Cambridge Univ. Press
2003, Lemma 4.2.1, number unverified); on roles, such an element maps $r$ to the position of $g(\rho_c(r))$ in
$\rho_{gc}$. The classes are the connected components of the nonempty even masks under these permutations
(union-find). Orbits are taken in order of their least check index (`--order bfs`: least BFS rank), classes
by decreasing weight, then increasing least mask (the order `num` of vsearch_orb.c), and $\pi_\Pi$ is the least
mask of $\Pi$. $\mathcal P_c=\{\emptyset\}$ on earlier orbits is realized by $X_0$, which makes it automatic.

**Lemma C (counting test).** Let $c$ be a check with role set $R$, $F$ a set of nonempty even subsets of $R$,
$\mathcal P_c=\{P\subseteq R:\lvert P\rvert\text{ even},\ P\notin F\}$, $I,O\subseteq R$ disjoint, $\Phi=R\setminus(I\cup O)$ and
$f=\lvert\Phi\rvert\ge2$. If $2^{f-2}>\lvert F\rvert$, the rule of Propagate for $c$ (with $T\cap c=I$, $X\cap c=O$) neither
fails nor puts a qubit into $T$ or $X$.

*Proof.* For $r\in\Phi$ there are $2^{f-2}$ even $P$ with $I\subseteq P\subseteq I\cup\Phi$ and $r\in P$ (the parity
of $P\cap(\Phi\setminus\{r\})$ is fixed), and $2^{f-2}$ with $r\notin P$; at most $\lvert F\rvert$ of each kind are
forbidden. So both kinds occur in $\mathcal P_c$: some compatible $P$ exists, and $r$ lies in some but not in all. $\square$

**Lemma D (plain checks).** For $\mathcal P_c$ the even subsets of $c$ the rule of Propagate fails iff $\Phi=\emptyset$
and $\lvert I\rvert$ is odd, puts the unique free qubit into $T$ (resp. $X$) iff $\lvert\Phi\rvert=1$ and $\lvert I\rvert$ is odd
(resp. even), and does nothing otherwise. *Proof.* Lemma C with $F=\emptyset$ for $f\ge2$; for $f\le1$ the only
compatible even set is $I\cup\Phi$ or $I$, according to the parity of $\lvert I\rvert$. $\square$

**Lemma E (table).** For $w\le12$ the program stores, for every pair $(I,O)$ of disjoint role sets (ternary
index $\sum_{r\in I}3^r+2\sum_{r\in\Phi}3^r$), whether a compatible $P\in\mathcal P_c$ exists, and the intersection
and union of all of them. The recursion on the least free role $j$ is correct because the compatible sets
for $(I,O)$ are the disjoint union of those for $(I,O\cup\{j\})$ and $(I\cup\{j\},O)$. The rule then puts
(intersection)$\setminus I$ into $T$ and $\Phi\setminus$(union) into $X$, which is the rule of Algorithm 1. For $w>12$
(or `--tabmax`) the program applies Lemma C and otherwise enumerates the $2^{f-1}$ compatible even sets; the
tests compare both paths (identical node counts).

**Remark F (split items).** With `--split 1` the item $(J,b)$ applies $T_0,X_0$ of job $J$, propagates, takes the
same branching check and candidates $q_1,\dots,q_s$ as the root call of Grow, and replays the iterations
$i<b$ (skip $q_i\in X$; stop if $q_i\in T$, where the original loop calls Grow and breaks; otherwise put $q_i$
into $X$ and stop if Propagate fails), then performs iteration $b$. The original loop reaches iteration $b$
iff this replay does not stop, and then in the same state (as sets; the order of the list of unsatisfied
checks may differ, which changes only the value of the packing bound, a valid lower bound either way). The
induction in the proof of Proposition 3.1 applies to the item $b=j$, $q_j$ the first candidate in $z$:
compatible states never fail. Without the packing bound the node counts of split and unsplit runs agree
exactly (tested).

**Remark H (early decision in $\beta$).** Propagate uses $\beta$ only in the test $\lvert T\rvert+\beta\le W$. With
budget $b=W-\lvert T\rvert\ge1$ the program returns as soon as the greedy packing exceeds $b$, and otherwise returns
a value $>b$ iff no free qubit of an unsatisfied check lies in $\lceil u/b\rceil$ or more unsatisfied checks, i.e.
iff $\lceil u/h\rceil>b$. So the test, and hence the search tree, is that of the full computation; node counts were
compared on all validation runs (identical).

**Lemma G (parity).** Every element of $\ker H_X$ has even weight iff $\mathbf 1\in\operatorname{rowsp}H_X$.
*Proof.* $\ker H_X\subseteq\mathbf 1^\perp$ iff $\mathbf 1\in(\ker H_X)^\perp=\operatorname{rowsp}H_X$. $\square$
So with parity, no word of weight $\le W$ (searched up to the largest even number $\le W$) gives $d_Z\ge W+1$,
and $d_Z\ge W+2$ if $W$ is even. (If $\mathbf 1\notin\operatorname{rowsp}H_X$, the nontrivial words may still all
have one parity, e.g. when $k=1$ and $\operatorname{rowsp}H_Z$ is even; the program does not use such cases.)

## Files

* `src/exactA.c`: the program (one file, C11, pthreads).
* `src/inst.py`: board code to instance (`board`), seeded relabelling (random qubit permutation, random order
  of the checks and of the qubits within each check; `inst/SLUG_rS.perm` is the permutation, `.wX`, `.wZ` the
  mapped board witnesses), Tanner-graph group by dreadnaut (`aut`, with a check that every generator maps the
  rows of $H_X$ and $H_Z$ onto themselves), independent GF(2) checks (`check`, `info`).
* `src/test_small.py`: tests against references written independently of exactA (kernel enumeration for
  random codes with mixed weights 2 to 10; one MILP per dual logical, scipy/HiGHS, for hypergraph products,
  both sides; disjoint unions with $k=81$).
* `src/validate.py`: validation runs on the board codes and the growth fits.
* `src/cmp_vsearch.py`: regression test against vsearch_orb.c (same codes, same role orders).
* `src/cmp_logs.py`: run-by-run comparison of verdicts and node counts of two log directories.
* `inst/`: instances (`SLUG.{hx,hz}` as on the board, `SLUG_r1`, `SLUG_r2` relabelled with seeds 1, 2, with
  `.gens`; the twelve codes of C1 item 6 with both seeds, the six mixed-weight codes with seed 1),
  `inst/cmp/` for the comparison with vsearch_orb (`cmp_vsearch.py` does not write bytecode into the package).
* `bin/` (git-ignored): `exactA`; `exactA_round1`, `exactA_round2` (first and second version, for the logs);
  `vsearch_orb` built from the package for the regression test.
* `Makefile`: `make`, `make test`, `make instances`, `make validate`, `make summary`.
* Logs: `new_codes/logs/T-E/exactA/` (`test_full/`, `test_quick/`; `val/` round 1 with the first version;
  `val2/`, `val_r2/`, `val_mixed/` with the second version; `val3*` with the final version; `cmp/`).

## Validation protocol (prompt C1 item 6)

For each board code (snapshot `8042f9da`, certified $d$ in `certs/SLUG.json`, tropidist except `144-2-12`,
CryptoMiniSat) and each side: the relabelled instance; $W=d-2$ if all logicals of the side are even (parity,
Lemma G), else $W=d-1$: expect `none`; $W=d$: expect a word, which `validate.py summary` rechecks with its own
GF(2) code (`inst.py`, Python integers, independent of exactA). Modes: P = Tanner-graph group with pattern
classes, A = Tanner-graph group with orbit anchors, T = trivial group (one anchor per qubit, BFS order).

## Results (computed 2026-11-10)

Final version: `src/exactA.c` sha256 `2d1f0d1a55d9b9d71eff102eab6f131b4a6ff2da0e973518651877c3c3187856`, built `gcc -O3 -march=native` (gcc 16.2.1); single
thread per run, two runs at a time, on the shared 8-core workstation. Every run below gave the expected verdict;
every word found at $W=d$ has weight $d$ and was rechecked by `inst.py` (in $\ker$, nontrivial). Runs: 236 with
the first version (`val/`), 348 with the second (`val2/`, `val_r2/`, `val_mixed/`), the same 348 with the final
version (`val3*`: identical verdicts and node counts, run by run, `src/cmp_logs.py`), plus the test suites
(`test_full/`: 45 instance-sides, `test_quick/`: 14, no failure). CPU of the track: about 1.9 hours.

**Board codes, relabelled with seed 1** (cell: verdict, nodes, CPU seconds). ¹ first version, under load 13 to
15 (verdict and nodes of the first version; its node counts differ from the final version in a few runs by
at most 0.04 % through the order of the greedy packing). The trivial group on `234-24-18` was skipped
(extrapolated above 30 CPU-minutes per side).

| code | w | side | W | P: patterns, Tanner group | A: orbit anchors, Tanner group | T: one anchor per qubit | W = d (P / A / T) |
|---|---|---|---|---|---|---|---|
| 144-2-12 | 4 | Z | 10 (parity) | none 280 (0.00 s) | none 536 (0.00 s) | none 6 454 (0.00 s) | 12 / 12 / 12 |
| 144-2-12 | 4 | X | 10 (parity) | none 280 (0.00 s) | none 536 (0.00 s) | none 6 510 (0.00 s) | 12 / 12 / 12 |
| 240-4-18 | 5 | Z | 17 | none 28 593 (0.03 s) | none 103 177 (0.09 s) | none 1 181 361 (1.02 s) | 18 / 18 / 18 |
| 240-4-18 | 5 | X | 17 | none 24 389 (0.02 s) | none 85 399 (0.08 s) | none 1 252 387 (1.06 s) | 18 / 18 / 18 |
| 780-64-10 | 5 | Z | 9 | none 2 494 (0.01 s) | none 3 076 (0.00 s) | none 21 436 (0.02 s) | 10 / 10 / 10 |
| 780-64-10 | 5 | X | 9 | none 2 555 (0.01 s) | none 3 106 (0.00 s) | none 21 444 (0.02 s) | 10 / 10 / 10 |
| 248-10-18 | 6 | Z | 16 (parity) | none 107 344 (0.14 s) | none 314 057 (0.40 s) | none 3 404 418 (3.95 s) | 18 / 18 / 18 |
| 248-10-18 | 6 | X | 16 (parity) | none 100 863 (0.14 s) | none 329 147 (0.40 s) | none 3 379 273 (3.93 s) | 18 / 18 / 18 |
| 224-12-16 | 6 | Z | 14 (parity) | none 15 053 (0.02 s) | none 30 860 (0.04 s) | none 451 245 (0.52 s) | 16 / 16 / 16 |
| 224-12-16 | 6 | X | 14 (parity) | none 16 441 (0.02 s) | none 37 610 (0.05 s) | none 452 710 (0.52 s) | 16 / 16 / 16 |
| 240-16-18 | 7 | Z | 17 | none 657 726 (1.32 s) | none 2 112 443 (3.97 s) | none 20 965 380 (118.31 s)¹ | 18 / 18 / 18 |
| 240-16-18 | 7 | X | 17 | none 676 304 (1.38 s) | none 1 898 002 (3.57 s) | none 21 045 834 (150.93 s)¹ | 18 / 18 / 18 |
| 224-18-16 | 7 | Z | 15 | none 52 768 (0.12 s) | none 164 037 (0.33 s) | none 2 367 532 (4.14 s) | 16 / 16 / 16 |
| 224-18-16 | 7 | X | 15 | none 64 181 (0.14 s) | none 194 163 (0.37 s) | none 2 247 572 (3.92 s) | 16 / 16 / 16 |
| 248-24-18 | 8 | Z | 17 | none 2 237 190 (6.61 s) | none 7 871 358 (19.60 s) | none 80 397 749 (636.53 s)¹ | 18 / 18 / 18 |
| 248-24-18 | 8 | X | 17 | none 2 194 541 (6.27 s) | none 7 885 453 (20.02 s) | none 83 260 476 (603.26 s)¹ | 18 / 18 / 18 |
| 232-62-12 | 8 | Z | 10 (parity) | none 9 332 (0.02 s) | none 11 274 (0.02 s) | none 166 834 (0.27 s) | 12 / 12 / 12 |
| 232-62-12 | 8 | X | 10 (parity) | none 9 581 (0.03 s) | none 11 266 (0.02 s) | none 165 075 (0.26 s) | 12 / 12 / 12 |
| 234-24-18 | 9 | Z | 17 | none 6 408 570 (26.16 s) | none 28 098 730 (94.91 s) | skipped | 18 / 18 / - |
| 234-24-18 | 9 | X | 17 | none 6 288 720 (22.86 s) | none 28 093 352 (100.54 s) | skipped | 18 / 18 / - |
| 254-42-14 | 10 | Z | 12 (parity) | none 21 392 (0.17 s) | none 51 940 (0.22 s) | none 682 612 (2.67 s) | 14 / 14 / 14 |
| 254-42-14 | 10 | X | 12 (parity) | none 22 922 (0.17 s) | none 54 097 (0.22 s) | none 698 742 (2.80 s) | 14 / 14 / 14 |
| 672-336-12 | 12 | Z | 10 (parity) | none 45 524 (3.93 s) | none 53 358 (0.18 s) | none 2 360 692 (7.59 s) | 12 / 12 / 12 |
| 672-336-12 | 12 | X | 10 (parity) | none 40 463 (3.93 s) | none 54 977 (0.19 s) | none 2 353 327 (7.62 s) | 12 / 12 / 12 |

Seed 2 (another random relabelling and its own Tanner group), primary mode only (P for $w\le9$, A above), both
sides, $W=d-1$ or $d-2$ and $W=d$: 48 runs, all as expected (`val3_r2/summary.md`).

**Certified board codes with mixed check weights**, seed 1 (no parity on either side; `230-16-12` has two
repeated $X$-rows, dropped, and its Tanner group acts trivially on the qubits, so P and A use the trivial group
in index order; `196-15-11` has the trivial group; `172-48-8` has checks of weight 13, above the table limit,
so its pattern runs use Lemma C):

| code | weights | side | W | P | A | T | W = d (P / A / T) |
|---|---|---|---|---|---|---|---|
| 132-32-8 | 7-11 | Z | 7 | none 1 127 (0.24 s) | none 1 505 (0.00 s) | none 4 838 (0.01 s) | 8 / 8 / 8 |
| 132-32-8 | 7-11 | X | 7 | none 632 (0.24 s) | none 1 503 (0.00 s) | none 5 050 (0.01 s) | 8 / 8 / 8 |
| 172-48-8 | 9-13 | Z | 7 | none 859 (0.37 s) | none 1 160 (0.00 s) | none 7 172 (0.02 s) | 8 / 8 / 8 |
| 172-48-8 | 9-13 | X | 7 | none 1 514 (0.37 s) | none 1 185 (0.01 s) | none 8 817 (0.02 s) | 8 / 8 / 8 |
| 210-56-8 | 7-10 | Z | 7 | none 1 083 (0.11 s) | none 1 827 (0.00 s) | none 5 294 (0.01 s) | 8 / 8 / 8 |
| 210-56-8 | 7-10 | X | 7 | none 1 996 (0.11 s) | none 2 486 (0.01 s) | none 6 160 (0.01 s) | 8 / 8 / 8 |
| 230-16-12 | 2-8 | Z | 11 | none 113 249 (0.29 s) | none 170 479 (0.31 s) | none 80 017 (0.15 s) | 12 / 12 / 12 |
| 230-16-12 | 2-8 | X | 11 | none 94 252 (0.26 s) | none 135 993 (0.23 s) | none 84 529 (0.15 s) | 12 / 12 / 12 |
| 196-15-11 | 2-8 | Z | 10 | none 48 779 (0.17 s) | none 45 866 (0.08 s) | none 30 960 (0.05 s) | 11 / 11 / 11 |
| 196-15-11 | 2-8 | X | 10 | none 42 707 (0.15 s) | none 77 857 (0.13 s) | none 40 123 (0.07 s) | 11 / 11 / 11 |
| 127-1-15 | 4, 8 | Z | 14 | none 700 007 (0.53 s) | none 912 847 (0.65 s) | none 684 223 (0.45 s) | 15 / 15 / 15 |
| 127-1-15 | 4, 8 | X | 14 | none 761 474 (0.53 s) | none 941 771 (0.64 s) | none 667 670 (0.41 s) | 15 / 15 / 15 |

**Growth (prompt C1 item 7, T13).** Side $Z$, primary mode, `--noparity`, node counts at $W=d-8,\dots,d-1$,
least-squares fit of $\log N=\log N_0+W\log\beta$ (zero counts dropped). Computed:

| code | w | mode | beta, W=d-8..d-1 | beta, even W | beta, W=d-4..d-1 | N(d-4), ..., N(d-1) |
|---|---|---|---|---|---|---|
| 144-2-12 | 4 | P | 2.17 | 2.26 | 2.02 | 69, 139, 280, 567 |
| 240-4-18 | 5 | P | 2.23 | 2.22 | 2.31 | 2322, 5320, 12385, 28593 |
| 780-64-10 | 5 | P | 2.55 | 2.75 | 2.23 | 225, 517, 1136, 2494 |
| 248-10-18 | 6 | P | 2.62 | 2.61 | 2.65 | 15346, 40568, 107344, 284044 |
| 224-12-16 | 6 | P | 2.65 | 2.62 | 2.74 | 2046, 5424, 15053, 42008 |
| 240-16-18 | 7 | P | 2.84 | 2.83 | 2.91 | 26605, 76699, 222845, 657726 |
| 224-18-16 | 7 | P | 2.88 | 2.79 | 3.09 | 1781, 5426, 16710, 52768 |
| 248-24-18 | 8 | P | 3.42 | 3.40 | 3.53 | 50772, 180549, 632813, 2237190 |
| 232-62-12 | 8 | P | 3.99 | 3.70 | 3.88 | 621, 2415, 9332, 36275 |
| 234-24-18 | 9 | P | 3.53 | 3.48 | 3.75 | 121991, 447354, 1674190, 6408570 |
| 254-42-14 | 10 | A | 3.88 | 3.85 | 4.04 | 3106, 11098, 51940, 195569 |
| 672-336-12 | 12 | A | 4.54 | 4.29 | 5.26 | 1910, 10258, 53358, 278793 |

Heuristic reading (two to three codes per weight, $n\le780$, $d\le18$): $\beta(w)$ grows about linearly, roughly
$0.4w+0.4$ on the upper half ($2.0$ at $w=4$; $2.2$ to $2.3$ at $5$; $2.65$ to $2.75$ at $6$; $2.9$ to $3.1$ at $7$;
$3.5$ to $3.9$ at $8$; $3.75$ at $9$; $4.0$ at $10$; $5.3$ at $12$), and depends on the code, not only on $w$. At $w=6$
this is above the factor $2.45$ per unit of $W$ that prompt C3 extrapolates from runs on $[[756,16]]$,
$[[558,14,28]]$ and $[[600,24,26]]$ ($d=26$ to $34$); whether $\beta$ decreases with $n$ or $d$ is not tested here.

**Mode comparison** (same code, side, $W$). Patterns need 1.2 to 4.5 times fewer nodes than orbit anchors for
$w\le10$ (prompt C3 quotes about 3, from Remark ae:rem:variants of another tex, not [ED]); at $w=12$
(`672-336-12`) 1.2 to 1.4 times fewer nodes but about 20 times more time, spent building the $3^{12}$ tables of
4094 jobs. One anchor per qubit costs 7 to 44 times the nodes of the orbit anchors of the Tanner group; with the trivial or a small group the BFS order beats the index
order (`230-16-12`: 80k against 170k nodes; `127-1-15`: 684k against 913k).

**Regression against vsearch_orb.c** (`src/cmp_vsearch.py`, `logs/.../cmp/summary_final.txt`): on the gross code
$[[144,12,12]]$ at $W=10,11,12$ and on B510 and C504 of [ED] at $W=18,20$, `exactA --mode pattern --keeporder`
with the translation group and `vsearch_orb -K num` give identical node counts (up to $2.6\cdot10^6$); the
final version needs 3.50 s against 6.53 s (B510, $W=20$) and 2.28 s against 4.20 s (C504), back to back.

**Throughput** of the final version (nodes per CPU-second, load average 4 to 7): 0.5 to 0.86 million at
$w=6,7$; 0.34 to 0.4 million at $w=8$; 0.24 to 0.31 million at $w=9$ to $12$. The first version was 3.6 to 5.6 times
slower per node in round 1. The incremental counts gave 1.3 to 1.7 in spot checks and the early decision 1.7 over
the whole validation set (`val2` against `val3`); the rest is the heavier machine load of round 1 (not separated).

## Limits

1. Check weight $\le16$ on the searched side (role masks of 16 bits); the dual side may have any weight
   (it enters only through the majority family and the dense echelon form). Restricted checks of the
   pattern mode use a $3^w$ table for $w\le12$ (rebuilt per job: at $w=12$ and 4094 jobs about 4 s, which
   dominates `672-336-12`; use the anchor mode there), and Lemma C with enumeration above.
2. CSS codes only (prompt C1 item 5: non-CSS codes need a symplectic search, not implemented).
3. Dense GF(2) eliminations: memory $n^2/8$ bytes per echelon form, time $O(n^3/64)$; fine to $n\approx2\cdot10^4$.
4. The group must map the rows of $H_X$ and of $H_Z$ onto themselves (checked); a group that only preserves
   row spaces is refused. The X/Z exchange of [ED], Lemma 2.3, is not used: both sides are searched.
5. Pattern mode: $\sum_i(\text{number of }S_i\text{-classes of nonempty even masks})$ jobs, i.e. $m(2^{w-1}-1)$ for the
   trivial group; useful with a group with few check orbits. With the trivial group the anchor mode is the
   sensible choice (one job per qubit).
6. Cost per node grows with $w$ (larger majority family, more checks per qubit), and $\beta(w)$ grows about
   linearly (Results). At $w=8$, $\beta\approx3.5$: each unit of $W$ costs 3.5 times more, so $d\approx26$ at $w=8$
   is about $3.5^{8}\approx2\cdot10^4$ times the $d=18$ run of `248-24-18` (7 s), i.e. some 40 CPU-hours
   (heuristic extrapolation).
7. Parallelism is over jobs, or with `--split 1` over the iterations of the root call of Grow (at most $w$ items
   per job); no deeper splitting, so one expensive job bounds the wall time.
8. With the stop at the first word, a found word has weight $\le W$, not necessarily $d$; the distance needs a
   complete run at $W=d-1$ (or $d-2$ with parity) and a word of weight $d$.
9. Not implemented from C1: block parities (item 3, e.g. $\lvert z\cap L\rvert$ even for type $3+4$) as extra
   constraints in Propagate and $\beta$; second-level pattern splitting of the first job (C3).
10. The result is a solver verdict: correctness rests on [ED], Proposition 3.1, Lemma A, Lemmas C to G and on the
   program; the second implementation of track T-E (exactB) is the independent check.
