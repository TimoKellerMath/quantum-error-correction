# exactB: exact search at every check weight (Track T-E, implementation B)

`exactB.c` generalizes `src/cluster2.c` of `quantum-error-correction/weight6-exact-distances` (implementation B of the note `exact_distances.tex`, cited as **ED** below) to

1. arbitrary and mixed check weights (up to 64; validated up to 12),
2. any number $k$ of logical qubits (classes are $k$-bit vectors in $\lceil k/64\rceil$ words),
3. any number of qubits (no block structure is assumed),
4. anchors at every qubit, or at one qubit per orbit of a supplied permutation group, with the earlier orbits excluded,
5. the parity reduction, decided by one rank computation ($\mathbf 1\in\operatorname{rowsp}H_X$).

It reads the common instance format, computes its own ranks and logical basis, verifies group generators, rechecks every word it reports, and can split the search over threads. It shares no code with other programs of the project; the search core (state, propagation, branching) follows `cluster2.c`.

## Interface

```
make                                   # bin/exactB (gcc -O3 -march=native)
bin/exactB NAME.hx NAME.hz --side Z|X --W W [--group GENS] [--threads t]
      [--order id|bfs] [--orbit-order min|size] [--split D] [--bound 1|2|3]
      [--majority 0|1|2] [--stab-wmax w] [--run r] [--node-limit N] [--out FILE] [-v]
```

* `NAME.hx`, `NAME.hz`: first line `n m`, then $m$ lines with the 0-based qubits of one check (repeated qubits and out-of-range indices are rejected).
* `--side Z`: search $z\in\ker H_X\setminus\operatorname{rowsp}H_Z$ of weight $\le W$; `--side X` swaps the matrices. Below, $A$ is the matrix whose kernel holds the logicals of the chosen side ($A=H_X$ for side $Z$) and $B$ the other one.
* `--group GENS`: one permutation of $0,\dots,n-1$ per line in image form ($p(0)\ \dots\ p(n-1)$; lines starting with `#` are comments). Without it, every qubit is an anchor (the trivial group), in the order `id` or `bfs` (breadth first in the qubit graph of $A$).
* Output: `RESULT side= W= searched= found=yes|none nodes= generated= leaves= wall= cpu= ...`; if a word is found, `WORD weight=... syndrome=zero rowsp-test=nontrivial : q_1 ... q_w`, rechecked without the logical basis (syndrome recomputed from the rows of $A$; non-membership in $\operatorname{rowsp}B$ by reduction against the RREF of $B$); `CONCLUSION: d_Z >= ...` or `d_Z <= w`. Exit code 3 if a recheck fails (never observed).
* `nodes` counts calls of GROW (ED, Algorithm alg:search); `generated` counts children tried (set_in plus propagation), the node count of B in ED's table.
* `--run r` runs only the $r$-th anchor run (the conclusion is then not a bound). `--node-limit` stops early and marks the result incomplete.

Scripts: `scripts/relabel.py` (board JSON to a relabelled instance: seeded random qubit permutation, random check order and order within lines; also writes the permutation and the board witnesses after relabelling), `scripts/autgroup.m` (Magma: generators of the automorphism group of the Tanner graph with the three vertex classes fixed, restricted to the qubits), `scripts/vrun.sh` (one validation run with a summary line), `scripts/batch1.sh` to `scripts/batch4.sh` (the validation batches), `scripts/growth.py` (fits), `scripts/table.py` (validation tables, with a recheck of every reported word against the board matrices). Tests: `tests/gen_tests.py`, `tests/gen_odd.py` (small codes with distances by enumeration, written independently of exactB), `tests/run_tests.sh`.

## Correspondence with ED

| part of exactB | statement of ED | change for general weights |
|---|---|---|
| `search` (GROW), `propagate` (PROPAGATE) | Algorithm alg:search; completeness: Proposition prop:exact | none; the proof of prop:exact uses no weight |
| first rule of PROPAGATE with $\mathcal P_c$ = all even subsets | Algorithm alg:search, Corollary cor:runs (anchor case) | reduces to XOR unit propagation: Lemma N1 |
| a zero-syndrome $T$ is a leaf, reported iff nontrivial | Lemma lem:minimal(a) | none |
| majority rule over $\mathcal S$ = rows of $B$ and sums of two rows of $B$ with a common qubit | Lemma lem:minimal(b) | none (any subset of $\operatorname{rowsp}B$ is allowed); default weight cap $2\max_r|B_r|$, i.e. all such sums |
| branching on an unsatisfied check, $\beta$ | Lemma lem:minimal(c) | $\beta$ = max of the degree-sequence bound (Lemma N2, which dominates $\lceil u/h\rceil$ of ED) and a greedy packing (ED); greedy order (free count, index) |
| anchors and exclusions | Lemma lem:anchors, Corollary cor:runs | orbits of a verified group (Lemma N5) or singletons |
| parity | Lemma lem:parity | hypothesis replaced by $\mathbf 1\in\operatorname{rowsp}A$ (Lemma N3) |
| nontriviality test | ED, Section 2, first paragraph | classes in $\lceil k/64\rceil$ words, tested at leaves only (Lemma N6) |
| threads | (not in ED) | Lemma N4 |

The pattern anchoring (ED, Lemma lem:patterns and the pattern-class lemma of the prompt, Section 8, C1, item 2) is not implemented in exactB.

## Statements used beyond ED

Throughout, $A\in\mathbb F_2^{m_A\times n}$, $B\in\mathbb F_2^{m_B\times n}$ with $AB^{\mathsf T}=0$; a logical is $z\in\ker A$, trivial iff $z\in\operatorname{rowsp}B$, and $d$ is the least weight of a nontrivial one. A state is a pair of disjoint sets $(T,X)$ of qubits; the free qubits are the others; a check (row of $A$) $c$ is unsatisfied if $|T\cap c|$ is odd.

**Lemma N1.** Let $c$ be a check, $\mathcal P_c$ the family of even subsets of $c$, $(T,X)$ a state and $F=c\setminus(T\cup X)$. The sets $P\in\mathcal P_c$ with $T\cap c\subseteq P\subseteq c\setminus X$ are the sets $(T\cap c)\cup G$ with $G\subseteq F$ and $|G|\equiv|T\cap c|\pmod 2$. Hence: (i) if $F=\emptyset$, such $P$ exists iff $c$ is satisfied; (ii) if $F=\{q\}$, $q$ lies in every such $P$ if $c$ is unsatisfied and in none otherwise; (iii) if $|F|\ge2$, every $q\in F$ lies in some but not in all such $P$. So the first rule of PROPAGATE is: fail at an unsatisfied check without free qubits; put the only free qubit of a check into $T$ if the check is unsatisfied and into $X$ otherwise; nothing else.

*Proof.* The description of the $P$ is immediate. (i) and (ii) follow with $G\subseteq F$. (iii) For $q\ne q'$ in $F$, among $G=\{q\}$ and $G=\{q,q'\}$ one has the required parity and contains $q$, and among $G=\emptyset$ and $G=\{q'\}$ one has it and avoids $q$. $\square$

**Lemma N2 (degree-sequence bound).** Let $(T,X)$ be a state, $U$ the set of unsatisfied checks, $u=|U|$, and for a free qubit $q$ let $h(q)$ be the number of checks of $U$ containing $q$; let $h_1\ge h_2\ge\dots$ be these numbers over all free qubits. If $E$ is a set of free qubits meeting every check of $U$, then $|E|\ge t^*=\min\{t: h_1+\dots+h_t\ge u\}$, and no such $E$ exists if $\sum_q h(q)<u$. Moreover $t^*\ge\lceil u/h_1\rceil$.

*Proof.* $u\le\sum_{c\in U}|E\cap c|=\sum_{q\in E}h(q)\le h_1+\dots+h_{|E|}$. The last claim: $h_1+\dots+h_t\le t h_1$. $\square$

In the proof of ED, Proposition prop:exact, the only property of $\beta$ used is that $\beta(T,X)\le|E|$ for $E=z\setminus T$, which consists of free qubits and meets every unsatisfied check (ED, Lemma lem:minimal(c)); Lemma N2 and the packing bound of ED both have it. exactB reads $h$ from a histogram kept up to date in `flip`, so the bound costs $O(\max_q\deg q)$. (`cluster2.c`, option `-b 3`, computes this bound with $h$ capped at $3$; for qubits in more than three unsatisfied checks that overestimates $t^*$ and is not a valid lower bound. The published runs used `-b 2`, column weight $3$, so they are not affected; see "Limits of cluster2.c" below.)

**Lemma N3 (parity).** Every element of $\ker A$ has even weight iff $\mathbf 1\in\operatorname{rowsp}A$. In that case $d$ is even, and if there is no nontrivial logical of weight $\le W'=2\lfloor W/2\rfloor$, there is none of weight $\le W$ and $d\ge 2\lfloor W/2\rfloor+2$.

*Proof.* $\ker A=(\operatorname{rowsp}A)^\perp$, and $z\perp\mathbf 1$ iff $|z|$ is even; conversely, if $\ker A\subseteq\mathbf 1^\perp$, then $\mathbf 1\in(\ker A)^\perp=\operatorname{rowsp}A$. The rest follows. $\square$

ED, Lemma lem:parity is the case where every column of $A$ has odd weight (then the sum of all rows is $\mathbf 1$). exactB decides $\mathbf 1\in\operatorname{rowsp}A$ by reducing $\mathbf 1$ against the RREF of $A$ and then searches $W'$ (printed as `searched=`). Block parities (prompt, Section 8, C1, item 3) are not used.

**Lemma N4 (thread distribution).** Run the search with $t\ge2$ threads as follows: every thread walks the search tree depth first; at depth $D$ (option `--split`) the nodes are numbered in walk order; a thread executes the subtree of node $j$ iff $j$ equals its current claim, where claims are the values of a shared counter incremented atomically, and the thread takes a new claim after executing its claimed subtree; it skips all other depth-$D$ subtrees. Assume that every decision of GROW and PROPAGATE (choice of the check, order of the candidates, outcome and result of a propagation, value of $\beta$, leaf test) depends only on the state $(T,X)$ and the operation applied to it. Then every depth-$D$ subtree of the single-thread search is executed by exactly one thread.

*Proof.* Executing or skipping a subtree restores the state $(T,X)$ exactly (`undo`), so by the assumption all threads walk the same nodes above depth $D$ in the same order and number the $N$ depth-$D$ nodes $0,\dots,N-1$ identically. Claims are pairwise distinct values of the counter, so no node is executed twice. Each thread holds one claim at a time; its first claim is taken before its walk starts, and each later claim is taken right after executing node $j$ and exceeds $j$ (the counter passed $j$ when $j$ was claimed), so the walk still reaches the claimed node, and a claim $<N$ is executed. If some $j<N$ were never claimed, the counter would stay $\le j$, so every thread's last claim would be $<N$ and executed, after which the thread claims again, a contradiction. $\square$

The assumption holds in exactB: the check is the minimum of (free count, index) over the unsatisfied checks; candidates are taken in the stored order of the check; the propagation stack is empty before every operation and is then filled in an order determined by the state and the operation; and the outcome of a propagation does not depend on the processing order: its rules are monotone (if a rule applies to a state $S$ with conclusion $q\mapsto v$, then every consistent extension of $S$ that has examined the same check or member of $\mathcal S$ assigns $v$ to $q$, and failure conditions persist under extension), so every fair processing order either fails or ends in the same least fixpoint above the state, and $\beta$ is then evaluated on that fixpoint. $\beta$ reads the histogram of $h$ and the unsatisfied checks in index order (bit set `ubits`); the leaf test is the XOR over $T$. The order of the list `ulist`, which depends on the history, is never used for a decision. With threads, nodes above depth $D$ are counted by thread 0 only, so the node count equals that of the single-thread search unless a word stops the search early (checked below).

**Lemma N5 (generators).** Let $\gamma$ be a permutation of the qubits with $\gamma(r)\in\operatorname{rowsp}A$ for every row $r$ of $A$ and $\gamma(r)\in\operatorname{rowsp}B$ for every row $r$ of $B$. Then $\gamma(\operatorname{rowsp}A)=\operatorname{rowsp}A$ and $\gamma(\operatorname{rowsp}B)=\operatorname{rowsp}B$, and the group generated by such permutations satisfies the hypothesis of ED, Lemma lem:anchors.

*Proof.* $\gamma$ acts on $\mathbb F_2^n$ as an invertible linear map, so $\gamma(\operatorname{rowsp}A)\subseteq\operatorname{rowsp}A$ is an equality by dimension; the set of permutations preserving both subspaces is a group. $\square$

exactB checks the hypothesis for every generator (reduction of $\gamma(r)$ against the RREF of $A$, resp. $B$) and refuses the group otherwise. Its orbits (union-find over the generators) are ordered by least element (or by size, `--orbit-order size`), $q_i$ is the least element of $O_i$ and $X_i=O_1\cup\dots\cup O_{i-1}$; run $i$ sets $X_i$ to $X$, propagates, puts $q_i$ into $T$, propagates and calls GROW. This is ED, Algorithm alg:search with $T_0=\{q_i\}$, $X_0=X_i$, the propagation done in two steps (the proof of prop:exact only needs that each applied rule preserves compatibility). By ED, Corollary cor:runs, if no run reports a set, $d\ge W+1$. Without a group the orbits are the single qubits in the order `id` or `bfs`.

**Lemma N6 (logical basis).** Let $L_1,\dots,L_k\in\ker B$ with $\operatorname{rowsp}A+\langle L_1,\dots,L_k\rangle=\ker B$. For $z\in\ker A$: $z\in\operatorname{rowsp}B$ iff $\langle z,L_j\rangle=0$ for all $j$.

*Proof.* $\operatorname{rowsp}B=(\ker B)^\perp$, and $z\perp\operatorname{rowsp}A$ since $z\in\ker A$. $\square$

exactB checks $AB^{\mathsf T}=0$ entrywise, computes the RREF of $B$, takes the standard kernel basis (one vector per non-pivot column, each checked against all rows of $B$), and inserts first the RREF rows of $A$ (which lie in $\ker B$) and then the kernel vectors into an echelon basis; the kernel vectors that are not reduced to $0$ are the $L_j$, and their number must equal $n-\operatorname{rk}A-\operatorname{rk}B$ (checked). The class of $T$ is $\bigl(\langle 1_T,L_j\rangle\bigr)_j$, the XOR of the rows of the qubits of $T$ in an $n\times\lceil k/64\rceil$ word table; it is computed only at leaves ($\sigma(T)=0$). Reported words are rechecked without the $L_j$.

## What changes for general weights

* ED, Section 2 and 3 do not depend on the check weight, except Lemma lem:patterns (not used here) and Lemma lem:parity (generalized in N3).
* Implementation: the candidate array holds up to 64 free qubits (maximum check weight 64, checked at input); the histogram of $h$ has $\max_q\deg q+2$ entries; the propagation stacks are sized by the number of nonzeros of $A$, resp. by $\sum_{s\in\mathcal S}|s|$ (each qubit changes state at most once between two drains); classes use $\lceil k/64\rceil$ words; no part of the program assumes two blocks or a group structure.
* The default majority family includes all sums of two overlapping rows of $B$ (cap $2\max|B_r|$; `cluster2.c` used the cap $6\cdot$`majo`, i.e. all such sums at weight 6).

## Limits of cluster2.c found while generalizing (no effect on the published runs)

All concern parameters outside the published use (weight 6, column weight 3, options `-m 2 -b 2`, $k\le64$):

* `int cntH[8]` is indexed by the number of unsatisfied checks through a qubit: out of bounds for column weight $\ge8$.
* `-b 3` caps $h$ at $3$ (see Lemma N2): invalid bound for column weight $\ge4$; `-b 4` uses `int saved[8]` (column weight $\le8$).
* `cstack` has $8n+64$ entries: overflow possible for average column weight above $8$.
* $k\le64$ (one 64-bit class word), enforced at input.

## Tests on small codes (computed)

`tests/gen_tests.py tests/small 1 60` and `tests/gen_odd.py tests/odd 2 20` (seeds 1 and 2): 62 + 20 codes with $n\le64$ (two-block cyclic codes of types $w_a+w_b$, $1\le w_a,w_b\le6$, including odd types with parity; the same with up to $m$ redundant rows, sums of two or three checks, appended, so that the check weights are mixed; hypergraph products of random classical codes with mixed row and column weights; two direct sums with $n=750$, $k=68$ and $k=70$ whose lightest summand comes last, so that its logicals lie in the second class word). Distances by enumeration of $\ker A$ (dimension $\le22$) or of all supports of weight $\le5$, with a logical basis computed in Python. `tests/run_tests.sh` runs both sides, $W=d-1$ and $W=d$, anchors `id`, `bfs`, two threads, and the translation group where known (also on the instances with redundant rows, where the generator maps added rows to sums of checks, so that the rank test of Lemma N5 is exercised): 984 + 400 runs per build, 0 failures with the `-O3` build and with an AddressSanitizer/UBSan build (2026-11-10). A random permutation given as `--group` is refused (exit code 2).

## Validation on board codes (computed, 2026-11-10)

*Instances.* Board snapshot `8042f9da` (`gross_research/session3/qecchallenge/data/qldpc-challenge/codes/<slug>.json`, SHA-256 prefixes in `instances/board_sha256.txt`), relabelled by `scripts/relabel.py BOARD.json instances/<slug> 20261110` (seeded random qubit permutation, check order and order within lines; `instances/<slug>.perm`). The board witnesses, mapped to the new labels, lie in the expected kernels (`instances/<slug>.wit`). Certified values: `certs/<slug>.json` (tropidist exhaustive search; for `144-2-12` a CryptoMiniSat/DRAT certificate), $d_X=d_Z=d$ for all codes below.

*Groups.* `magma -b fhx:=... fhz:=... fout:=instances/<slug>.gens scripts/autgroup.m` (automorphisms of the Tanner graph with qubits, $X$-checks and $Z$-checks fixed as vertex classes); every generator was verified by exactB (Lemma N5). Orders on the qubits and orbits: `144-2-12` 144 (72, 72); `240-4-18` 120 (120, 120); `780-64-10` 12 (66 orbits: 6, 6 and 64 of size 12); `248-10-18` 124 (124, 124); `224-12-16` 56 (4 orbits of 56); `240-16-18` 60 (4 of 60); `224-18-16` 112 (112, 112); `248-24-18` 124 (124, 124); `232-62-12` 29 (8 of 29); `234-24-18` 117 (117, 117); `254-42-14` 127 (127, 127); `672-336-12` 84 (8 of 84); `252-26-18` 126 (126, 126); `240-22-18` 120 (120, 120).

*Parity* ($\mathbf 1\in\operatorname{rowsp}$, decided by exactB; a separate rank computation in Python agrees for every code and side that exactB ran): both sides of `144-2-12`, `248-10-18`, `224-12-16`, `232-62-12`, `254-42-14`, `672-336-12`, `240-22-18`; not for `240-4-18`, `780-64-10`, `240-16-18`, `224-18-16`, `248-24-18`, `252-26-18`, `234-24-18`. With parity, $W=d-1$ searches $d-2$.

*Protocol.* `scripts/batch1.sh` ($w=4,\dots,8$), `scripts/batch2.sh` ($w=9,10,12$): both sides, group anchors, $W=d-8,\dots,d$; `scripts/batch3.sh`: anchors at every qubit (no group, order `bfs`), side $Z$, $W=d-1,d$; two-thread runs; `scripts/batch4.sh`: the weight-8 codes `252-26-18`, `240-22-18` of T13, side $Z$. Default options (bound 3, majority family 2). Summary lines in `new_codes/logs/T-E-B/batch*.out`, full logs in `runs/`. `scripts/table.py` rebuilds the tables and rechecks every reported word against the board matrices in the original labelling (inverse permutation, syndrome and $\operatorname{rk}[H;z]>\operatorname{rk}H$ in Python): all words pass. Single-thread CPU seconds on the Ryzen 7 PRO 6850U, with the machine shared by about ten jobs.

**Result.** Every run gives the expected verdict: nothing at $W=d-1$ (searched $d-2$ with parity) and a nontrivial word of weight exactly $d$ at $W=d$, on both sides of all twelve codes of the brief (including $w=9,10,12$ and $k=336$), with group anchors, with anchors at every qubit (side $Z$, eight codes) and with two threads. Full table (also $W=d-8,\dots,d-2$): `runs/validation_table_full.md`.

| code | w | side | anchors | W (searched) | verdict | nodes | CPU s |
|---|---|---|---|---|---|---|---|
| 144-2-12 | 4 | X | group | 11 (10) | none | 536 | 0.00 |
| 144-2-12 | 4 | X | group | 12 (12) | word of weight 12 (board-recheck ok) | 24 | 0.00 |
| 144-2-12 | 4 | Z | group | 11 (10) | none | 536 | 0.00 |
| 144-2-12 | 4 | Z | group | 12 (12) | word of weight 12 (board-recheck ok) | 12 | 0.00 |
| 144-2-12 | 4 | Z | qubits | 11 (10) | none | 6414 | 0.01 |
| 144-2-12 | 4 | Z | qubits | 12 (12) | word of weight 12 (board-recheck ok) | 12 | 0.00 |
| 240-4-18 | 5 | X | group | 17 (17) | none | 79604 | 0.17 |
| 240-4-18 | 5 | X | group | 18 (18) | word of weight 18 (board-recheck ok) | 1514 | 0.00 |
| 240-4-18 | 5 | Z | group | 17 (17) | none | 99432 | 0.20 |
| 240-4-18 | 5 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 1705 | 0.00 |
| 240-4-18 | 5 | Z | qubits | 17 (17) | none | 1177896 | 2.12 |
| 240-4-18 | 5 | Z | qubits | 18 (18) | word of weight 18 (board-recheck ok) | 1705 | 0.00 |
| 780-64-10 | 5 | X | group | 9 (9) | none | 3053 | 0.01 |
| 780-64-10 | 5 | X | group | 10 (10) | word of weight 10 (board-recheck ok) | 207 | 0.00 |
| 780-64-10 | 5 | Z | group | 9 (9) | none | 3230 | 0.01 |
| 780-64-10 | 5 | Z | group | 10 (10) | word of weight 10 (board-recheck ok) | 45 | 0.00 |
| 224-12-16 | 6 | X | group | 15 (14) | none | 35160 | 0.11 |
| 224-12-16 | 6 | X | group | 16 (16) | word of weight 16 (board-recheck ok) | 10353 | 0.04 |
| 224-12-16 | 6 | Z | group | 15 (14) | none | 32349 | 0.11 |
| 224-12-16 | 6 | Z | group | 16 (16) | word of weight 16 (board-recheck ok) | 4300 | 0.02 |
| 248-10-18 | 6 | X | group | 17 (16) | none | 322334 | 1.12 |
| 248-10-18 | 6 | X | group | 18 (18) | word of weight 18 (board-recheck ok) | 12836 | 0.05 |
| 248-10-18 | 6 | Z | group | 17 (16) | none | 291456 | 0.90 |
| 248-10-18 | 6 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 29425 | 0.10 |
| 248-10-18 | 6 | Z | qubits | 17 (16) | none | 3550408 | 12.18 |
| 248-10-18 | 6 | Z | qubits | 18 (18) | word of weight 18 (board-recheck ok) | 29425 | 0.17 |
| 224-18-16 | 7 | X | group | 15 (15) | none | 154871 | 0.90 |
| 224-18-16 | 7 | X | group | 16 (16) | word of weight 16 (board-recheck ok) | 13460 | 0.07 |
| 224-18-16 | 7 | Z | group | 15 (15) | none | 121885 | 0.64 |
| 224-18-16 | 7 | Z | group | 16 (16) | word of weight 16 (board-recheck ok) | 7604 | 0.04 |
| 224-18-16 | 7 | Z | qubits | 15 (15) | none | 1988095 | 9.09 |
| 224-18-16 | 7 | Z | qubits | 16 (16) | word of weight 16 (board-recheck ok) | 7604 | 0.04 |
| 240-16-18 | 7 | X | group | 17 (17) | none | 1499981 | 8.66 |
| 240-16-18 | 7 | X | group | 18 (18) | word of weight 18 (board-recheck ok) | 337022 | 1.90 |
| 240-16-18 | 7 | X | threads2_split2 | 17 (17) | none | 1499981 | 10.98 |
| 240-16-18 | 7 | X | threads2_split4 | 17 (17) | none | 1499981 | 10.76 |
| 240-16-18 | 7 | Z | group | 17 (17) | none | 1650996 | 8.04 |
| 240-16-18 | 7 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 34961 | 0.22 |
| 232-62-12 | 8 | X | group | 11 (10) | none | 10703 | 0.06 |
| 232-62-12 | 8 | X | group | 12 (12) | word of weight 12 (board-recheck ok) | 528 | 0.00 |
| 232-62-12 | 8 | Z | group | 11 (10) | none | 10549 | 0.06 |
| 232-62-12 | 8 | Z | group | 12 (12) | word of weight 12 (board-recheck ok) | 685 | 0.00 |
| 232-62-12 | 8 | Z | qubits | 11 (10) | none | 151069 | 0.69 |
| 232-62-12 | 8 | Z | qubits | 12 (12) | word of weight 12 (board-recheck ok) | 685 | 0.00 |
| 240-22-18 | 8 | Z | group | 17 (16) | none | 1076696 | 10.05 |
| 240-22-18 | 8 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 1687193 | 14.12 |
| 248-24-18 | 8 | X | group | 17 (17) | none | 3800632 | 30.81 |
| 248-24-18 | 8 | X | group | 18 (18) | word of weight 18 (board-recheck ok) | 216689 | 1.91 |
| 248-24-18 | 8 | Z | group | 17 (17) | none | 3730228 | 27.44 |
| 248-24-18 | 8 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 1623371 | 12.92 |
| 248-24-18 | 8 | Z | qubits | 17 (17) | none | 44642759 | 309.69 |
| 248-24-18 | 8 | Z | qubits | 18 (18) | word of weight 18 (board-recheck ok) | 1623371 | 13.70 |
| 252-26-18 | 8 | Z | group | 17 (17) | none | 3379662 | 27.77 |
| 252-26-18 | 8 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 48827 | 0.41 |
| 234-24-18 | 9 | X | group | 17 (17) | none | 11140137 | 91.07 |
| 234-24-18 | 9 | X | group | 18 (18) | word of weight 18 (board-recheck ok) | 1232554 | 13.14 |
| 234-24-18 | 9 | Z | group | 17 (17) | none | 13067792 | 124.02 |
| 234-24-18 | 9 | Z | group | 18 (18) | word of weight 18 (board-recheck ok) | 359609 | 3.26 |
| 254-42-14 | 10 | X | group | 13 (12) | none | 25563 | 0.28 |
| 254-42-14 | 10 | X | group | 14 (14) | word of weight 14 (board-recheck ok) | 518 | 0.01 |
| 254-42-14 | 10 | Z | group | 13 (12) | none | 23756 | 0.27 |
| 254-42-14 | 10 | Z | group | 14 (14) | word of weight 14 (board-recheck ok) | 3497 | 0.04 |
| 254-42-14 | 10 | Z | qubits | 13 (12) | none | 348090 | 3.88 |
| 254-42-14 | 10 | Z | qubits | 14 (14) | word of weight 14 (board-recheck ok) | 3497 | 0.05 |
| 672-336-12 | 12 | X | group | 11 (10) | none | 55011 | 0.55 |
| 672-336-12 | 12 | X | group | 12 (12) | word of weight 12 (board-recheck ok) | 13124 | 0.13 |
| 672-336-12 | 12 | Z | group | 11 (10) | none | 50602 | 0.46 |
| 672-336-12 | 12 | Z | group | 12 (12) | word of weight 12 (board-recheck ok) | 4343 | 0.04 |
| 672-336-12 | 12 | Z | qubits | 11 (10) | none | 2165642 | 18.68 |
| 672-336-12 | 12 | Z | qubits | 12 (12) | word of weight 12 (board-recheck ok) | 4343 | 0.05 |

*Two threads* (`--threads 2`, split depth 2 and 4) reproduce the single-thread node and generated counts exactly (`248-24-18` $Z$, $W=16$: 1012335 nodes; `240-16-18` $X$, $W=17$: 1499981 nodes), as Lemma N4 predicts; wall time 4.8 s against 8.0 s, resp. 5.8 s against 8.9 s (CPU 9.2 to 11.0 s against 7.8, resp. 8.7 s: the part above the split depth is walked twice).

*Anchors at every qubit* cost 12 to 43 times the group runs (`248-24-18`, $W=17$: $4.46\cdot10^7$ against $3.73\cdot10^6$ nodes, 310 s against 27 s); order `bfs` needs about 0.65 times the nodes of order `id` (`248-10-18`, $W=12,14$).

## Growth of the node count (computed)

Least-squares fit of $\log N(W)=\log N_0+W\log\beta$ over the distinct searched values $W\in[d-8,d-1]$ (`scripts/growth.py`; with parity only even $W$ occur, the fit is still per unit of $W$). $N$ = calls of GROW, group anchors.

| $w$ | code | side | searched $W$ | $\beta$ |
|---|---|---|---|---|
| 4 | 144-2-12 | X, Z | 4..10 (even) | 2.25, 2.25 |
| 5 | 240-4-18 | X, Z | 10..17 | 2.39, 2.31 |
| 5 | 780-64-10 | X, Z | 3..9 | 2.27, 2.31 |
| 6 | 224-12-16 | X, Z | 8..14 (even) | 2.65, 2.73 |
| 6 | 248-10-18 | X, Z | 10..16 (even) | 2.79, 2.80 |
| 7 | 224-18-16 | X, Z | 8..15 | 3.10, 3.21 |
| 7 | 240-16-18 | X, Z | 10..17 | 3.08, 3.03 |
| 8 | 232-62-12 | X, Z | 4..10 (even) | 3.51, 3.48 |
| 8 | 248-24-18 | X, Z | 10..17 | 3.52, 3.53 |
| 8 | 252-26-18 | Z | 10..17 | 3.45 |
| 8 | 240-22-18 | Z | 10..16 (even) | 3.46 |
| 9 | 234-24-18 | X, Z | 10..17 | 3.90, 3.86 |
| 10 | 254-42-14 | X, Z | 6..12 (even) | 3.44, 3.43 |
| 12 | 672-336-12 | X, Z | 4..10 (even) | 4.34, 4.27 |

The local ratio $N(W)/N(W-1)$ grows slowly with $W$ (e.g. `248-24-18` $Z$: 3.52, 3.33, 3.52, 3.51, 3.59, 3.64, 3.68 for $W=11,\dots,17$). T13 (prompt, Section 7) asks for $\beta(8)$ over $W=13,\dots,17$ on `248-24-18`, `252-26-18`, `240-22-18`: 3.61 ($Z$ and $X$), 3.54, 3.55 (parity code, searched 12, 14, 16); so $\beta(8)\ge3$ holds for this implementation and T13 is not refuted. Time per node: about 7 to 10 $\mu$s at $w=8,9$ (4 generated children per node), so $W=17$ costs about 30 s at $w=8$ and 2 min at $w=9$ for $n\approx240$; extrapolated with $\beta=3.6$, $W=23$ at $w=8$ would cost about $28\cdot3.6^6\approx6\cdot10^4$ s (heuristic).

## Options (computed, side $Z$, group anchors)

| code, $W$ | bound 1 | bound 2 | bound 3 | majority 0 | majority 1 | majority 2 (default) |
|---|---|---|---|---|---|---|
| 248-10-18, 14 | 246810 (0.64 s) | 35716 (0.13 s) | 35716 (0.13 s) | 77523 (0.26 s) | 39392 (0.13 s) | 35716 (0.13 s) |
| 248-24-18, 15 | 3515130 (21.9 s) | 278080 (2.61 s) | 278080 (2.49 s) | 305666 (2.19 s) | 278328 (1.90 s) | 278080 (2.49 s) |
| 672-336-12, 10 | 347564 (2.68 s) | 50602 (0.59 s) | 50602 (0.57 s) | 50610 (0.36 s) | 50602 (0.45 s) | 50602 (0.57 s) |

Nodes (CPU seconds). The packing bound is what prunes (factor 7 to 13); the degree-sequence bound adds nothing on top of it here. Sums of two checks in $\mathcal S$ save few nodes; at $w=8$ `--majority 1` saves about 24 per cent of the time, at $w=12$ `--majority 0` about 37 per cent. Defaults were kept for the validation.

## Limits

* Check weight $\le64$ (`MAXW`); validated on board codes up to $w=12$ and on test codes with up to eight distinct check weights in one matrix and weights up to 18 (`tests/small/tbmix*`, `tests/small/hgp*`). Column weight and $k$ unbounded.
* Memory and setup: dense GF(2) matrices with $m_A+m_B+n$ rows of $n$ bits (RREF, kernel, logical basis), the CSS check costs $O(m_A\cdot\operatorname{nnz}B)$; fine for $n\lesssim3\cdot10^4$. Anchors at every qubit add an $O(n^2\max\deg)$ setup.
* Search cost grows like $\beta(w)^W$ with the factors above; the first orbit run dominates (later runs exclude earlier orbits and are cheap), so the cost is about that of one unrestricted anchor; anchors at every qubit cost 12 to 43 times more. Only the group of the instance file is used; the stabilizer of an anchor is not used for further symmetry breaking, and the pattern anchoring of ED, Lemma lem:patterns is not implemented.
* Threads: dynamic distribution of the subtrees at depth `--split` (default 4); the part above that depth is walked by every thread.
* Not handled: non-CSS codes, block parities (prompt, Section 8, C1, item 3), checkpointing and resumption of long runs (cluster2's `-C`), pair mode (cluster2's `-A`).

## CPU used

About 25 CPU-minutes in total (923 s in the logged runs of `runs/`, the rest for the test suites, including two passes with the sanitizer build, and probes); at most two cores at a time.
