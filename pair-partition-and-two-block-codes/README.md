# Pair-partition codes of check weight six and a two-block code over C28 x Dic3

Six CSS codes submitted to the [qLDPC challenge board](https://github.com/unitaryfoundation/qldpc-challenge), one
pull request per code: five quasi-cyclic codes lifted from sparse $6\times14$ pair-partition templates with checks
of weight 6 (Part A), and a weight-6 two-block group algebra code over $C_{28}\times\mathrm{Dic}_3$ (Part B). For
every code this directory holds its definition, the check matrices, two witnesses, the records of two exact
distance searches on both sides, the runs of two third-party exact programs, and a script that rebuilds and checks
all of it.

| board slug | code | construction | $w$ | qubit degrees X / Z | $k/n$ | $kd^2/n$ | distance |
|---|---|---|---|---|---|---|---|
| `294-42-11` | $[[294,42,11]]$ | pair-partition template, $P=21$ | 6 | 2 or 3 / 2 or 3 | 1/7 | 17.29 | exact |
| `322-46-12` | $[[322,46,12]]$ | pair-partition template, $P=23$ | 6 | 2 or 3 / 2 or 3 | 1/7 | 20.57 | exact |
| `434-62-13` | $[[434,62,13]]$ | the template of `322-46-12`, $P=31$ | 6 | 2 or 3 / 2 or 3 | 1/7 | 24.14 | exact |
| `518-74-14` | $[[518,74,14]]$ | pair-partition template, $P=37$ | 6 | 2 or 3 / 2 or 3 | 1/7 | 28.00 | exact |
| `742-106-15` | $[[742,106,15]]$ | pair-partition template, $P=53$ | 6 | 2 or 3 / 2 or 3 | 1/7 | 32.14 | exact |
| `672-28-24` | $[[672,28,24]]$ | two-block group algebra code over $C_{12}\rtimes_5C_{28}\cong C_{28}\times\mathrm{Dic}_3$ | 6 | 3 / 3 | 1/24 | 24.00 | exact |

$w$ is the largest check weight. "exact" means $d_X=d_Z=d$ by Theorem 5 below: two exhaustive search programs of
ours prove the lower bound on both sides, the witnesses give the upper bound. The distances are computed, not proved
by hand. Two third-party programs, dist-m4ri and tropidist (tropidist was not run on `672-28-24`), give the same
values without any symmetry reduction; Theorem 5 does not use them.

## Reproduce

Requirements: Python 3.8 or later (standard library only), `make` and `sha256sum`; for `make rerun` also a C
compiler and the directory `two-block-and-amplified-codes/` of this repository, which contains the two search
programs.

```
make test     # the checksums of all files (SHA256SUMS), then src/rebuild.py (a few seconds);
              # our output: logs/rebuild.log
make rerun    # compiles exactA and exactB from ../two-block-and-amplified-codes/src after checking their sha256
              # against src/programs.sha256, reruns the 24 exact searches of the six records/ directories and
              # compares them with our logs (about 45 CPU-seconds); our output: logs/rerun.log
make related  # reruns 672-28-24/related/k_covers.py (Part B, Lemma 8) and compares with k_covers.out (seconds)
```

`python3 src/rebuild.py SLUG` checks one code. `python3 src/rebuild.py --entry FILE`, with `FILE` a board entry
`codes/<slug>.json`, checks that its check matrices equal `hx.txt` and `hz.txt` of `<slug>/` row by row, its
parameters, and that its witnesses are nontrivial logicals of the claimed weights. `make rerun PROGS=DIR` takes
`exactA.c` and `exactB.c` from `DIR` (they must have the checksums of `src/programs.sha256`). `src/rerun.py`
compares a rerun with our log through `src/same_results.py` (identical to the file of that name in
`two-block-and-amplified-codes/src/`): the `RESULT` line without its timing fields, every `WORD` line and the
`CONCLUSION` line must agree; both programs are deterministic for a fixed instance, generator file and number of
threads. Our runs used gcc with `-O3 -march=native` and one thread per run.

## Part A: pair-partition codes from a sparse $6\times14$ template

### Construction

**Definition 1.** Let $J=6$, $L=14$ and $P\ge2$. A *template* is a pair of $J\times L$ binary arrays $T_X$, $T_Z$,
given by their row supports $X_0,\dots,X_5$ and $Z_0,\dots,Z_5\subseteq\{0,\dots,13\}$, such that

1. every row has weight 6;
2. every column has degree 2 or 3 in $T_X$ and degree 2 or 3 in $T_Z$, and exactly six columns have degree 2 in
   each array;
3. the six columns of degree 2 of $T_X$, each read as an edge between the two rows of $T_X$ containing it, form a
   connected unicyclic graph on the six rows, and the same holds for $T_Z$;
4. every overlap $X_i\cap Z_j$ has even size and carries a partition $\Pi_{ij}$ into pairs, and the pairing array
   $(\Pi_{ij})$ satisfies the girth rule of Okada and Kasai in the form used by our search: no pair $\{u,v\}$
   occurs in two cells of one row $i$, or of one column $j$, of the array.

Exponents $e_{iu}\in\mathbb Z/P$ ($u\in X_i$) and $d_{ju}\in\mathbb Z/P$ ($u\in Z_j$) are required to satisfy the
*paired-difference equations* $e_{iu}-d_{ju}=e_{iv}-d_{jv}$ for every pair $\{u,v\}\in\Pi_{ij}$. The code has the
$n=LP$ qubits $uP+s$ ($u\in\{0,\dots,13\}$ a block column, $s\in\mathbb Z/P$). The X-check $(i,r)$, $r\in\mathbb Z/P$,
is row $iP+r$ of $H_X$ and acts on the qubits $uP+((r+e_{iu})\bmod P)$, $u\in X_i$; the Z-check $(j,r)$ is row
$jP+r$ of $H_Z$ and acts on the qubits $uP+((r+d_{ju})\bmod P)$, $u\in Z_j$. Thus $H_X$ is the $6\times14$ array of
$P\times P$ blocks with the circulant permutation matrix $S^{e_{iu}}$ at the positions of $T_X$ and zero elsewhere
($S$ the cyclic shift), and likewise $H_Z$. Every check has weight 6, a qubit of block column $u$ lies in
$\deg_{T_X}(u)\in\{2,3\}$ X-checks and in $\deg_{T_Z}(u)\in\{2,3\}$ Z-checks, and the design rate is
$(L-2J)/L=1/7$.

**Lemma 2 (orthogonality).** In the situation of Definition 1, $H_XH_Z^{\mathsf T}=0$.

*Proof.* The X-check $(i,r)$ and the Z-check $(j,s)$ meet in the qubits $uP+t$ with $u\in X_i\cap Z_j$ and
$t\equiv r+e_{iu}\equiv s+d_{ju}$, that is, in one qubit of each block column $u\in X_i\cap Z_j$ with
$e_{iu}-d_{ju}\equiv s-r$, and in no other. By the paired-difference equations this condition holds for $u$ iff it
holds for the partner of $u$ in $\Pi_{ij}$, so the overlap consists of whole pairs and has even size. $\square$

This is the orthogonality argument of Okada and Kasai (arXiv:2607.14091) for their all-one protographs; it does not
use that the template is all-one. The exponents modulo the relabellings $s\mapsto s+c_u$ of the block columns and
$r\mapsto r+c_i$, $r\mapsto r+c_j$ of the check rows are the 1-cocycles modulo coboundaries of the *pairing square
complex* (a vertex for every row of $T_X$, every row of $T_Z$ and every column, an edge for every incidence, a
square $(i,u,j,v)$ for every pair $\{u,v\}\in\Pi_{ij}$): the paired-difference equation of a pair is the cocycle
condition on its square. Our search therefore organized the exponents of a template by $H^1$ of this complex with
coefficients in $\mathbb Z/P$ (job labels in `definition.json`).

**Lemma 3.** $k\ge(L-2J)P=2P$, with equality iff $H_X$ and $H_Z$ have full row rank $6P$.

*Proof.* $H_X$ and $H_Z$ have $JP$ rows each, and $k=n-\operatorname{rank}H_X-\operatorname{rank}H_Z$. $\square$

### The five codes

Template rows and exponents, as `column:exponent` (the same data are in `<slug>/definition.json`, fields
`template`, `exponents` and `pairing`; `src/rebuild.py` builds $H_X$, $H_Z$ from them and compares with `hx.txt`,
`hz.txt` row by row). `322-46-12` and `434-62-13` have the same template.

`294-42-11`, $P=21$:

| row | columns : exponents | row | columns : exponents |
|---|---|---|---|
| $X_0$ | 0:0 1:0 2:0 4:0 5:0 12:0 | $Z_0$ | 0:0 1:0 2:14 3:15 5:14 6:0 |
| $X_1$ | 0:0 2:20 3:0 5:14 7:0 13:0 | $Z_1$ | 0:0 2:0 3:0 7:1 8:0 11:0 |
| $X_2$ | 1:0 2:14 4:19 8:14 9:19 10:0 | $Z_2$ | 0:0 3:11 4:0 7:11 9:0 13:0 |
| $X_3$ | 3:10 4:0 6:16 7:11 9:20 10:12 | $Z_3$ | 1:0 4:19 5:19 6:14 12:0 13:5 |
| $X_4$ | 5:0 6:7 7:9 11:8 12:2 13:19 | $Z_4$ | 4:0 6:5 8:16 9:20 10:1 12:0 |
| $X_5$ | 8:1 9:5 10:1 11:1 12:0 13:5 | $Z_5$ | 5:0 7:7 8:1 10:8 11:8 12:0 |

`322-46-12`, $P=23$:

| row | columns : exponents | row | columns : exponents |
|---|---|---|---|
| $X_0$ | 0:0 1:0 2:0 3:0 5:0 10:0 | $Z_0$ | 0:0 1:0 4:0 6:0 12:0 13:0 |
| $X_1$ | 0:0 2:22 4:0 5:3 9:0 11:0 | $Z_1$ | 0:0 2:0 3:4 5:4 9:0 13:16 |
| $X_2$ | 1:0 2:9 6:0 7:0 8:16 9:9 | $Z_2$ | 0:0 6:7 8:0 10:0 11:0 13:13 |
| $X_3$ | 3:0 6:12 7:19 10:5 11:22 13:12 | $Z_3$ | 1:0 3:0 4:5 7:19 8:16 9:5 |
| $X_4$ | 4:12 5:0 6:6 7:3 12:6 13:12 | $Z_4$ | 2:0 4:20 5:0 6:14 11:1 12:6 |
| $X_5$ | 8:8 9:20 10:0 11:8 12:13 13:13 | $Z_5$ | 2:0 7:14 8:8 9:1 10:0 12:17 |

`434-62-13`, $P=31$:

| row | columns : exponents | row | columns : exponents |
|---|---|---|---|
| $X_0$ | 0:0 1:0 2:0 3:0 5:0 10:0 | $Z_0$ | 0:0 1:0 4:0 6:0 12:0 13:0 |
| $X_1$ | 0:0 2:30 4:0 5:24 9:0 11:0 | $Z_1$ | 0:0 2:0 3:25 5:25 9:0 13:10 |
| $X_2$ | 1:0 2:13 6:0 7:0 8:3 9:13 | $Z_2$ | 0:0 6:28 8:0 10:0 11:0 13:17 |
| $X_3$ | 3:0 6:16 7:6 10:19 11:30 13:16 | $Z_3$ | 1:0 3:0 4:19 7:6 8:3 9:19 |
| $X_4$ | 4:16 5:0 6:27 7:3 12:27 13:16 | $Z_4$ | 2:0 4:7 5:0 6:18 11:1 12:27 |
| $X_5$ | 8:22 9:7 10:0 11:22 12:17 13:17 | $Z_5$ | 2:0 7:18 8:22 9:1 10:0 12:11 |

`518-74-14`, $P=37$:

| row | columns : exponents | row | columns : exponents |
|---|---|---|---|
| $X_0$ | 0:0 1:0 2:0 5:0 9:0 10:0 | $Z_0$ | 0:0 1:24 2:24 7:0 9:0 13:0 |
| $X_1$ | 0:0 2:24 3:0 4:0 5:22 12:0 | $Z_1$ | 0:0 1:0 5:22 6:0 8:0 9:22 |
| $X_2$ | 1:0 2:14 4:20 7:27 8:0 13:13 | $Z_2$ | 1:0 3:20 4:20 6:8 10:0 11:0 |
| $X_3$ | 3:27 4:29 6:15 7:0 9:0 10:9 | $Z_3$ | 2:0 3:13 4:15 5:0 8:23 13:8 |
| $X_4$ | 5:0 6:15 7:8 11:7 12:7 13:8 | $Z_4$ | 2:0 6:6 10:0 11:13 12:13 13:36 |
| $X_5$ | 8:15 9:0 10:14 11:14 12:27 13:0 | $Z_5$ | 3:27 7:28 8:1 9:0 10:0 12:27 |

`742-106-15`, $P=53$:

| row | columns : exponents | row | columns : exponents |
|---|---|---|---|
| $X_0$ | 0:0 1:0 2:0 6:0 10:0 12:0 | $Z_0$ | 0:0 1:0 3:0 5:0 9:0 11:0 |
| $X_1$ | 0:0 2:20 3:0 4:0 8:0 13:0 | $Z_1$ | 0:0 4:0 8:50 9:27 12:0 13:50 |
| $X_2$ | 1:0 2:38 4:32 5:0 7:0 8:29 | $Z_2$ | 0:0 6:0 10:14 11:15 12:14 13:0 |
| $X_3$ | 3:6 4:32 6:0 7:15 9:6 10:14 | $Z_3$ | 1:0 2:0 3:44 7:0 8:44 13:33 |
| $X_4$ | 5:4 6:0 7:20 11:4 12:3 13:0 | $Z_4$ | 1:0 2:38 4:18 5:39 10:0 12:38 |
| $X_5$ | 8:35 9:1 10:0 11:1 12:38 13:24 | $Z_5$ | 2:0 5:4 6:0 7:15 8:33 11:52 |

Computed by `src/rebuild.py`, for all five codes: the template satisfies items 1 to 4 of Definition 1 and the
exponents the paired-difference equations; $\operatorname{rank}H_X=\operatorname{rank}H_Z=6P$, so $k=2P=n/7$
(Lemma 3); every check has weight 6; $6P$ qubits have X-degree 2 and $8P$ have X-degree 3, and the same for Z; no two
rows of $H_X$, and no two rows of $H_Z$, share two qubits (no 4-cycle within $H_X$ or within $H_Z$); the Tanner graph
is connected; $\mathbf 1$ lies in $\operatorname{rowsp}H_X$ and in $\operatorname{rowsp}H_Z$ for `518-74-14` and in
neither for the other four codes.

### How the codes were found

By a computer search of ours (job `l3_templates`, 2026-11-10; computed, labels in `definition.json`, field `job`):
all template pairs satisfying items 1 to 3 of Definition 1, up to permutations of rows and columns and the exchange
of X and Z, with the girth rule as a pruning condition (19756 single-type templates, 4755087 labelled pairs that
admit pairing arrays); 12954 classes of template pairs admit a pairing array whose pairing square complex has
first Betti number $b_1\ge3$, and a seeded sample of 600 of these classes was processed (4056 array orbits); 3835
orbits were dropped because a cofactor word of weight below 16 is a logical in every lift whose $P$ is prime to the
torsion of $H_1$ of the pairing square complex (3208 of them have torsion-free $H_1$; the other lifts of the remaining
627 were not examined), and 221 were kept; for the kept arrays and every $P$ from 5 to 71 a sample of cohomology
classes (356432 lifts; the job's control lifts on arrays with $b_1\le2$ are not counted), upper bounds by cofactor
words and by an information-set search, and the exact certification of the best lifts. No lift screened by the
information-set search had an upper bound above 15.
The template of `742-106-15` has $b_1=4$, the other three templates $b_1=3$; $H_1$ of all four complexes is
torsion-free.

### Context on the board

At board commit `82a52ee51` (2026-11-11, 3194 entries), no entry dominates or equals any of the six codes of this
directory in $(n,k,d,w)$, every one of them is admissible (`742-106-15` has $n>700$, check weight $6\le8$ and
$d\le40$), and the largest claimed distance of a CSS entry with $n'\le n$, $k'\ge k$ and check weight at most 6 is 6
for the four smallest codes, 8 for `742-106-15` and 20 for `672-28-24` (computed from the board's entries with the
board's domination relation).

## The distances of all six codes

**Lemma 4 (parity).** If $AB^{\mathsf T}=0$ and $\mathbf 1\in\operatorname{rowsp}A$, every element of $\ker A$ has
even weight. (Proposition 4(2) of `board-certificates-2026-11/README.md`: $\mathbf 1$ is a sum of rows of $A$, each
orthogonal to $x\in\ker A$, and $\langle\mathbf 1,x\rangle\equiv|x|\pmod2$.)

**Theorem 5.** For each of the six codes of this directory (Part A above, Part B below), $d_X=d_Z=d$, where $d$ is
the third number of its slug (computed).

*Proof.* Upper bound: `witness_X.txt` (in $\ker H_Z$, outside $\operatorname{rowsp}H_X$) and `witness_Z.txt` (in
$\ker H_X$, outside $\operatorname{rowsp}H_Z$) have weight $d$ (checked by `src/rebuild.py`). Lower bound: exactA and
exactB are exhaustive searches. Given $H_X$, $H_Z$, a side, a bound $W$ and generators of a group of qubit
permutations that maps the rows of $H_X$ onto rows of $H_X$ and the rows of $H_Z$ onto rows of $H_Z$, a run that
reports no word proves that the side has no nontrivial logical of weight at most $W$; when $\mathbf 1$ lies in the row
space of the side's check matrix, the program searches only the even weights up to $W$, which suffices by Lemma 4
(each program decides this by a rank computation and prints the decision). For one anchor per qubit orbit (exactB,
and exactA on `672-28-24`) this is Proposition 4 of `board-certificates-2026-11/README.md` in this repository; for
exactA's pattern classes of check orbits (exactA on the five codes of Part A) it is Lemma A with Corollary B of
`two-block-and-amplified-codes/src/exactA_NOTES.md`. The completeness proofs of both programs are in
`exactA_NOTES.md` and `exactB_NOTES.md` there (Lemma A, Corollary B, Lemmas C to E and G; Lemmas N1 to N6), on top
of Proposition 3.1 and Corollary 3.2 of `weight6-exact-distances/exact_distances.pdf`. The programs that produced
the records below are byte-identical to the copies in `two-block-and-amplified-codes/src/` and
`board-certificates-2026-11/src/` (sha256 in `src/programs.sha256`). For every code both programs ran on both
sides at $W=d-1$ and reported no word (tables below, logs in `<slug>/records/`). Both programs verify that the
generators they are given preserve the rows, and `src/rebuild.py` checks this again and computes the generated
group. $\square$

The groups: for Part A the translations $uP+s\mapsto uP+(s+1)$ of order $P$, whose orbits are the 14 block columns
(for `742-106-15` the generator files also contain the generator of Magma's Tanner-graph group, which is this
translation group); for Part B the Tanner-graph group of order 168 (Part B below). Where the runs ran:
`742-106-15` and `672-28-24` on our compute server (jobs `l3_templates`, stage 2, and `f1`, stage 5, 2026-11-10; the
jobs' JSON records are beside the logs); the four smaller codes of Part A on a laptop (2026-11-11), because the job's
stage 2 scheduled only one code. In our verification of 2026-11-11 all 24 searches were rerun on byte-identical copies
of the instance files with the same programs and gave the same node counts; `make rerun` repeats this.

Columns: program and its decomposition, side, $W$ (in parentheses the largest weight searched after the parity
reduction), result, nodes of the search tree, CPU seconds, log.

`294-42-11`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, pattern classes | X | 10 (10) | none | 4,284 | 0.01 | `records/A_X_W10.log` |
| exactA, pattern classes | Z | 10 (10) | none | 4,069 | 0.01 | `records/A_Z_W10.log` |
| exactB, orbit anchors | X | 10 (10) | none | 4,721 | 0.00 | `records/B_X_W10.log` |
| exactB, orbit anchors | Z | 10 (10) | none | 4,873 | 0.01 | `records/B_Z_W10.log` |

`322-46-12`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, pattern classes | X | 11 (11) | none | 10,045 | 0.01 | `records/A_X_W11.log` |
| exactA, pattern classes | Z | 11 (11) | none | 13,305 | 0.01 | `records/A_Z_W11.log` |
| exactB, orbit anchors | X | 11 (11) | none | 11,722 | 0.01 | `records/B_X_W11.log` |
| exactB, orbit anchors | Z | 11 (11) | none | 15,797 | 0.02 | `records/B_Z_W11.log` |

`434-62-13`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, pattern classes | X | 12 (12) | none | 25,003 | 0.03 | `records/A_X_W12.log` |
| exactA, pattern classes | Z | 12 (12) | none | 32,523 | 0.03 | `records/A_Z_W12.log` |
| exactB, orbit anchors | X | 12 (12) | none | 28,822 | 0.03 | `records/B_X_W12.log` |
| exactB, orbit anchors | Z | 12 (12) | none | 40,114 | 0.04 | `records/B_Z_W12.log` |

`518-74-14`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, pattern classes | X | 13 (12) | none | 30,004 | 0.03 | `records/A_X_W13.log` |
| exactA, pattern classes | Z | 13 (12) | none | 25,174 | 0.03 | `records/A_Z_W13.log` |
| exactB, orbit anchors | X | 13 (12) | none | 34,536 | 0.04 | `records/B_X_W13.log` |
| exactB, orbit anchors | Z | 13 (12) | none | 31,399 | 0.03 | `records/B_Z_W13.log` |

`742-106-15`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, pattern classes | X | 14 (14) | none | 156,789 | 0.23 | `records/A_X_W14.log` |
| exactA, pattern classes | Z | 14 (14) | none | 140,548 | 0.21 | `records/A_Z_W14.log` |
| exactB, orbit anchors | X | 14 (14) | none | 175,740 | 0.23 | `records/B_X_W14.log` |
| exactB, orbit anchors | Z | 14 (14) | none | 185,616 | 0.28 | `records/B_Z_W14.log` |

`672-28-24`:

| program | side | $W$ (searched) | result | nodes | CPU s | log |
|---|---|---|---|---|---|---|
| exactA, orbit anchors | X | 23 (22) | none | 7,199,307 | 14.92 | `records/A_X_W23.log` |
| exactA, orbit anchors | Z | 23 (22) | none | 7,345,978 | 15.29 | `records/A_Z_W23.log` |
| exactB, orbit anchors | X | 23 (22) | none | 7,370,526 | 17.04 | `records/B_X_W23.log` |
| exactB, orbit anchors | Z | 23 (22) | none | 7,507,311 | 17.36 | `records/B_Z_W23.log` |

### Third-party exact programs (not used in Theorem 5)

dist-m4ri (L. P. Pryadko, W. Zeng; github.com/QEC-pages/dist-m4ri, version 0.16.1, commit `4cc9c126`; `method=2`, the
connected-cluster search of Dumer, Kovalev and Pryadko, IEEE Trans. Inf. Theory 63 (2017) 4675-4686,
arXiv:1611.07164) and tropidist (Cheng, Guo, Pan, arXiv:2611.08247; github.com/Fanerst/tropidist
at commit `b03ec2d4077f`, the commit of the board's tropidist certificates; `python -m cdist --matrices H L`), each
run single-threaded on the MatrixMarket files in `<slug>/third_tools/mtx/` (side Z: $H=H_X$, side X: $H=H_Z$;
dist-m4ri computes the logicals from the other check matrix or reads the basis $L$, tropidist reads $L$). Every run
without a symmetry searches from every qubit. On `672-28-24` dist-m4ri also ran with `start=0,168,336,504
start_scan=1`, one start column per orbit of its Tanner-graph group (Part B); without symmetry it first excluded all
weights up to 21 and then, starting from the proved bound 22, the weight 22 (odd weights are excluded by parity,
which dist-m4ri decides itself). tropidist was not run on `672-28-24` (our estimate of more than 30 CPU-minutes for a
symmetry-free run, not measured). dist-m4ri's result line is "proved lower bound, weight of the words found (0 if
none), 0"; its `.nz` files list codewords with 1-based columns, tropidist's records give 0-based witnesses.
`src/rebuild.py` checks the inputs (`HX.mtx`, `HZ.mtx` equal $H_X$, $H_Z$; `LX.mtx`, `LZ.mtx` are bases of the
logical classes), every result line and record against `third_tools/runs.json`, and every returned word (a nontrivial
logical of the right side and weight $d$).

`294-42-11`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=0, weights up to 11; finH=HZ.mtx finL=LZ.mtx | $d_X\ge11$, words of weight 11 found (399 listed) | `third_tools/dist-m4ri/dm0L_X_w11.out` |
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 11; finH=HZ.mtx finG=HX.mtx | $d_X\ge11$, words of weight 11 found (399 listed) | `third_tools/dist-m4ri/dm_X_w11.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=0, weights up to 9; finH=HX.mtx finG=HZ.mtx | $d_Z\ge10$ | `third_tools/dist-m4ri/dm0_Z_w9.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 9; finH=HX.mtx finG=HZ.mtx | $d_Z\ge10$ | `third_tools/dist-m4ri/dm_Z_w9.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=0, weights up to 11; finH=HX.mtx finL=LX.mtx | $d_Z\ge11$, words of weight 11 found (399 listed) | `third_tools/dist-m4ri/dm0L_Z_w11.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 11; finH=HX.mtx finG=HZ.mtx | $d_Z\ge11$, words of weight 11 found (399 listed) | `third_tools/dist-m4ri/dm_Z_w11.out` |
| tropidist | X | no symmetry; `--matrices` with HZ.mtx LZ.mtx | exact, $d_X=11$ (0.2 s) | `third_tools/tropidist/rec_X.json` |
| tropidist | Z | no symmetry; `--matrices` with HX.mtx LX.mtx | exact, $d_Z=11$ (1.5 s) | `third_tools/tropidist/rec_Z.json` |

`322-46-12`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 12; finH=HZ.mtx finG=HX.mtx | $d_X\ge12$, words of weight 12 found | `third_tools/dist-m4ri/dX_G.out` |
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 12; finH=HZ.mtx finL=LZ.mtx | $d_X\ge12$, words of weight 12 found | `third_tools/dist-m4ri/dX_L.out` |
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 12; finH=HZ.mtx finG=HX.mtx | $d_X\ge12$, words of weight 12 found (322 listed) | `third_tools/dist-m4ri/dX_enum.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 9; finH=HX.mtx finG=HZ.mtx | $d_Z\ge10$ | `third_tools/dist-m4ri/dZ_wmax9.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge12$, words of weight 12 found | `third_tools/dist-m4ri/dZ_G.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finL=LX.mtx | $d_Z\ge12$, words of weight 12 found | `third_tools/dist-m4ri/dZ_L.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge12$, words of weight 12 found (322 listed) | `third_tools/dist-m4ri/dZ_enum.out` |
| tropidist | X | no symmetry; `--matrices` with HZ.mtx LZ.mtx | exact, $d_X=12$ (0.3 s) | `third_tools/tropidist/dX.json` |
| tropidist | Z | no symmetry; `--matrices` with HX.mtx LX.mtx | exact, $d_Z=12$ (1.5 s) | `third_tools/tropidist/dZ.json` |

`434-62-13`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 12; finH=HZ.mtx finG=HX.mtx | $d_X\ge13$ | `third_tools/dist-m4ri/m4ri_X_w12.log` |
| dist-m4ri | X | no symmetry; CC, cc_prune=0, weights up to 12; finH=HZ.mtx finG=HX.mtx | $d_X\ge13$ | `third_tools/dist-m4ri/m4ri_X_w12_prune0.log` |
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 13; finH=HZ.mtx finG=HX.mtx | $d_X\ge13$, words of weight 13 found (248 listed) | `third_tools/dist-m4ri/m4ri_X_w13.log` |
| dist-m4ri | Z | no symmetry; CC, weights up to 9; finH=HX.mtx finG=HZ.mtx | $d_Z\ge10$ | `third_tools/dist-m4ri/probe_Z_w9.log` |
| dist-m4ri | Z | no symmetry; CC, weights up to 10; finH=HX.mtx finG=HZ.mtx | $d_Z\ge11$ | `third_tools/dist-m4ri/probe_Z_w10.log` |
| dist-m4ri | Z | no symmetry; CC, weights up to 11; finH=HX.mtx finG=HZ.mtx | $d_Z\ge12$ | `third_tools/dist-m4ri/probe_Z_w11.log` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge13$ | `third_tools/dist-m4ri/m4ri_Z_w12.log` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=0, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge13$ | `third_tools/dist-m4ri/m4ri_Z_w12_prune0.log` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 13; finH=HX.mtx finG=HZ.mtx | $d_Z\ge13$, words of weight 13 found (248 listed) | `third_tools/dist-m4ri/m4ri_Z_w13.log` |
| tropidist | X | no symmetry; `--matrices` with HZ.mtx LZ.mtx | exact, $d_X=13$ (0.8 s) | `third_tools/tropidist/rec_X.json` |
| tropidist | Z | no symmetry; `--matrices` with HX.mtx LX.mtx | exact, $d_Z=13$ (2.1 s) | `third_tools/tropidist/rec_Z.json` |

`518-74-14`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 14; finH=HZ.mtx finG=HX.mtx | $d_X\ge14$, words of weight 14 found (3182 listed) | `third_tools/dist-m4ri/dm_X_w14.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 10; finH=HX.mtx finG=HZ.mtx | $d_Z\ge12$ | `third_tools/dist-m4ri/dm_Z_w10.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 11; finH=HX.mtx finG=HZ.mtx | $d_Z\ge12$ | `third_tools/dist-m4ri/dm_Z_w11.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge14$ | `third_tools/dist-m4ri/dm_Z_w12.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 14; finH=HX.mtx finG=HZ.mtx | $d_Z\ge14$, words of weight 14 found (3182 listed) | `third_tools/dist-m4ri/dm_Z_w14.out` |
| tropidist | X | no symmetry; `--matrices` with HZ.mtx LZ.mtx | exact, $d_X=14$ (3.0 s) | `third_tools/tropidist/trop_X.json` |
| tropidist | Z | no symmetry; `--matrices` with HX.mtx LX.mtx | exact, $d_Z=14$ (3.3 s) | `third_tools/tropidist/trop_Z.json` |

`742-106-15`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 15; finH=HZ.mtx finG=HX.mtx | $d_X\ge15$, words of weight 15 found (371 listed) | `third_tools/dist-m4ri/dm_X_w15.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 6; finH=HX.mtx finG=HZ.mtx | $d_Z\ge7$ | `third_tools/dist-m4ri/dm_Z_w6.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 7; finH=HX.mtx finG=HZ.mtx | $d_Z\ge8$ | `third_tools/dist-m4ri/dm_Z_w7.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 8; finH=HX.mtx finG=HZ.mtx | $d_Z\ge9$ | `third_tools/dist-m4ri/dm_Z_w8.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 12; finH=HX.mtx finG=HZ.mtx | $d_Z\ge13$ | `third_tools/dist-m4ri/dm_Z_w12.out` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 15; finH=HX.mtx finG=HZ.mtx | $d_Z\ge15$, words of weight 15 found (371 listed) | `third_tools/dist-m4ri/dm_Z_w15.out` |
| tropidist | X | no symmetry; `--matrices` with HZ.mtx LZ.mtx | exact, $d_X=15$ (9.6 s) | `third_tools/tropidist/trop_X.json` |
| tropidist | Z | no symmetry; `--matrices` with HX.mtx LX.mtx | exact, $d_Z=15$ (10.3 s) | `third_tools/tropidist/trop_Z.json` |

`672-28-24`:

| program | side | run | result | file |
|---|---|---|---|---|
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights up to 20; finH=HZ.mtx finG=HX.mtx | $d_X\ge22$ | `third_tools/dist-m4ri/full_X_w20.log` |
| dist-m4ri | X | no symmetry; CC, cc_prune=1, weights 22 to 22; finH=HZ.mtx finG=HX.mtx | $d_X\ge24$ (starting from $d_X\ge22$, proved by the run above) | `third_tools/dist-m4ri/full_X_w22.log` |
| dist-m4ri | X | start columns 0,168,336,504 (one per orbit); CC, cc_prune=1, weights up to 24; finH=HZ.mtx finG=HX.mtx | $d_X\ge24$, words of weight 24 found | `third_tools/dist-m4ri/sym_X_w24.log` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights up to 20; finH=HX.mtx finG=HZ.mtx | $d_Z\ge22$ | `third_tools/dist-m4ri/full_Z_w20.log` |
| dist-m4ri | Z | no symmetry; CC, cc_prune=1, weights 22 to 22; finH=HX.mtx finG=HZ.mtx | $d_Z\ge24$ (starting from $d_Z\ge22$, proved by the run above) | `third_tools/dist-m4ri/full_Z_w22.log` |
| dist-m4ri | Z | start columns 0,168,336,504 (one per orbit); CC, cc_prune=1, weights up to 24; finH=HX.mtx finG=HZ.mtx | $d_Z\ge24$, words of weight 24 found | `third_tools/dist-m4ri/sym_Z_w24.log` |

## Part B: the $[[672,28,24]]$ code

### Construction

$G=C_{12}\rtimes_5C_{28}=\langle x,y\mid x^{12}=y^{28}=1,\ yxy^{-1}=x^5\rangle$. Every element is $x^iy^j$ with
$0\le i<12$, $0\le j<28$, and $x^iy^j\cdot x^ky^l=x^{i+5^jk}y^{j+l}$ (well defined because $5^2\equiv1\pmod{12}$
and $28$ is even). The element $x^iy^j$ has index $(j\bmod2)\cdot168+14i+\lfloor j/2\rfloor$ (the numbering of our
search program), and qubit $(e,g)$ of block $e\in\{0,1\}$ has index $336e+\mathrm{index}(g)$. With
$$a=1+x^2y^9+x^{10}y^{11},\qquad b=1+x^{10}y^{22}+x^7y^{27}$$
(indices 0, 200, 313 and 0, 151, 279), the code is $C(G;a,b)$: $H_X=(L_a\mid R_b)$,
$H_Z=(R_b^{\mathsf T}\mid L_a^{\mathsf T})$ with $(L_a)_{c,g}=1$ iff $cg^{-1}\in a$ and $(R_b)_{c,h}=1$ iff
$h^{-1}c\in b$. Explicitly, the X-check at $c\in G$ (row $\mathrm{index}(c)$) consists of the block-0 qubits
$s^{-1}c$ ($s\in a$) and the block-1 qubits $ct^{-1}$ ($t\in b$), the Z-check at $c$ of the block-0 qubits $ct$
($t\in b$) and the block-1 qubits $sc$ ($s\in a$). Left and right multiplications commute, so
$H_XH_Z^{\mathsf T}=0$. These are the two-block group algebra codes of Lin and Pryadko (arXiv:2306.16400), in the
convention of `two-block-and-amplified-codes/` (Part A there). The group, in this presentation, is that of the board
entry `672-20-32` (a two-block code of vprusso over the same group, which numbers $x^iy^j$ as $28i+j$).

**Lemma 6.** $G\cong C_{28}\times\mathrm{Dic}_3\cong\mathbb Z_{84}\rtimes_{29}\mathbb Z_4$, where
$\mathrm{Dic}_3=C_3\rtimes C_4$ with $C_4$ acting on $C_3$ by inversion, and $\mathbb Z_{84}\rtimes_{29}\mathbb Z_4$
is $\langle\xi,\sigma\mid\xi^{84}=\sigma^4=1,\ \sigma\xi\sigma^{-1}=\xi^{29}\rangle$.

*Proof.* Put $u=x^4$, $v=x^3$, $s=y^{21}$ and $t=y^8$, of orders 3, 4, 4 and 7; then $x=uv^{-1}$ and $y=st$ (as
$4-3=1$ and $21+8\equiv1\pmod{28}$). Conjugation by $y$ raises $u$ to the power $5\equiv-1\pmod3$ and fixes $v$
($5\equiv1\pmod4$), and $t$ acts by $5^8\equiv1\pmod{12}$; so $v$ and $t$ are central, $s$ inverts $u$, and
$s^2=y^{42}$ centralizes $u$. Hence $H=\langle u,s\rangle\cong C_3\rtimes C_4=\mathrm{Dic}_3$, of order 12 (as
$\langle x\rangle\cap\langle y\rangle=1$), and $\langle v\rangle$, $\langle t\rangle$, $H$ commute pairwise and
generate $G$. The orders multiply to $4\cdot7\cdot12=336=|G|$, so $G=\langle v\rangle\times\langle t\rangle\times H
\cong C_{28}\times\mathrm{Dic}_3$. For the second group, $29\equiv-1\pmod3$, $29\equiv1\pmod4$, $29\equiv1\pmod7$ and
$29^2\equiv1\pmod{84}$, and the same argument with $\xi^{28}$, $\xi^{21}$, $\xi^{12}$ and $\sigma$ gives
$C_4\times C_7\times(C_3\rtimes C_4)$. $\square$

GAP's `IdGroup` and Magma's `IdentifyGroup` both identify $G$ as SmallGroup(336,81) (computed in our verification
of 2026-11-11, not by `src/rebuild.py`).

**Lemma 7 (parity).** Every column of $H_X$ and of $H_Z$ has weight 3, so the sum of the rows of $H_X$, and the sum
of the rows of $H_Z$, is $\mathbf 1$; hence every logical operator of either side has even weight. (Lemma 1 of
`two-block-and-amplified-codes/README.md`; Lemma 4 above.)

Computed (`src/rebuild.py`): the group law (associativity by Light's test on the generators $x$, $y$, identity,
inverses, the relations of the presentation), $\operatorname{rank}H_X=\operatorname{rank}H_Z=322$, so $k=28$; every
check has weight $6=3+3$, every qubit lies in 3 X-checks and 3 Z-checks; the Tanner graph is connected; the
witnesses. The automorphism group of the Tanner graph with qubits, X-checks and Z-checks as fixed vertex classes
(Magma, in the search job) has order 168; Magma identifies it as SmallGroup(168,32), which GAP describes as
$C_{14}\times\mathrm{Dic}_3$ (computed in our verification); it acts semiregularly on the qubits with the four orbits
$\{0,\dots,167\}$, $\{168,\dots,335\}$, $\{336,\dots,503\}$, $\{504,\dots,671\}$, two in each block. Its three
generators are in `records/` (files `.gensA`, `.gensB`, `.gensM`); `src/rebuild.py` checks that each maps the rows
of $H_X$ onto rows of $H_X$ and the rows of $H_Z$ onto rows of $H_Z$, and that they generate a group of order 168
with these orbits.

### Relation to the coset codes of Qian and Li

The group $C_{28}\times\mathrm{Dic}_3$ (isomorphic to $\mathbb Z_{84}\rtimes_{29}\mathbb Z_4$, Lemma 6) is also the
host group of Qian and Li, *Multi-agent discovery of practical quantum LDPC codes*, arXiv:2608.08996, who build
coset codes over it, among them $[[336,12,20]]$ ($w=6$, $A=e+xs+x^3$, $B=s^3+x^2+x^4$, $K=\langle x^{42}\rangle$).
Our code is not the double cover ($K$ trivial) of their $[[336,12,20]]$, nor a double cover of any two-block code
with $k=12$ (Lemma 8 and the computation below).

**Lemma 8.** Let $G$ be a finite group, $a,b\in\mathbb F_2[G]$ and $K$ a central subgroup of $G$ whose order is a
power of 2. Let $k(G)$ be the number of logical qubits of $C(G;a,b)$ (in the convention of Part B) and $k(G/K)$ that
of $C(G/K;\bar a,\bar b)$, where $\bar a,\bar b$ are the images of $a,b$ in $\mathbb F_2[G/K]$. Then
$k(G/K)\le k(G)\le|K|\,k(G/K)$.

*Proof.* Put $R=\mathbb F_2[G]$ and let ${}^*$ be the anti-automorphism of $R$ induced by $g\mapsto g^{-1}$ (a
linear bijection with $(xy)^*=y^*x^*$). $L_a$ is the left multiplication by $a$ and $R_b$ the right multiplication by
$b$, so $\operatorname{rank}H_X=\dim(aR+Rb)$. The transposes $R_b^{\mathsf T}$, $L_a^{\mathsf T}$ are the right
multiplication by $b^*$ and the left multiplication by $a^*$, so $\operatorname{rank}H_Z=\dim(Rb^*+a^*R)$, and
applying ${}^*$ once more gives $\operatorname{rank}H_Z=\dim(bR+Ra)$. Hence $k(G)=\dim M_1+\dim M_2$ with
$M_1=R/(aR+Rb)$ and $M_2=R/(bR+Ra)$. As $K$ is central, $M_1$ and $M_2$ are modules over $A=\mathbb F_2[K]$; let
$I$ be its augmentation ideal. Since
$R/IR=\mathbb F_2[G/K]$, $M_1/IM_1$ and $M_2/IM_2$ are the two corresponding quotients of $\mathbb F_2[G/K]$, so
$k(G/K)=\dim M_1/IM_1+\dim M_2/IM_2$, which gives the first inequality. $K$ is a 2-group, so $I$ is nilpotent and $A$
is local with residue field $\mathbb F_2$; by Nakayama's lemma $M_i$ is generated by $\dim M_i/IM_i$ elements, so
$\dim M_i\le|K|\dim M_i/IM_i$, which gives the second. $\square$

Computed (`672-28-24/related/k_covers.py`, standard library only, GF(2) ranks; its output is
`672-28-24/related/k_covers.out`, and `make related` reruns it and compares): the two-block code of $A$, $B$ over
$\mathbb Z_{84}\rtimes_{29}\mathbb Z_4$ has $k=12$ for $K=\langle x^{42}\rangle$ ($n=336$, as Qian and Li report)
and also $k=12$ for $K$ trivial ($n=672$), with $A$ acting from the left or from the right; our code has $k=28$, and
its quotients by the three central involutions $y^{14}$, $x^6$, $x^6y^{14}$ of $G$ have $k=20$. So the double cover
of their code has $k=12\ne28$, and by Lemma 8 with $|K|=2$ every double cover of a two-block code with $k=12$ has
$k\le24$.

## Novelty and attribution

Construction and method: each code of Part A is a CPM lift of a sparse $6\times14$ template (check weight 6, column
degrees 2 and 3; Definition 1, whose six columns of degree 2 form a connected unicyclic graph on the rows of each
side) whose X/Z row overlaps are made even by pair partitions with the girth rule of Okada and Kasai,
*Pair-Partition Constructions for CPM-Based Quantum LDPC Codes*, arXiv:2607.14091; their codes and Kasai's catalogue
(kasai.ict.eng.isct.ac.jp/pair_partition_cpm_css_codes_20260714.html, 44 codes) use all-one $J\times L$ templates
with check weight $L\ge8$, as do the follow-ups arXiv:2609.30069 and arXiv:2609.35601. Quasi-cyclic CSS codes from
circulant permutation matrices go back to Hagiwara and Imai, *Quantum Quasi-Cyclic LDPC Codes* (ISIT 2007,
arXiv:quant-ph/0701020). Appendix A of arXiv:2609.30069 applies to the codes of Okada and Kasai the distance bounds
of MacKay and Davey and of Smarandache and Vontobel (IEEE Trans. Inf. Theory 58 (2012) 585-607; Theorems 7 and 8 of
arXiv:0901.4129v2); the cofactor words of the latter are the structural upper bounds of our search.

Relation to the board, export `82a52ee51`: CPM lifts whose X/Z row overlaps split into pairs of equal exponent
difference are common on the board at check weights 4 to 8, among them two-block and bicycle codes, the
$[[6N,2N+2,6]]$ family of arXiv:2610.19451 (Nourozi, Smarandache, Mitchell; Construction I, rate about 1/3; all-one
$2\times6$ template, $d=6$; on the board the eight entries `186-64-6` to `768-258-6`, whose check matrices we
checked) and irregular templates such as `260-58-6` and the $14\times28$ lifts `280-20-8` to `700-50-8` at check
weight 6; entries tagged `pair-partition-cpm` occur only at check weight $\ge7$ (e.g. `990-198-17`, `994-142-40`,
`972-166-34`). In the layouts tested (consecutive or interleaved blocks of the stored qubit order), no entry of
check weight 6 is a lift of a template with column degrees 2 and 3 (computed with `board_scan/cpm_scan.py` on the
`codes/` directory of the board at `82a52ee51`, as `python3 board_scan/cpm_scan.py <codes directory>`; its output,
one line per entry with a CPM layout, 204 entries, is `board_scan/cpm_scan.out`). `990-198-17` (MicheldeMessieres)
enumerates one-side templates with columns of degree 2 or 3 (exhaustively for $J=3,4$, sampled for $J=5$), ranks
them by its generalization of Appendix A of arXiv:2609.30069 to sparse templates, records that its pairings with 5
free dimensions carry a forced weight-5 logical through the two degree-2 columns, and finds no weight-6 template with
$J\le4$ that reaches its target; `994-142-40` uses $6\times14$ templates of rate 1/7 and asks for at least 3 free
exponent directions, and `972-166-34` (jacksongquera) uses the square complex of rows, columns and pairs and its free
directions (our $b_1$). Ours: the census of the $6\times14$ templates of check weight 6, the five templates,
pairings, cocycle classes and lifts, and the exact distances.

Two-block group algebra codes: Lin and Pryadko, *Quantum two-block group algebra codes*, arXiv:2306.16400
(Phys. Rev. A 109, 022407, 2024). For `672-28-24` see "Relation to the coset codes of Qian and Li" above.

Novelty claim of the six entries: new parameters only. No published code and no board entry has these parameters,
and each code dominates the board entries stated in its note and pull request. A literature search on 2026-11-11
(arXiv, with the papers above and their follow-ups; Kasai's catalogue; our catalogue of the code tables of 76 papers;
the board at `82a52ee51`, its research notes and open pull requests; web searches for the parameters) found no source
for any of the six parameter sets. Not checked: QECdb, the full text of Google Scholar, journal-only publications.

## Files

Per code, `<slug>/`:

- `definition.json`: the definition (Part A: `template`, `exponents`, `pairing`, `P`, the lift rule; Part B: the
  group, `a`, `b`, the convention), the parameters, the identifiers of our search job.
- `hx.txt`, `hz.txt`: the check matrices, one check per line, 0-based qubits, rows in the order of Definition 1
  (Part A) or of $\mathrm{index}(c)$ (Part B). These are the checks of the board entry.
- `witness_X.txt`: an X-logical of weight $d$ (in $\ker H_Z$, outside $\operatorname{rowsp}H_X$);
  `witness_Z.txt`: a Z-logical of weight $d$.
- `records/`: the instance files of the exact searches (`<id>.hx`, `<id>.hz`: a header "n m", then one check per
  line), the generator files (`.gensA` in exactA's format, a count line and one permutation per line; `.gensB` and,
  where present, `.gensM` in exactB's format, with a comment line), the witnesses `.wX`, `.wZ`, the four logs
  `{A,B}_{X,Z}_W<W>.log`, for `742-106-15` and `672-28-24` the jobs' records `{A,B}_{X,Z}.json`, `prep.json` and
  `autgroup.log`; `searches.json` summarizes the four runs (made from the logs by our assembly script and checked
  against them by `src/rebuild.py`). `<id>` is the identifier of the code in our search job.
- `third_tools/`: `mtx/` (inputs), `dist-m4ri/` (`.out` result line, `.err` diagnostics, or one `.log` with both;
  `.nz` codeword lists), `tropidist/` (`.json` records, console output), `runs.json` (one entry per run: program,
  side, files, result, whether a symmetry was used, the first weight searched; dist-m4ri's timing runs on
  `672-28-24` without result lines are listed under `calibration`).
- `672-28-24/related/`: `k_covers.py` and its output `k_covers.out` (Lemma 8 and the computation after it).

`src/rebuild.py` (the checks of `make test`), `src/rerun.py` and `src/same_results.py` (`make rerun`),
`src/programs.sha256`; `logs/rebuild.log`, `logs/rerun.log` (our outputs of `make test` and `make rerun`) and
`logs/environment.txt`; `SHA256SUMS` (all files).

In the copied logs, the absolute path of exactB in its `cmd:` line was replaced by `bin/exactB`, and the absolute
paths of our verification work directories were cut to paths relative to the work directory (whose `mtx/` is
`third_tools/mtx/` here); apart from that the logs, records and tool outputs are unchanged copies.

## Credits and licenses

The template search, the six codes, the programs, the computations and the texts of this directory were made with
Claude Opus 5.5 (Anthropic) in Claude Code. Pair partitions and the girth rule are due to Okada and Kasai
(arXiv:2607.14091), the two-block group algebra codes to Lin and Pryadko (arXiv:2306.16400). dist-m4ri is by
Pryadko and Zeng, tropidist by Cheng, Guo and Pan; only inputs written by us and the outputs of these programs are
included. Everything in this directory is under the MIT license of the repository.
