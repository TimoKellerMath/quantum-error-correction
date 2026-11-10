# Distance certificates for entries of the qLDPC challenge board (November 2026)

Evidence for one pull request to the qLDPC challenge board (https://github.com/unitaryfoundation/qldpc-challenge,
entries `codes/<slug>.json` at commit `8042f9da`, 2026-11-09): 163 distance certificates `certs/<slug>.json` for
entries of other authors (160 exact, 3 partial), the correction of an over-claimed side value in four entries,
and the correction of one name field. This directory holds the four programs, the proofs of their methods, the
recorded runs and a driver that rebuilds every instance from the board's `codes/<slug>.json` and reruns the
checks. Programs, computations and first drafts of the texts were made with Claude Opus 5.5 (Anthropic) in
Claude Code.

Conventions (those of the board): a CSS code is a pair $H_X$, $H_Z$ over $\mathbb F_2$ with $H_XH_Z^{\mathsf T}=0$;
an X-type logical operator is a vector of $\ker H_Z\setminus\operatorname{rowsp}H_X$, a Z-type one a vector of
$\ker H_X\setminus\operatorname{rowsp}H_Z$; $d_X$, $d_Z$ are their least weights and $d=\min(d_X,d_Z)$. For a pair
$(A,B)$ with $AB^{\mathsf T}=0$ put $\delta(A,B)=\min\{|x|:x\in\ker A\setminus\operatorname{rowsp}B\}$, so
$d_Z=\delta(H_X,H_Z)$ and $d_X=\delta(H_Z,H_X)$.

## What is certified

| kind | entries | method | statement used |
|---|---|---|---|
| graph-like | 126 (124 with this method only) | `graphlike_A` and `graphlike_B` | Propositions 1 and 2 |
| exhaustive search | 7: 96-17-6, 126-18-14, 182-6-18, 185-64-5, 185-65-5, 192-12-14, 336-20-16 | `exactA` and `exactB` (185-64-5, 185-65-5 also both graph-like programs) | Proposition 4 |
| amplified | 29: 904-12-36, 910-12-36 ($d=36$), 948-24-28, 612-36-28, 630-36-28 ($d=28$), and 24 entries whose base carries a board certificate | audit `src/amp_audit.py` and the base distance | Proposition 3 |
| partial | 3: 128-12-16 ($14\le d\le16$), 628-24-32 ($28\le d\le32$), 672-20-32 ($26\le d\le32$) | `exactA` and `exactB` below the claim; for 628-24-32 Corollary 3' | Propositions 3, 4 |

`CERTIFICATES.md` lists every entry with its exact side values, the method, and the claims on the board.
In every exact case the claimed $d$ of the entry is the exact one. Four entries claim a side value above the exact
one (their $d$ comes from the other side and is unaffected): 714-22-1 ($d_Z=2$, claimed 7), 719-305-1 ($d_Z=3$,
claimed 4), 749-155-3 ($d_Z=4$, claimed 7), 749-156-3 ($d_X=4$, claimed 5); the pull request lowers these side
values with witnesses (in `data/entries.json`, field `correction`). The name field of 910-12-36 read
"[[910,12,39]]"; the entry is the $[[4,2,2]]$ amplification of 182-6-18, so $d=36$ (Proposition 3) and the name is
corrected to "[[910,12,36]]".

Nine certificates already on the board store $d$ as the value of the side with the larger distance (96-35-3,
96-36-3, 119-35-3, 119-36-3, 120-35-4, 120-37-4, 336-58-6, 896-194-6, 960-258-6), so the site, which compares
each side value with the entry's, does not honour them; both graph-like programs give the exact side values
equal to the entries' claims. Those files are not changed by the pull request.

## Statements and proofs

### Graph-like sides

**Lemma 1.** Let $AB^{\mathsf T}=0$ and assume every column of $A$ has weight at most 2. Let $G$ be the multigraph
with vertices the rows of $A$ and a vertex $\infty$, and one edge $e$ per column: it joins the two rows containing
$e$, or that row and $\infty$ if the column has weight 1, and is a loop at $\infty$ if it has weight 0. Then
$\ker A$ is the cycle space of $G$ (edge sets in which every vertex has even degree), and every
$x\in\ker A\setminus\operatorname{rowsp}B$ with $|x|=\delta(A,B)$ is a circuit (the edge set of a closed walk
without repeated vertices; a loop, two parallel edges, or a cycle).

*Proof.* For a row $i$, $(Ax)_i$ is the parity of the degree of $i$ in $x$; the degree of $\infty$ is
$2|x|-\sum_i\deg_x(i)$, so it is even when all the others are. By Veblen's theorem a set with all degrees even is
a disjoint union of circuits $C_1,\dots,C_r$, each in $\ker A$. As $x\notin\operatorname{rowsp}B$, some $C_i$ is
not in $\operatorname{rowsp}B$, so $|C_i|\ge\delta(A,B)=|x|=\sum_j|C_j|$ and $r=1$. $\square$

**Lemma 2.** Let $AB^{\mathsf T}=0$ and $y_1,\dots,y_k\in\ker B$ with
$\operatorname{rowsp}A+\langle y_1,\dots,y_k\rangle=\ker B$. For $x\in\ker A$: $x\in\operatorname{rowsp}B$ iff
$\langle x,y_j\rangle=0$ for all $j$.

*Proof.* The standard form is nondegenerate, so $\operatorname{rowsp}B=(\ker B)^\perp$; and $x\perp\operatorname{rowsp}A$
because $x\in\ker A$. $\square$

Both programs compute such $y_j$ by Gaussian elimination and check at run time that they lie in $\ker B$ and that
their number is $n-\operatorname{rk}A-\operatorname{rk}B$; write $a(e)=(\langle e,y_j\rangle)_j\in\mathbb F_2^k$ for an
edge $e$ and $a(x)=\sum_{e\in x}a(e)$.

**Proposition 1 (graphlike_A).** In the situation of Lemmas 1 and 2 with $k\ge1$, fix a total order on the
vertices. For a vertex $v$ let $G_v$ be the subgraph induced on the vertices $\ge v$ (without loops), $T_v$ a
breadth-first tree of the component of $v$ in $G_v$ with depths $\delta_v$ and tree paths $P_v(\cdot)$ from $v$,
and for an edge $e=xy$ of $G_v$ with both ends in $T_v$ put $C(v,e)=P_v(x)+e+P_v(y)$ and
$\ell(v,e)=\delta_v(x)+1+\delta_v(y)$. If a loop $e$ has $a(e)\ne0$, then $\delta(A,B)=1$; otherwise
$\delta(A,B)=\min\{\ell(v,e):a(C(v,e))\ne0\}$.

*Proof.* A loop lies in $\ker A$, and it is nontrivial iff $a(e)\ne0$ (Lemma 2). Otherwise: ($\ge$) $C(v,e)$ has
zero boundary, so it lies in $\ker A$ (Lemma 1), it is nontrivial when $a(C(v,e))\ne0$ (Lemma 2), and
$|C(v,e)|\le\ell(v,e)$. ($\le$) By Lemma 1 there is a nontrivial circuit $c$ of length $d=\delta(A,B)$, not a loop.
Let $v$ be its least vertex and $v=u_0,e_1,u_1,\dots,e_d,u_d=v$ the closed walk around it; it lies in $G_v$, so
all $u_i$ lie in $T_v$. Then $\sum_iC(v,e_i)=c$, because every $P_v(u_i)$ with $0<i<d$ occurs twice and
$P_v(v)=\emptyset$. Since $a(c)\ne0$, some $a(C(v,e_i))\ne0$, and walking along $c$ in both directions gives
$\delta_v(u_{i-1})\le i-1$, $\delta_v(u_i)\le d-i$, so $\ell(v,e_i)\le d$. $\square$

The program stops each breadth-first search at the depth that can still improve the current minimum (the
bookkeeping is Lemma 7 of `src/graphlike_A_NOTES.md`, where this proposition is Proposition 5) and rechecks the
word it returns: weight, $Ax=0$, and $x\notin\operatorname{rowsp}B$ by elimination against $B$ (independent of the
$y_j$).

**Proposition 2 (graphlike_B).** In the situation of Lemmas 1 and 2, let $S_j=\operatorname{supp}y_j$ and let
$\widetilde G_j$ be the double cover of $G$ with vertices $(u,b)$, $b\in\mathbb F_2$, in which every edge $e$ from
$u$ to $w$ has the two lifts from $(u,b)$ to $(w,b+[e\in S_j])$. For any ordering $v_1,\dots,v_N$ of the vertices
let $\widetilde G_{j,i}$ be the double cover of the subgraph induced on $\{v_i,\dots,v_N\}$. Then
$$\delta(A,B)=\min_{i,j}\operatorname{dist}_{\widetilde G_{j,i}}\big((v_i,0),(v_i,1)\big).$$

*Proof.* A path from $(v,0)$ to $(v,1)$ of length $L$ projects to a closed walk at $v$ that uses edges of $S_j$ an
odd number of times; the set $x$ of edges it traverses an odd number of times lies in $\ker A$, has
$\langle x,y_j\rangle=1$, hence is nontrivial (Lemma 2), and $|x|\le L$. Conversely, a nontrivial $x$ of weight
$\delta(A,B)$ is a circuit $C$ (Lemma 1) with $\langle C,y_j\rangle=1$ for some $j$ (Lemma 2); if $v_i$ is its first
vertex in the ordering, $C$ lies in the subgraph on $\{v_i,\dots,v_N\}$, and traversing it from $v_i$ lifts to a
path of length $|C|$ from $(v_i,0)$ to $(v_i,1)$. $\square$

The program runs the $k$ breadth-first searches of one root at once, one bit per cover, and rechecks its word as
above (`src/graphlike_B_NOTES.md`, Corollary 4 and Propositions 8, 9, Lemma 10). The two programs share no code:
graphlike_A follows minimal circuits through annotated trees, graphlike_B searches the covers.

### Amplified entries

**Definition.** For a CSS code $(H_X,H_Z)$ with $m_X$ X-checks and $m_Z$ Z-checks (rows not necessarily independent)
and $N\ge2$ even, $\mathrm{Amp}_N(H_X,H_Z)$ has the central qubits $(q,a)$, $q\in[n]$, $a\in[N]$, a side qubit
$s_c$ for every X-check $c$ and a side qubit $t_z$ for every Z-check $z$; its X-checks are
$(c,a)=\{(q,a):q\in H_X[c]\}\cup\{s_c\}$ and $(q)^X=\{(q,a):a\in[N]\}\cup\{t_z:q\in H_Z[z]\}$, its Z-checks
$(z,a)=\{(q,a):q\in H_Z[z]\}\cup\{t_z\}$ and $(q)^Z=\{(q,a):a\in[N]\}\cup\{s_c:q\in H_X[c]\}$. This is the
central truncation of the tensor product of the base complex with the $[[N,N-2,2]]$ complex
$\mathbb F_2\xrightarrow{\mathbf 1^{\mathsf T}}\mathbb F_2^N\xrightarrow{\mathbf 1}\mathbb F_2$ (Z. Liang, B. Gu, Y.-A. Chen,
J. Eisert, Z. Wang, *Ultra-high-distance quantum memories from amplified qLDPC codes*, arXiv:2609.37231, Methods A,
Eqs. (14), (15)). A central vector $z$ is written $\sum_az_a\otimes e_a$ with $z_a\in\mathbb F_2^n$.

**Lemma 3.** $\mathrm{Amp}_N(H_X,H_Z)$ is a CSS code, and for $N\ge4$ even,
$[u]\otimes[\ell]\mapsto[u\otimes\ell]$ (the vector $u\otimes\ell$ placed on the central qubits) is an isomorphism
$(\ker H_X/\operatorname{rowsp}H_Z)\otimes(\mathbf 1^\perp/\langle\mathbf 1\rangle)\to\ker\widetilde H_X/\operatorname{rowsp}\widetilde H_Z$,
and dually for X-type operators.

*Proof.* Commutation is a count: $(c,a)$ and $(z,a')$ overlap in $|H_X[c]\cap H_Z[z]|$ qubits, an even number, if
$a=a'$ and not at all otherwise; $(c,a)$ and $(q)^Z$ share $(q,a)$ and $s_c$ if $q\in H_X[c]$ and nothing otherwise;
$(q)^X$ and $(z,a)$ likewise with $t_z$; $(q)^X$ and $(q')^Z$ share the $N$ central qubits $(q,a)$ if $q=q'$, an even
number, and nothing otherwise. The quotient is the homology in the middle degree of the truncated complex, which equals
that of the full tensor product; by the Künneth formula over a field it is
$\bigoplus_{i+j=2}H_i\otimes H_j$, and the amplifier complex has homology only in the middle,
$\mathbf 1^\perp/\langle\mathbf 1\rangle$, of dimension $N-2$, because $\mathbf 1\ne0$ and $\langle\mathbf 1,\mathbf 1\rangle=N=0$. $\square$

**Lemma 4.** For $\omega\in\mathbb F_2^N$ of even weight put $c_\omega(z)=\sum_{a:\omega_a=1}z_a$. Then
$c_\omega(\ker\widetilde H_X)\subseteq\ker H_X$, $c_\omega(\operatorname{rowsp}\widetilde H_Z)\subseteq\operatorname{rowsp}H_Z$
and $c_\omega(u\otimes\ell)=\langle\omega,\ell\rangle u$.

*Proof.* If $\widetilde H_Xz=0$, the check $(c,a)$ gives $(H_Xz_a)_c=z(s_c)$ for every $a$; summing over the
$|\omega|$ indices with $\omega_a=1$ gives $(H_Xc_\omega(z))_c=|\omega|\,z(s_c)=0$. The check $(z,a)$ contracts to
$H_Z[z]$ or $0$, and $(q)^Z$ to $|\omega|e_q=0$. The last claim is the definition. $\square$

**Proposition 3.** Let $(H_X,H_Z)$ be a CSS code with $k\ge1$ and $N\ge4$ even. Then
$d_Z(\mathrm{Amp}_N(H_X,H_Z))=2d_Z$ and $d_X(\mathrm{Amp}_N(H_X,H_Z))=2d_X$. In particular, if a code $E$ has the
row spaces of $\mathrm{Amp}_N(S_X,S_Z)$ for some $S_X$, $S_Z$ with $\operatorname{rowsp}S_P=\operatorname{rowsp}H_P(B)$
($P=X,Z$), then $d_P(E)=2d_P(B)$.

*Proof.* Lower bound: let $z$ be a nontrivial Z-type logical of the amplified code. Fix representatives
$\ell_1,\dots,\ell_{N-2}$ of a basis of $\mathbf 1^\perp/\langle\mathbf 1\rangle$; by Lemma 3,
$[z]=\sum_j[u_j\otimes\ell_j]$ with some $[u_i]\ne0$. The form induces a perfect pairing on
$\mathbf 1^\perp/\langle\mathbf 1\rangle$, so there is $\omega\in\mathbf 1^\perp$ with $\langle\omega,\ell_j\rangle=\delta_{ij}$;
$\omega'=\omega+\mathbf 1$ also has even weight and the same pairings. By Lemma 4, $c_\omega(z)$ and $c_{\omega'}(z)$
both represent $[u_i]\ne0$, so each has weight at least $d_Z$; as the supports of $\omega$ and $\omega'$ partition
$[N]$, $2d_Z\le|c_\omega(z)|+|c_{\omega'}(z)|\le\sum_a|z_a|\le|z|$. Upper bound: for $u$ of weight $d_Z$ and
$\ell=e_0+e_1\notin\langle\mathbf 1\rangle$ (as $N\ge4$), $u\otimes\ell$ is nontrivial by Lemma 3 and has weight
$2d_Z$. The X statement follows by exchanging the roles of X and Z, which the definition allows. The last
sentence holds because distances depend only on row spaces. $\square$

This is Lemma 4 of arXiv:2609.37231 for this amplifier (two disjoint representatives per logical class), proved
there by the same contraction. **Corollary 3'.** If $a\le d_P(B)\le b$, then $2a\le d_P(E)\le2b$; with
$d_P(128\text{-}12\text{-}16)\in\{14,16\}$ this gives $d_P(628\text{-}24\text{-}32)\in\{28,32\}$.

The audit `src/amp_audit.py` reads, from the rows of the entry $E$ that consist of one copy of a base support and
one side qubit, the supports of $S_X$ and $S_Z$ with their side qubits, checks
$\operatorname{rowsp}S_P=\operatorname{rowsp}H_P(B)$ by ranks, rebuilds $\mathrm{Amp}_N(S_X,S_Z)$ by the definition
with those side qubits and compares row sets. For all 30 amplified entries here (29 certified, 628-24-32 partial)
the rows are identical.

### Exhaustive search

[ED] is `weight6-exact-distances/exact_distances.pdf` in this repository (T. Keller, *Exact distances of five
weight-six two-block codes*), whose Algorithm 1, Proposition 3.1 and Corollary 3.2 state and prove the search for
check weight 6. `exactA` generalizes its implementation A, `exactB` its implementation B, to any check weights,
any $k$ and any group of qubit permutations; `src/exactA_NOTES.md` (Lemma A, Corollary B, Lemmas C to E and G) and `src/exactB_NOTES.md`
(Lemmas N1 to N6) give the changes and their proofs. The search for the side with check matrix $A$ (logicals in
$\ker A\setminus\operatorname{rowsp}B$) rests on four facts.

**Proposition 4.** Let $AB^{\mathsf T}=0$, $k\ge1$, $d=\delta(A,B)$, and let $\Gamma$ be a group of permutations of
the qubits mapping the set of rows of $A$ onto itself and the set of rows of $B$ onto itself, with orbits
$O_1,\dots,O_r$ and representatives $q_i\in O_i$.

1. If $x$ is nontrivial of weight $d$ and $s\in\operatorname{rowsp}B$, then $|x\cap s|\le|s|/2$ (the majority rule).
2. If $\mathbf 1\in\operatorname{rowsp}A$, every $x\in\ker A$ has even weight (the parity reduction).
3. Some nontrivial $x$ of weight $d$ contains some $q_i$ and avoids $O_1\cup\dots\cup O_{i-1}$ (orbit anchors).
4. ([ED], Proposition 3.1, Corollary 3.2) If $d\le W$ and some nontrivial $x$ of weight $d$ contains
   $T_0=\{q_i\}$ and avoids $X_0=O_1\cup\dots\cup O_{i-1}$, the search started from $(T_0,X_0)$ with bound $W$
   reports a nontrivial word of weight at most $W$. Hence, by item 3, if none of the $r$ runs reports anything,
   then $d>W$.

*Proof.* (1) $x+s$ is nontrivial and $|x+s|=|x|+|s|-2|x\cap s|\ge|x|$. (2) $\mathbf 1$ is a sum of rows of $A$, each
orthogonal to $x$, and $\langle\mathbf 1,x\rangle=|x|\bmod2$. (3) Take $x_0$ nontrivial of weight $d$, $i$ least with
$x_0\cap O_i\ne\emptyset$, $q\in x_0\cap O_i$ and $\gamma\in\Gamma$ with $\gamma(q)=q_i$; $\gamma$ preserves
$\ker A$, $\operatorname{rowsp}B$, weights and every orbit. (4) The proof of [ED], Proposition 3.1 uses only that
every row of $A$ meets $x$ evenly, item 1 for the members of the majority family (rows of $B$ and sums of two
overlapping rows of $B$), and that the lower bound used for pruning never exceeds $|x\setminus T|$; none of this
depends on the check weight, on $k$ or on the group (`src/exactA_NOTES.md`, "What changes in the proof";
`src/exactB_NOTES.md`, Lemmas N1, N2, N5). $\square$

Both programs verify that the generators they are given preserve the row sets, recheck every word they report
without the logical basis, and print the parity decision; a run at $W=d$ that finds a word of weight $d$ supplies
the upper bound, the entry's own witness supplies it as well.

The results are computer-assisted: two independently written programs agree on every certified value (four on
185-64-5 and 185-65-5), and their methods are proved above. This makes an undetected programming error unlikely
but does not exclude it; none of the certificates is a replayable proof object.

## What was run

- `src/certify.py --record` on the board at `8042f9da` (log `logs/record_run.jsonl`, 23 s on one core): both
  graph-like programs on the 126 graph-like entries; exactA and exactB at $W=d_P-1$ and $W=d_P$ on both sides of
  the 7 entries of the exhaustive-search kind (among them the bases 182-6-18, 192-12-14, 126-18-14 of the
  round-two amplified entries) and of the 22 bases of the other amplified entries, which carry board
  certificates (all 22 agree with them); the audit of the 30 amplified entries. Every recorded result is in `data/results/<slug>.json`.
- The partial certificates of 128-12-16 and 672-20-32 rest on longer searches, recorded by the tracks that found
  them (`logs/tracks/`, about 38 CPU-minutes in total): for 128-12-16 at $W=12$ on relabelled copies of the board
  matrices (exactA 428 s and 410 s, exactB 285 s and 272 s for the X and Z sides), confirmed by a third run of
  exactA in another mode by an independent check; for 672-20-32 at $W=25$ (weights $\le24$ searched by parity) on
  the board matrices (exactA 202 s and 205 s, exactB 247 s and 229 s). `logs/long_run.jsonl` is the rerun of
  both by the driver of this package (`--long`) on the board labelling with the groups of `data/groups/` (31
  minutes, same verdicts; for 672-20-32 with the same node counts as the original runs, which used the same
  labelling and generators).
- Validation of the programs before these runs: graphlike_A on 33 surface codes of known distance, in 14223
  comparisons with exhaustive enumeration on random graph-like codes and against all 302 exact board certificates
  of the board entries with two graph-like sides; graphlike_B on toric, planar and rotated surface codes, on 800
  random graph-like codes against a minimum-cycle-basis check (706 also against exhaustive enumeration) and with
  an unpruned reference search on all 456 board entries with a graph-like side; exactA and exactB on twelve
  board codes with certificates (check weights 4 to 12, relabelled, both sides: nothing at $W=d-1$, resp. $d-2$
  under parity, a word of weight $d$ at $W=d$) and on small codes against exhaustive enumeration.

## Reproduce

Requirements: a C compiler, `make`, Python 3.8 or later (standard library only); nauty's `dreadnaut` only to
recompute the group generators of `data/groups/` (`src/tanner_group.py`).

```
make                                                    # bin/graphlike_A, graphlike_B, exactA, exactB
git clone https://github.com/unitaryfoundation/qldpc-challenge board
python3 src/certify.py --board board                    # every entry except the long searches, under a minute
python3 src/certify.py --board board --slug 182-6-18 904-12-36
python3 src/certify.py --board board --long --slug 128-12-16 672-20-32   # about 40 CPU-minutes
```

The driver identifies each entry by the board's binding hash of its matrices (`checks_sha256` in
`data/entries.json`), so it also runs on later commits of the board as long as the matrices are unchanged; it
compares verdicts and values (not node counts or times) with `data/results/<slug>.json` and exits with status 1
on any difference. `make check BOARD=board` runs the second command.

## Layout

- `src/`: the four programs (`graphlike_A.c`, `graphlike_B.c`, `exactA.c`, `exactB.c`, each with its design notes
  `*_NOTES.md`, written while the programs were developed; they cite a research brief and paths of the working
  tree that are not public, and the proofs needed here are restated above), `certify.py` (the driver),
  `amp_audit.py` (the audit of amplified entries), `tanner_group.py` (group generators with dreadnaut).
- `data/entries.json`: every entry: category, methods, $n$, $k$, binding hash, exact side values or $d_{\rm lower}$,
  base and $N$ for amplified entries, the side corrections with their witnesses; `file` is false for the 22
  bases that are only rechecked. `data/results/<slug>.json`: the recorded outputs of both programs per entry
  (JSON lines of the graph-like programs with witnesses, parsed RESULT lines and words of the searches, the
  audits). `data/groups/<slug>.gens`: Tanner-graph group generators.
- `CERTIFICATES.md`: the table of all entries.
- `logs/`: `record_run.jsonl` (the run that wrote `data/results/`), `long_run.jsonl`, `cleantest.log` (the driver
  on five entries in a clean copy of this directory), `match_board.log` (the tower and quotient codes of our
  other computations matched against the board by canonical forms of their Tanner graphs), `environment.txt`,
  `tracks/` (logs of the long searches).
- `SHA256SUMS`: checksums of all files.

## Credits and license

The amplification and its distance lemma: Liang, Gu, Chen, Eisert and Wang (arXiv:2609.37231), the overlap bound
there after Audoux and Couvreur (arXiv:1512.07081). The exhaustive search refines the complete searches of Dumer,
Kovalev and Pryadko (arXiv:1611.07164) and Okada and Kasai (arXiv:2607.14091); [ED], Section 3 credits further
ingredients. Shortest nontrivial cycles: the idea of Horton's candidate set (J. D. Horton, SIAM J. Comput. 16
(1987)) and of double covers for shortest odd cycles. The board entries and their matrices are those of their
authors (see `CERTIFICATES.md`). MIT license (see `../LICENSE`).
