# Folded stabilizer codes for the qLDPC challenge board

Five non-CSS stabilizer codes submitted to the [qLDPC challenge board](https://github.com/unitaryfoundation/qldpc-challenge)
(one pull request per code), with their generators, the proofs used, programs and logs. Each is the fold
(symplectic halving) of a two-block CSS board entry along a fixed-point-free ZX-duality $h\mapsto c-h$ that swaps
the two blocks (Burton, Durso-Sabina, Brown, [arXiv:2406.09951](https://arxiv.org/abs/2406.09951), cited as
[BDB24]; Lee, Okada, Maskara, Kasai, Zhou, [arXiv:2609.30069](https://arxiv.org/abs/2609.30069), Def. 3, cited as
[LOMK26]); equivalently, each is a mirror code of Khesin and Lu ([arXiv:2603.05496](https://arxiv.org/abs/2603.05496),
cited as [KL26]). Every distance is exact, by three independent exhaustive computations: tropidist in the board's
Pauli encoding, a meet-in-the-middle search of ours, and Magma.

| board slug | code | parent (board entry) | $c$ | $w$ | qubit degree | $Y$ factors | $kd^2/n$ | stabilizer entries it dominates |
|---|---|---|---|---|---|---|---|---|
| `60-9-8` | $[[60,9,8]]$ | `120-18-12` (G. Lee) | $(1,0,0,0)$ | 8 | 8 | 12 | 9.60 | `63-9-7`, `81-9-8` |
| `56-6-8-s` | $[[56,6,8]]$ | `112-12-8` (@vprusso) | $(0,1)$ | 6 | 6 | 0 | 6.86 | `65-5-8`, `72-6-8-b` |
| `32-9-6` | $[[32,9,6]]$ | `64-18-8` (@MathysRennela) | $(1,1,0)$ | 8 | 8 | 0 | 10.12 | `36-7-6`, `40-9-6` |
| `48-6-8-s` | $[[48,6,8]]$ | `96-12-9` (@MathysRennela) | $(0,1,0)$ | 7 | 7 | 16 | 8.00 | `50-6-8` |
| `33-2-7` | $[[33,2,7]]$ | `66-4-10` (@MathysRennela) | $(0,0)$ | 6 | 6 | 9 | 2.97 | `34-2-7-s` |

$w$ is the largest Pauli weight of a generator (a $Y$ counts once), the qubit degree the largest number of
generators acting on a qubit. Domination is the board's rule among stabilizer entries in the unrestricted locality
class, with the claimed distances of board commit `8042f9da`; the dominated entries `63-9-7`, `81-9-8`, `65-5-8`,
`36-7-6`, `40-9-6`, `50-6-8` and `34-2-7-s` carry layouts (`locality`) of the local-2d-bilayer class and keep their
places in those cells. The slugs `56-6-8` and `48-6-8` belong to CSS entries, so these two codes carry the suffix
`-s`, as the board's stabilizer entries `42-6-7-s` and `34-2-7-s` do. Stabilizer and CSS codes are ranked on
separate boards; none of the five codes is locally Clifford equivalent to a CSS code (Lemma 6).

## Reproduce

Requirements: Python 3.8 or later (standard library only), a C compiler, `make`. Optional: Magma (for `make magma`),
GAP 4 (for `make gap`).

```
make check              # rebuild the five codes from their parents and check them (under 1 s); log: logs/rebuild.log
make quick              # compile src/mitm.c and rerun the 8 meet-in-the-middle searches, compared with logs/mitm/
make magma              # Magma's MinimumDistance on the five codes (about 90 CPU-seconds), compared with logs/magma/
make gap                # GAP's numbering of the groups of four parents, compared with logs/gap_elements.log
make tropidist-inputs   # the inputs of our tropidist runs, compared with logs/tropidist/inputs.sha256
```

`python3 src/rebuild.py --board DIR`, with `DIR` a checkout of the board repository, also compares the parents in
`data/parents/` with the board's files, and the five codes with their board entries once these exist. Our runs used
gcc 16.2.1 with `-O3 -march=native` on x86-64, Magma V2.30-1 and GAP 4.16.1, one thread per job.

## Construction

**Definition 1 (parents).** For a finite abelian group $G=\mathbb Z/o_1\times\dots\times\mathbb Z/o_r$ and subsets
$F_A,F_B\subset G$ the two-block code has the qubits $(1,h)$ and $(2,h)$, $h\in G$, with indices $\iota(h)$ and
$|G|+\iota(h)$, where $\iota$ numbers $G$ in lexicographic order of the coordinates (for $\mathbb Z/2\times\mathbb
Z/2\times\mathbb Z/3\times\mathbb Z/5$: $\iota(a,b,e,f)=30a+15b+5e+f$). The X-check $g\in G$ acts on
$(1,g+F_A)\cup(2,g+F_B)$, the Z-check $g$ on $(1,g-F_B)\cup(2,g-F_A)$; this is the board's convention
$H_X=(A\mid B)$, $H_Z=(B^{\mathsf T}\mid A^{\mathsf T})$. The board files of the parents number the group elements
differently: `112-12-8` as $(i,j)\mapsto2i+j$ in $\mathbb Z/28\times\mathbb Z/2$ (which is $\iota$), the other four
in the order of GAP's `Elements(SmallGroup(N,id))`. `data/parents/<parent>.json` holds a copy of the board's checks
and witnesses and the list `labels.elements`: the group element (coordinates with respect to GAP's
`IndependentGeneratorsOfAbelianGroup`, written by `src/gap_elements.g`, `logs/gap_elements.log`) of every board
index. With this relabelling the board's rows and those of Definition 1 coincide as sets (`src/rebuild.py`).

| parent | $G$ | $F_A$ | $F_B$ | board numbering |
|---|---|---|---|---|
| `120-18-12` | $\mathbb Z/2\times\mathbb Z/2\times\mathbb Z/3\times\mathbb Z/5$ | $(0,0,0,0),(0,0,1,1),(0,1,0,4),(1,0,2,4)$ | $(0,0,0,0),(0,1,0,2),(0,1,2,4),(1,0,1,0)$ | `SmallGroup(60,13)` |
| `112-12-8` | $\mathbb Z/28\times\mathbb Z/2$ | $(0,0),(25,1),(27,0)$ | $(0,0),(8,0),(26,0)$ | $2i+j$ |
| `64-18-8` | $\mathbb Z/2\times\mathbb Z/4\times\mathbb Z/4$ | $(0,0,0),(0,0,3),(0,3,0),(1,3,3)$ | $(0,0,0),(0,1,2),(0,3,3),(1,2,1)$ | `SmallGroup(32,21)` |
| `96-12-9` | $\mathbb Z/3\times\mathbb Z/4\times\mathbb Z/4$ | $(0,0,0),(1,3,0),(2,3,3)$ | $(0,0,0),(0,3,0),(1,2,3),(1,3,2)$ | `SmallGroup(48,20)` |
| `66-4-10` | $\mathbb Z/3\times\mathbb Z/11$ | $(0,0),(1,7),(2,0)$ | $(0,0),(1,6),(2,7)$ | `SmallGroup(33,1)` |

**Definition 2 (the dualities $\Psi_c$).** For $c\in G$ let $\Psi_c$ be the permutation of the parent's qubits that
exchanges $(1,h)$ and $(2,c-h)$ for every $h\in G$. It is a fixed-point-free involution, and it maps the support of
the X-check $g$, $(1,g+F_A)\cup(2,g+F_B)$, onto $(2,c-g-F_A)\cup(1,c-g-F_B)$, the support of the Z-check $c-g$. So
$\Psi_c$ is a ZX-duality for every $c$.

**Definition 3 (fold).** Let $C=(H_X,H_Z)$ be a CSS code on $2m$ qubits and $\Psi$ a fixed-point-free involution of
the qubits that maps the set of supports of the X-checks onto the set of supports of the Z-checks. Let
$\sigma_0<\dots<\sigma_{m-1}$ be the qubits $q$ with $q<\Psi(q)$. For $a\in\mathbb F_2^{2m}$ let $\varphi(a)$ be the
Pauli operator on $m$ qubits with $X$ on $\{i:a_{\sigma_i}=1\}$ and $Z$ on $\{i:a_{\Psi(\sigma_i)}=1\}$ (so $Y$ on
their intersection). The fold $S_\Psi$ is the stabilizer group generated by the $\varphi(a)$, $a$ a row of $H_X$; its
generators in our data are these, in the order of the X-checks.

This is the symplectic halving of [LOMK26, Def. 3] with the pairing $\pi=\Psi$; [LOMK26] credit the halving to
[BDB24], who show that a CSS code on $2m$ qubits with a fixed-point-free involutory ZX-duality is the double of a
code on $m$ qubits.

**Lemma 4.** In the situation of Definition 3, for a Pauli operator $p=(x\mid z)$ on $m$ qubits let
$\lambda(p)\in\mathbb F_2^{2m}$ be given by $\lambda(p)_{\sigma_i}=z_i$ and $\lambda(p)_{\Psi(\sigma_i)}=x_i$. Then
$\lambda$ is a linear bijection, the symplectic product satisfies $\langle\varphi(a),p\rangle=a\cdot\lambda(p)$,
and $\lambda(\varphi(a))=a\circ\Psi$. Consequently $p$ commutes with $S_\Psi$ iff $\lambda(p)\in\ker H_X$, and
$p\in S_\Psi$ iff $\lambda(p)\in\operatorname{rowsp}H_Z$; $\lambda$ maps the nontrivial logical operators of
$S_\Psi$ bijectively onto the nontrivial Z-logicals of $C$, $k(S_\Psi)=k(C)/2$, and $|\lambda(p)|=2y(p)+t(p)$,
where $y(p)$ is the number of $Y$ factors of $p$ and $t(p)$ the number of $X$ or $Z$ factors, while the Pauli weight
is $y(p)+t(p)$. In particular $d(S_\Psi)\ge\lceil d_Z(C)/2\rceil$.

*Proof.* $\langle\varphi(a),(x\mid z)\rangle=\sum_i(a_{\sigma_i}z_i+a_{\Psi(\sigma_i)}x_i)=a\cdot\lambda(p)$, as
every qubit is exactly one of $\sigma_i$, $\Psi(\sigma_i)$. By definition
$\lambda(\varphi(a))_{\sigma_i}=a_{\Psi(\sigma_i)}$ and $\lambda(\varphi(a))_{\Psi(\sigma_i)}=a_{\sigma_i}$, i.e.
$\lambda(\varphi(a))=a\circ\Psi$. The first equivalence follows as $S_\Psi$ is generated by the $\varphi(a)$ for the
rows $a$ of $H_X$. The permutation $\Psi$ maps $\operatorname{rowsp}H_X$ onto $\operatorname{rowsp}H_Z$, so
$\lambda(S_\Psi)=\{a\circ\Psi:a\in\operatorname{rowsp}H_X\}=\operatorname{rowsp}H_Z$, and the second equivalence
follows from the bijectivity of $\lambda$. Then $2k(S_\Psi)=\dim\ker H_X-\dim\operatorname{rowsp}H_Z=k(C)$; the
weights are read off the definition. $\square$

**Remark 5 (explicit generators; mirror codes).** For $\Psi=\Psi_c$ of Definition 2 the qubits $q<\Psi_c(q)$ are
those of the first block, so fold qubit $i$ is $(1,h)$ with $\iota(h)=i$, and the X-check $g$ gives the generator
$$S(g)=X^{\,g+F_A}\,Z^{\,c-g-F_B},\qquad g\in G,$$
with a $Y$ wherever $g+a=c-g-b$ for some $a\in F_A$, $b\in F_B$. Up to the order of the generators this is the
mirror code of [KL26] with the generators $Z(A+g)X(B-g)$, $g\in G$, for $A=c-F_B$ and $B=F_A$ (substitute $-g$ for
$g$). Translating the fold qubits by $u\in G$ turns the fold along $\Psi_c$ into the fold along $\Psi_{c+2u}$, so up
to qubit permutations only the class of $c$ in $G/2G$ matters.

**Lemma 6 (no single-qubit Clifford makes the code CSS).** Let $s_1,\dots,s_r$ generate a stabilizer group $S$ and
let $I$ be the set of indices $i$ such that (i) the only elements of $S$ supported in $\operatorname{supp}s_i$ are
$1$ and $s_i$ (up to phase). If a product $U$ of single-qubit Clifford operators maps $S$ onto a CSS group, then $U$
maps every $s_i$, $i\in I$, to a pure $X$- or pure $Z$-type operator.

*Proof.* Write $Us_iU^\dagger=ab$ with $a$ of $X$-type and $b$ of $Z$-type, both in $USU^\dagger$ (a CSS group is
the product of its $X$-type and its $Z$-type subgroup). Qubitwise $a$ is $I$ or $X$ and $b$ is $I$ or $Z$, so
$\operatorname{supp}(ab)=\operatorname{supp}a\cup\operatorname{supp}b$. Single-qubit Cliffords preserve supports,
so $U^\dagger aU$ and $U^\dagger bU$ are elements of $S$ supported in $\operatorname{supp}s_i$, hence in
$\{1,s_i\}$ by (i); as $ab=Us_iU^\dagger\ne1$, exactly one of $a$, $b$ is $Us_iU^\dagger$ and the other is $1$.
$\square$

Whether single-qubit Cliffords exist that make every $s_i$, $i\in I$, pure is a parity system with one unknown per
generator (equal letters on a common qubit force equal types, different letters different types, three letters on
one qubit are impossible); for $I$ the set of all generators it is the generator-level test of the board's verifier
(`stabilizer_code_is_not_locally_css`). For all five codes `src/rebuild.py` finds that (i) holds for every generator
and that the parity system has no solution; by Lemma 6 no product of single-qubit Cliffords maps the stabilizer
group onto a CSS group.

**Lemma 7 (meet in the middle).** Let $S$ be a stabilizer group on $n$ qubits and $a,b\ge0$. Every Pauli operator
$P$ that commutes with $S$ and has weight at most $a+b$ is a product $P=AB$ with $\operatorname{wt}A\le a$,
$\operatorname{wt}B\le b$ and equal syndromes.

*Proof.* Let $q_1<\dots<q_w$ be the support of $P$; let $A$ be the restriction of $P$ to $q_1,\dots,q_{\min(w,a)}$
and $B=PA^{-1}$ (up to phase), supported on the remaining at most $b$ qubits. The syndrome of $P$ is the sum of those
of $A$ and $B$ and vanishes. $\square$

`src/mitm.c` tabulates all $B$ of weight at most $b$ by syndrome, enumerates all $A$ of weight at most $a$, and tests
every product of equal syndrome for $AB\ne1$ and $AB\notin S$ (reduction against the reduced row echelon form of
$S$). If none is found, every nontrivial logical has weight at least $a+b+1$.

## The five codes

All five are folds along $\Psi_c$ (Definition 2) of the parents of Definition 1, with the generators of Remark 5;
`data/codes/<slug>.json` holds the generators (`checks.S`, exactly as in the board entry), the witness and the
construction data, `data/witnesses/<slug>.json` all the nontrivial logicals of weight $d$ that an
information-set search of ours found (2000 trials), the first of which is the witness of the entry.

| slug | parent, $c$ | $n$, rank, $k$ | generator weights | mirror data $(A,B)$ of [KL26] |
|---|---|---|---|---|
| `60-9-8` | `120-18-12`, $(1,0,0,0)$ | 60, 51, 9 | 48 of weight 8, 12 of weight 7 (one $Y$ each) | $A=\{(0,0,2,0),(1,0,0,0),(1,1,0,3),(1,1,1,1)\}$, $B=F_A$ |
| `56-6-8-s` | `112-12-8`, $(0,1)$ | 56, 50, 6 | 56 of weight 6, no $Y$ | $A=\{(0,1),(2,1),(20,1)\}$, $B=F_A$ |
| `32-9-6` | `64-18-8`, $(1,1,0)$ | 32, 23, 9 | 32 of weight 8, no $Y$ | $A=\{(0,3,3),(1,0,2),(1,1,0),(1,2,1)\}$, $B=F_A$ |
| `48-6-8-s` | `96-12-9`, $(0,1,0)$ | 48, 42, 6 | 32 of weight 7, 16 of weight 6 (one $Y$ each) | $A=\{(0,1,0),(0,2,0),(2,2,2),(2,3,1)\}$, $B=F_A$ |
| `33-2-7` | `66-4-10`, $(0,0)$ | 33, 31, 2 | 25 of weight 6, 7 of weight 5, 1 of weight 4 (9 $Y$ factors in 8 generators) | $A=\{(0,0),(1,4),(2,5)\}$, $B=F_A$ |

For `120-18-12` the classes of $c$ in $G/2G\cong(\mathbb Z/2)^2$ are represented by $(0,0,0,0)$, $(0,1,0,0)$,
$(1,0,0,0)$, $(1,1,0,0)$; all four folds are $[[60,9,8]]$ codes with $w=8$ (for each, the meet-in-the-middle
search of our verification excludes logicals of weight at most 7, and a witness has weight 8). The first two have
20 $Y$ factors (40 generators of weight 8, 20 of weight 7), the last two 12 (48 of weight 8, 12 of weight 7). The
four codes are pairwise inequivalent under qubit permutations and single-qubit Cliffords, which preserve the number
of nontrivial logicals of weight 8 (Pauli operators up to phase): in the order above it is 8, 4, 6 and 10 (an
exhaustive count by a meet-in-the-middle search with $a=b=4$ in our check of these submissions, a program other than
`src/mitm.c` whose logs are not part of this directory; our information-set search with 2000 trials had found the
same numbers of distinct logicals of weight 8). We submit $c=(1,0,0,0)$, the first of the two folds with 12 $Y$
factors in our enumeration; of these two it has fewer logicals of weight 8 (6 against 10). The fewest over all four
classes, 4, belong to $c=(0,1,0,0)$, which has 20 $Y$ factors.

**Distances (proved, computer-assisted).** Upper bounds: the witnesses. Lower bounds, each by three independent
exhaustive computations, one thread each:

| slug | $d$ | tropidist (exact, s) | `src/mitm.c`: no logical of weight $\le d-1$ ($a+b$, s) | positive control $a+b=d$ | Magma `MinimumDistance` (s) |
|---|---|---|---|---|---|
| `60-9-8` | 8 | 1.6 | $4+3$, 2.1 | not run | 8 (79.0) |
| `56-6-8-s` | 8 | 0.2 | $4+3$, 1.4 | not run | 8 (6.5) |
| `32-9-6` | 6 | 0.1 | $3+2$, 0.0 | weight 6 found | 6 (0.0) |
| `48-6-8-s` | 8 | 0.2 | $4+3$, 0.4 | weight 8 found | 8 (2.7) |
| `33-2-7` | 7 | 0.2 | $3+3$, 0.0 | weight 7 found | 7 (0.0) |

- tropidist (Cheng, Guo, Pan, [arXiv:2611.08247](https://arxiv.org/abs/2611.08247)), github.com/Fanerst/tropidist at
  `b03ec2d4077f`, the commit of the board's tropidist certificates, run as the board's
  `research/certs/tropidist_certify.py` runs it (`python -m cdist --matrices H.mtx L.mtx --threads 1 --time-limit
  600`) on the board's encoding for stabilizer codes ("pauli-mechanisms": three columns $X_i,Z_i,Y_i$ per qubit,
  the generators and a basis of the logicals as symplectic products, `src/tropidist_inputs.py`). Each run returned
  `exact` with lower bound equal to upper bound $d$, and its vector was rechecked by `src/tropidist_inputs.py check`
  (a nontrivial logical of weight $d$). The records are in `logs/tropidist/`.
- `src/mitm.c` (Lemma 7), written in our verification of the fold search, where it was validated against plain
  enumeration on 12 small folds and found a logical of weight exactly $d$ in 329 positive controls. Logs in
  `logs/mitm/` (`<slug>_a<a>_b<b>.log`); the runs here reproduce the statistics of the verification's runs line
  for line. For `60-9-8` and `56-6-8-s` the controls $a+b=8$ would need more than $2\cdot10^7$ table entries or
  $10^8$ queries (the limits of the verification) and were not run; there the witnesses give the upper bound.
- Magma V2.30-1, `MinimumDistance` of `QuantumCode(M : ExtendedFormat := true)` (`magma/<slug>.m`, logs in
  `logs/magma/`).

**Not CSS under single-qubit Cliffords:** Lemma 6 and `src/rebuild.py` (for `60-9-8`, `48-6-8-s`, `33-2-7` some
qubit already carries $X$, $Y$ and $Z$; for `56-6-8-s` and `32-9-6`, whose generators have no $Y$, the parity system
has an odd cycle).

**Literature (our check, November 2026).** The five codes are mirror codes of [KL26] (Remark 5) and none of them
is in the code table of [KL26]; that table has a CSS code with the parameters $[[56,6,8]]$ (over
$\mathbb Z_8\times\mathbb Z_7$, marked CSS there), which is not equivalent to `56-6-8-s` by Lemma 6. No published
code equivalent to one of the five was found in [KL26], [BDB24], [LOMK26] and their data, the genon codes of Hsin
(arXiv:2611.09538), or further non-CSS constructions of 2026. The parameters $[[56,6,8]]$ and $[[48,6,8]]$ are
known for CSS codes (the board's CSS entries `56-6-8` and `48-6-8`, the latter of weight 9; Panteleev and Kalachev,
arXiv:1904.02703, for $[[48,6,8]]$); for $[[60,9,8]]$, $[[32,9,6]]$ and $[[33,2,7]]$ we found no published code with
these parameters. The tables of codetables.de list stabilizer codes without weight restriction with larger
distances for all five $(n,k)$ ($[[60,9]]$: 12, $[[56,6]]$: 12, $[[32,9]]$: 7, $[[48,6]]$: 11, $[[33,2]]$: 10, lower
bounds, accessed 2026-11-10), so none of the parameters is a record outside the LDPC setting.

## Files

- `src/rebuild.py`: rebuilds and checks the five codes and their parents (standard library only);
  `src/tropidist_inputs.py` writes the tropidist inputs and rechecks a record; `src/same_results.py` compares
  reruns with our logs; `src/mitm.c` the meet-in-the-middle search (usage in its head); `src/gap_elements.g` GAP's
  numbering of the groups.
- `data/codes/<slug>.json`: the submitted codes (`checks.S`, a list of generators `{"X": [...], "Z": [...]}`, qubits
  from 0, as in the board entry), `witness`, `w`, `y_factors`, `qubit_degree`, `construction` (parent, group,
  $F_A$, $F_B$, $c$; `search_tag` is the name in our search: class `m0` to `m3` of the dualities $h\mapsto c-h$ of the
  parent).
- `data/parents/<parent>.json`: checks, witnesses and the label list of the five parents, copied from the board.
- `data/witnesses/<slug>.json`: the logicals of weight $d$ found by our information-set search.
- `instances/<slug>.txt`: the input of `src/mitm.c` (`n m`, then per generator the $X$ support and the $Z$
  support, each preceded by its size); `magma/<slug>.m`: the Magma computation.
- `logs/rebuild.log` (output of `make check`), `logs/environment.txt` (the system it ran on), `logs/mitm/`,
  `logs/magma/`, `logs/gap_elements.log`, `logs/tropidist/` (records written by the tool, wall times in `runs.tsv`,
  checksums of the inputs in `inputs.sha256`).

## Credits and licenses

The search for these folds, the programs, the computations and the texts of this directory were made with
Claude Opus 5.5 (Anthropic) in Claude Code. The construction is the halving of [BDB24] and [LOMK26]; the codes are
mirror codes in the sense of [KL26]. The parents are board entries by G. Lee (`120-18-12`, from arXiv:2610.13880),
@vprusso (`112-12-8`) and @MathysRennela (`64-18-8`, `96-12-9`, `66-4-10`). The files in `data/parents/` are derived
from files of github.com/unitaryfoundation/qldpc-challenge (commit `8042f9da`), which is licensed under the Apache
License 2.0; everything else is under the MIT license of this repository. tropidist is by Cheng, Guo and Pan
(arXiv:2611.08247).
