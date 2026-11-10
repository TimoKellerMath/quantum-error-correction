# A two-block code over C12 x Dic5 and two amplified codes for the qLDPC challenge board

Three CSS codes submitted to the [qLDPC challenge board](https://github.com/unitaryfoundation/qldpc-challenge)
(one pull request per code), with their check matrices, the proofs of their distances, programs and logs.

| board slug | code | construction | $w$ | qubit degree X / Z / total | $kd^2/n$ | distance |
|---|---|---|---|---|---|---|
| `480-16-26` | $[[480,16,26]]$ | two-block group algebra code over $C_{12}\times\mathrm{Dic}_5$ | 6 | 3 / 3 / 6 | 22.53 | exact (A) |
| `860-130-10` | $[[860,130,10]]$ | $[[4,2,2]]$-amplified `185-65-5` | 7 | 6 / 6 / 10 | 15.12 | exact (B) |
| `463-34-12` | $[[463,34,12]]$ | $[[4,2,2]]$-amplified `96-17-6` | 7 | 6 / 6 / 10 | 10.57 | exact (B) |

$w$ is the largest check weight. (A): Theorem 2, by two exhaustive search programs of ours (`src/exactA.c`,
`src/exactB.c`) and a parity argument. (B): Proposition 9 (the distance amplification doubles the distance exactly)
together with the exact base distances, which the same two programs certify; independently, tropidist and both
programs compute the distances of the amplified codes directly.

## Reproduce

Requirements: Python 3.8 or later (standard library only), a C compiler (gcc or clang), `make`, `sha256sum`.
Optional: Magma for `src/blocks.m`.

```
make check            # rebuild the three codes from their definitions and check them (about 2 s); log: logs/rebuild.log
make quick            # compile exactA, exactB; rerun the 32 searches of Part B and compare them with logs/ (seconds)
make long             # the second set of exact searches of Theorem 2 (about 11 CPU-minutes), compared with logs/
make long-first       # the first set (about 25 CPU-minutes), compared with logs/
make tropidist-inputs # the inputs of all our tropidist runs, compared with logs/tropidist/inputs.sha256
make blocks           # Magma: src/blocks.m (Proposition 5); our output is logs/480-16-26/blocks.log
```

`src/same_results.py` compares a rerun with our log: the `RESULT` line without its timing fields, every reported
word and the `CONCLUSION` line must agree (both programs are deterministic for a fixed instance, generator file and
number of threads). `python3 src/rebuild.py --board DIR`, with `DIR` a checkout of the board repository, also
compares `data/bases/` with the board's files and the three codes with their board entries once these exist. Our runs
used gcc 16.2.1 with `-O3 -march=native` on x86-64, one thread per job.

## Part A: the $[[480,16,26]]$ code

### Construction

$G=C_{12}\times\mathrm{Dic}_5=\langle z\rangle\times\langle x,y\mid x^5,y^4,yxy^{-1}=x^{-1}\rangle$ with $|z|=12$, the group
SmallGroup(240,40). Every element is $z^cx^iy^j$ with $c\in\mathbb Z/12$, $i\in\mathbb Z/5$, $j\in\mathbb Z/4$, and
$(z^cx^iy^j)(z^{c'}x^{i'}y^{j'})=z^{c+c'}x^{i+(-1)^ji'}y^{j+j'}$. The element $z^cx^iy^j$ has index $20c+4i+j$, and
qubit $(e,g)$ of block $e\in\{0,1\}$ has index $240e+\mathrm{index}(g)$. With
$$a=1+z^8x^3y^2+z^{11}y^3,\qquad b=1+z^2xy^2+y,$$
the code is $C(G;a,b)$: $H_X=(L_a\mid R_b)$, $H_Z=(R_b^{\mathsf T}\mid L_a^{\mathsf T})$, where $L_a[g,h]=1$ iff
$g=sh$ for some $s\in a$ and $R_b[g,h]=1$ iff $g=ht$ for some $t\in b$. Explicitly, the X-check at $c\in G$ consists
of the block-0 qubits $s^{-1}c$ ($s\in a$) and the block-1 qubits $ct^{-1}$ ($t\in b$), the Z-check at $c$ of the
block-0 qubits $ct$ ($t\in b$) and the block-1 qubits $sc$ ($s\in a$), and the checks are ordered by the index of $c$.
This is the convention and the indexing of the board entries `720-16-36` and `960-24-38` (both over
$C_m\times\mathrm{Dic}_5$), and the convention of our $[[600,24,26]]$ code in `weight6-exact-distances/`. Identifying $\mathbb F_2^G$ with $R=\mathbb F_2[G]$,
$H_X(u,v)=au+vb$ and $H_Z^{\mathsf T}w=(wb,aw)$, so $H_XH_Z^{\mathsf T}=0$ and the logical operators are the middle
homology of $R\to R^2\to R$. These are the two-block group algebra codes of Lin and Pryadko (arXiv:2306.16400).

Computed (`src/rebuild.py`; again with Magma's group arithmetic in `src/blocks.m`): $\operatorname{rank}H_X=
\operatorname{rank}H_Z=232$, so $k=16$; every check has weight 6 (three qubits in each block), every qubit lies in
three X-checks and three Z-checks, and the Tanner graph is connected.

**Lemma 1 (parity).** The sum of the rows of $H_X$ and the sum of the rows of $H_Z$ are the all-one vector
$\mathbf 1$. Hence every element of $\ker H_X$ and of $\ker H_Z$ has even weight, and $d_X$, $d_Z$ are even.

*Proof.* Every column of $H_X$ has weight 3, so the sum of the rows of $H_X$ is $\mathbf 1$; for $u\in\ker H_X$,
$|u|\equiv\langle u,\mathbf 1\rangle=\sum_r\langle u,r\rangle=0\pmod 2$, the sum over the rows $r$. Likewise for
$H_Z$. $\square$

### Distance

**Theorem 2.** $d_X=d_Z=26$.

*Proof.* Upper bound: the witnesses in `data/codes/480-16-26.json` (and in the board entry) are nontrivial logical
operators of weight 26 on both sides (checked by `src/rebuild.py`: in the kernel of the other check matrix, outside
the row space of their own). Lower bound: by Lemma 1 it suffices to exclude nontrivial logicals of weight at most 24.
`src/exactA.c` and `src/exactB.c` are exhaustive searches: given $H_X$, $H_Z$, a side, a bound $W$ and generators of
a group of qubit permutations that maps the rows of $H_X$ onto themselves and the rows of $H_Z$ onto themselves, a run
that reports no word proves that the side has no nontrivial logical of weight at most $W$ (if $W$ is odd and the
parity of Lemma 1 holds, the program searches only the weights up to $W-1$, which suffices; each program tests the
parity by a rank computation). The completeness proofs are in
`src/exactA_NOTES.md` (Lemma A, Corollary B and the lemmas and remarks after them) and `src/exactB_NOTES.md`
(Lemmas N1 to N6), on top of
Proposition 3.1 and Corollary 3.2 of `weight6-exact-distances/exact_distances.pdf` in this repository, whose two
implementations (`vsearch_orb.c`, `cluster2.c`) the programs generalize to any check weight, any $k$ and any group.
Both programs verify that the given generators preserve the rows, and `src/rebuild.py` checks it again. The eight
runs below report no word. $\square$

| set | instance | group | program, decomposition | side | $W$ (searched) | nodes | CPU s |
|---|---|---|---|---|---|---|---|
| first | `c480_16` | order 240, 4 generators (nauty) | exactA, one anchor per qubit orbit, earlier orbits excluded | Z | 25 (24) | 106,226,591 | 313.6 |
| | | | | X | 25 (24) | 105,374,138 | 365.6 |
| | | | exactB, the same anchors | Z | 25 (24) | 110,404,866 | 349.0 |
| | | | | X | 25 (24) | 108,721,700 | 402.2 |
| second | `c480v` | order 240, 3 generators (Magma) | exactA, pattern classes of the single check orbit (31 jobs) | Z | 24 (24) | 30,042,832 | 49.9 |
| | | | | X | 24 (24) | 28,171,810 | 47.2 |
| | | | exactB, one anchor per qubit orbit, orbits in the other order | Z | 24 (24) | 122,361,101 | 268.7 |
| | | | | X | 24 (24) | 119,022,416 | 261.4 |

The instances: `instances/c480_16.*` is the code in another labelling ($z^cx^i(y^2)^{j_2}y^e$ at
$120e+2(5c+i)+j_2$ in each block, the function `rep_index` of `src/rebuild.py`); `instances/c480v.*` is the entry's
labelling relabelled by a random permutation $q\mapsto p[q]$ (`c480v.perm`, seed 4802616) with the checks and the
qubits within each check shuffled. The groups are automorphism groups of the Tanner graph with the qubits, the
X-checks and the Z-checks as fixed vertex classes, computed independently with nauty's dreadnaut and with Magma; both
have order 240 and two qubit orbits of size 240, and the first orbit of `c480v` is the second of `c480_16`, so the
anchors and exclusions of the two sets differ (`logs/480-16-26/second/gens_check.log`). Logs:
`logs/480-16-26/first/`, `logs/480-16-26/second/` (also two calibration runs of the pattern mode at $W=18,20$ and the
Magma construction `cross.log`).

Further evidence, not used in the proof:
- tropidist (Cheng, Guo, Pan, arXiv:2611.08247; github.com/Fanerst/tropidist at commit `b03ec2d4077f`, the commit of
  the board's tropidist certificates), one thread, 1200 s per side, on the board's css-sides problems
  (`src/tropidist_inputs.py`; side X: $H_Z$ against a basis of the Z-logicals, side Z mirrored): status timeout on
  both sides, with the proved lower bounds $d_Z\ge22$ and $d_X\ge22$ (all thresholds up to 21 exhausted; the tool
  used no symmetry and no parity). Records: `logs/480-16-26/tropidist/c480.{X,Z}.json`. The time per threshold grew by
  a factor of about 2.8 per unit, so reaching 26 this way would take of the order of 15 to 20 CPU-hours per side
  (an extrapolation).
- Information-set search (our program, not included), $2\cdot10^4$ trials per side (seed 21): nothing lighter than
  26, and 31 resp. 16 distinct words of weight 26 on the sides Z and X (`logs/480-16-26/first/witness_search.log`).

### Where the 16 logical qubits sit

$A=\langle z,x,y^2\rangle\cong C_{12}\times C_5\times C_2$ is an abelian normal subgroup of index 2. Its 2-part
$A_2=\langle z^3,y^2\rangle\cong C_4\times C_2$ is central in $G$, its odd part $A_{2'}=\langle t\rangle$,
$t=z^4x$, is cyclic of order 15, and conjugation by $y$ is $\tau$: $z\mapsto z$, $x\mapsto x^{-1}$, $y^2\mapsto y^2$,
so $yty^{-1}=t^4$. The primitive idempotents of $\mathbb F_2[\langle t\rangle]=\mathbb F_2[t]/(t^{15}+1)$ belong to
the factors $t+1$, $t^2+t+1$, $t^4+t+1$, $t^4+t^3+1$, $t^4+t^3+t^2+t+1$; they are fixed by $t\mapsto t^2$, hence by
$t\mapsto t^4$, so they are central in $\mathbb F_2[G]$ and split the code into five summands. Over
$K=\mathbb F_{16}$ a character $\chi$ of $A_{2'}$ of order 15 has the $\tau$-orbit $\{\chi,\chi^\tau\}$ with
$\chi^\tau=\chi^4$, inside its Galois orbit of size 4; the block of $\mathbb F_2[G]$ of such an orbit is
$M_2(\mathbb F_4[A_2])$, of dimension 64. Write $a=\alpha+y\beta$ and $b=\gamma+\delta y$ with
$\alpha,\beta,\gamma,\delta\in\mathbb F_2[A]$ (here $\alpha=1+z^8x^3y^2$, $\beta=z^{11}y^2$, $\gamma=1+z^2xy^2$,
$\delta=1$), extend $\chi$ to the ring map $\mathbb F_2[A]\to\mathcal O=K[A_2]$,
$z^cx^iy^{2j}\mapsto\chi(z^{4c}x^i)\,z^{9c}y^{2j}$ ($z^{4c}$ and $z^{9c}$ are the odd and the 2-part of $z^c$), and
write $u(\chi)$ for the image of $u$. In the block of $\chi$, left multiplication by $a$ and $b$ on
$K[G]e_\chi=\mathcal Oe_\chi\oplus\mathcal Oye_\chi$ is given by the Clifford symbols
$$P_a=\begin{pmatrix}\alpha(\chi)&y^2\beta(\chi^\tau)\\\beta(\chi)&\alpha(\chi^\tau)\end{pmatrix},\qquad
P_b=\begin{pmatrix}\gamma(\chi)&y^2\delta(\chi)\\\delta(\chi^\tau)&\gamma(\chi^\tau)\end{pmatrix},$$
whose determinants are $N(a)(\chi)$ and $N(b)(\chi)$ for the norms $N(a)=\alpha\tau(\alpha)+y^2\beta\tau(\beta)$,
$N(b)=\gamma\tau(\gamma)+y^2\delta\tau(\delta)$ in $\mathbb F_2[A]$.

*Norm formula.* $k=h_0+h_2$ with $h_0=\dim R/(aR+Rb)$ and $h_2=\dim R/(Ra+bR)$ (G. Lee, arXiv:2610.13880,
Proposition B.1), and both split over the central idempotents. If the residue matrices of $P_a$ and $P_b$ modulo the
maximal ideal $\mathfrak m$ of the local ring $\mathcal O$ both have rank at least 1, then the $\tau$-orbit
$\{\chi,\chi^\tau\}$ contributes $\dim_K\mathcal O/(N(a)(\chi),N(b)(\chi))$ to $\dim_K$ of each of $h_0$, $h_2$. Indeed,
over a local ring a $2\times2$ matrix with a unit entry is equivalent to $\mathrm{diag}(1,\det)$; for $P_a=UDV$,
$P_b=U'D'V'$ the bijection $X\mapsto U^{-1}XV'^{-1}$ maps $P_aM_2(\mathcal O)+M_2(\mathcal O)P_b$ onto
$DM_2(\mathcal O)+M_2(\mathcal O)D'$, whose quotient is $\mathcal O/(N(a)(\chi),N(b)(\chi))$ in the entry $(2,2)$ and
zero elsewhere; the same holds for $h_2$. Summed over the Galois orbit this gives
$2\,[\mathbb F_4:\mathbb F_2]\dim_K\mathcal O/(N(a)(\chi),N(b)(\chi))$ for $k$.

**Proposition 3 (computed by `src/rebuild.py`).**
1. All 16 logical qubits lie in the summand of the factor $t^4+t+1$, that is at the Galois orbit $O_1$ of the
   characters $\chi$ of order 15 with $\chi(t)$ a root of $t^4+t+1$; the other four summands have $k=0$ (computed
   by ranks on $e\,\mathbb F_2[G]$ for each idempotent $e$).
2. For $\chi(z^4)=g^5$, $\chi(x)=g^3$ ($g^4=g+1$ in $\mathbb F_{16}$, so $\chi(t)=g^8$ is a root of $t^4+t+1$), the
   residue matrices of $P_a$ and $P_b$ have rank 1, and $\mathcal O/(N(a)(\chi),N(b)(\chi))$ has length 4 and Hilbert
   function $(1,1,1,1)$, so it is $K[\varepsilon]/(\varepsilon^4)$: the two norms cut out a curvilinear point of
   length 4 in the fat point $\operatorname{Spec}K[C_4\times C_2]$ of length 8 and embedding dimension 2. By the
   norm formula, $O_1$ contributes $2\cdot2\cdot4=16$, in agreement with item 1.

**Proposition 4 (computed by `src/rebuild.py`).** The two-block codes with the same supports over
$C_m\times\mathrm{Dic}_5$ ($z$-exponents read modulo $m$; for $m=24$ the same exponents, a double cover along
$z^{12}$) have $k=4,8,16,16$ for $m=3,6,12,24$ ($n=120,240,480,960$). Over $G/\langle z^4\rangle\cong
C_4\times\mathrm{Dic}_5$ and $G/\langle x\rangle\cong C_{12}\times C_4$ they have $k=0$, over
$G/\langle y^2\rangle\cong C_{12}\times D_5$ they have $k=8$. So along the 2-tower $k=4\ell$ with
$\ell=1,2,4,4$: if only the orbit $O_1$ carries in each case (we checked this for $m=12$ only), $\ell$ is the
length of the contact scheme, which grows with the 2-part of $\langle z\rangle$ up to 4 and then stays at 4.

**Proposition 5 (Brauer sandwich).** Let $G$ be a finite group, $K$ a splitting field of characteristic 2,
$S_1,\dots,S_r$ the simple $K[G]$-modules, $(c_{ij})$ the Cartan matrix, $a,b\in\mathbb F_2[G]$, and
$\kappa_i(u)=\dim_K\ker(u|_{S_i})$. Then $k_{ss}\le k\le k_C$ for $C(G;a,b)$, with
$k_{ss}=2\sum_i\kappa_i(a)\kappa_i(b)$ and $k_C=2\sum_{i,j}c_{ij}\kappa_i(a)\kappa_j(b)$; more precisely
$k_{ss}/2\le h_0,h_2\le k_C/2$.

*Proof.* Let $R=K[G]$ with radical $J$. Lower bound: $R/(aR+Rb)$ maps onto
$\prod_i\operatorname{End}(S_i)/(a\operatorname{End}(S_i)+\operatorname{End}(S_i)b)$, and for linear maps $A,B$ of an
$n$-dimensional space, $\dim\operatorname{End}/(A\operatorname{End}+\operatorname{End}B)=(n-\operatorname{rk}A)
(n-\operatorname{rk}B)$. Upper bound: $X\mapsto X/(aX+Xb)$ is right exact on $R$-bimodules, so its dimension is
subadditive along a bimodule composition series of $R$; the factors are $\operatorname{Hom}_K(S_j,S_i)$ with
multiplicity $\dim e_iRe_j=c_{ij}$, and each contributes $\kappa_i(a)\kappa_j(b)$. For $h_2$ use $X/(Xa+bX)$ and
$c_{ij}=c_{ji}$. $\square$

For this code `src/blocks.m` (Magma; `logs/480-16-26/blocks.log`) finds 9 simple modules over $\mathbb F_{16}$, of
which exactly two (of dimension 2, the two $\tau$-orbits of $O_1$) have $\kappa(a)=\kappa(b)=1$ (the rank-1 residues
of Proposition 3), with Cartan entries $c_{ii}=8=|A_2|$; so $k_{ss}=4$, $k_C=32$, and $k=16$ ($h_0=h_2=8$). The
semisimple quotient accounts for 4 of the 16 logical qubits; the other 12 come from the central defect group $A_2$,
through the contact of length 4 inside the local ring of length 8.

### How the code was found

By a design from the structure above (our "norm design"), not by a uniform search. Write $a=1+u+yv$ and $b=1+u'+v'y$
with $u,u',v,v'\in A$ (every weight-3 element with two terms in one coset of $A$ has this form up to two-sided
translations, which give equivalent codes). At a character of order 15 the residue of $N(a)(\chi)$ vanishes for
exactly one value of the $C_3$-component of $v$ per value of the odd part of $u\ne1$, which lists the pairs singular
at one Galois orbit of order-15 characters; we enumerated all such normal-form pairs on $C_{12}\times\mathrm{Dic}_5$
by the ideal $(N(a)(\chi),N(b)(\chi))$ (78,848 pairs with contact length 4 at the orbit $O_1$ and nothing at the
other orbits, hence $k=16$; as many at $O_2$), drew 320 of them uniformly, and screened each with an information-set
search ($10^3$ trials per side and searches restricted to words invariant under central subgroups). Ten of the 320
pairs reached a screening value of at least 26 (one of them 28, lowered to 26 by a deeper search), and these ten are
all equivalent to this code (canonical forms of the Tanner graphs); none of 80 uniform weight-3 pairs with $k=16$ on
the same group exceeded 20. Pairs with non-curvilinear contact (Hilbert function $(1,2,1)$) gave at most 12 in the
same screening.

### Novelty and related searches

At board commit `8042f9da` no entry has $(n,k)=(480,16)$; no CSS entry dominates or equals the code in
$(n,k,d,w)$, and it dominates 50 entries, among them `504-16-26` and `510-16-26` on $n$ alone. The 61 entries with
$n\le480$, $k\ge16$ and $d\ge26$ all have check weight at least 8.

The board contributor @MicheldeMessieres reports in the board's `notes/720-16-36.md` (entries `720-16-36` and
`960-24-38`) an exhaustive census, up to equivalence, of the weight-6 (3+3) two-block codes over the groups
$C_m\times H$ with $H$ among a list of small non-abelian groups including $\mathrm{Dic}_5$, for group orders up to 480,
in the same convention and indexing. $C_{12}\times\mathrm{Dic}_5$ has order 240, so that census very likely contains
this code; the code is not submitted or described there (the note names near misses at $n=720$ and $n=800$). We
found it independently, by the design above (a single carrying orbit of order 15, a curvilinear point
$\operatorname{Spec}\mathbb F_4[\varepsilon]/(\varepsilon^4)$); our contribution is its identification and the exact
distance.

Literature checked for $[[480,16,d\ge26]]$ codes and for codes over $C_{12}\times\mathrm{Dic}_5$, with no hit:
H.-K. Lin and L. P. Pryadko, arXiv:2306.16400 (enumeration up to $n=200$), and R. Wang, H.-K. Lin, L. P. Pryadko,
arXiv:2305.06890 (abstract only); G. Lee, arXiv:2610.13880 (tables up to $n=320$); A. Aydin, I. Tamo, A. Barg,
arXiv:2606.17268 (tables read up to $n=306$, without the appendix table 6); R. P. Rajpurohit, V. Lalitha, S. S. Garani, arXiv:2610.12782 ($120\le n\le400$);
arXiv:2606.24808, arXiv:2608.08996, arXiv:2608.12509 (no $n=480$ code with $k=16$); the Error Correction Zoo. QECDB
was not checked.

## Part B: amplified codes

### Construction and distance

Notation: a CSS code is a pair $(H_X,H_Z)$ of $\mathbb F_2$-matrices with $n$ columns and $H_XH_Z^{\mathsf T}=0$,
whose rows need not be independent; $k=n-\operatorname{rk}H_X-\operatorname{rk}H_Z$,
$d_Z=\min\{|u|:u\in\ker H_X\setminus\operatorname{rowsp}H_Z\}$,
$d_X=\min\{|v|:v\in\ker H_Z\setminus\operatorname{rowsp}H_X\}$. For a qubit $q$, $\deg_X(q)$ and $\deg_Z(q)$ are the
numbers of rows of $H_X$, $H_Z$ containing $q$; $e_a$ are the unit vectors of $\mathbb F_2^N$ and $\mathbf 1$ the
all-one vector.

**Definition 6.** Let $(H_X,H_Z)$ be a CSS code on $n$ qubits with rows $H_X[0],\dots,H_X[m_X-1]$ and
$H_Z[0],\dots,H_Z[m_Z-1]$, and let $N\ge2$ be even. The code $\mathrm{Amp}_N(H_X,H_Z)$ has the qubits
- central $(q,a)$, $0\le q<n$, $0\le a<N$, with index $Nq+a$;
- $t_j$ for each Z-check $j$, index $Nn+j$;
- $s_i$ for each X-check $i$, index $Nn+m_Z+i$;

and the checks (in this order)
- X-checks $(c_i,a)=\{(q,a):q\in H_X[i]\}\cup\{s_i\}$ ($i$, then $a$) and
  $(q)^X=\{(q,a):0\le a<N\}\cup\{t_j:q\in H_Z[j]\}$;
- Z-checks $(z_j,a)=\{(q,a):q\in H_Z[j]\}\cup\{t_j\}$ and $(q)^Z=\{(q,a):0\le a<N\}\cup\{s_i:q\in H_X[i]\}$.

This is the code of Liang, Gu, Chen, Eisert and Wang (arXiv:2609.37231, cited as [LGCE26], Eqs. (14), (15)) for
the amplifier $[[N,N-2,2]]$ with $G_X=G_Z=\mathbf 1^{\mathsf T}$: the truncation $T_3\to T_2\to T_1$ of the tensor
product of the base complex $Q_2\xrightarrow{H_Z^{\mathsf T}}Q_1\xrightarrow{H_X}Q_0$ with
$\mathcal A_N\colon\mathbb F_2\xrightarrow{\mathbf 1}\mathbb F_2^N\xrightarrow{\mathbf 1^{\mathsf T}}\mathbb F_2$. $N=4$
is the $[[4,2,2]]$ amplifier, used for both codes here.

**Lemma 7.** In the situation of Definition 6 with $N\ge4$:
1. $\mathrm{Amp}_N(H_X,H_Z)$ is a CSS code with $n'=Nn+m_X+m_Z$ and $k'=(N-2)k$. If the rows of $H_X$ and of $H_Z$
   are linearly independent, $n'=(N+1)n-k$.
2. $(c_i,a)$ has weight $|H_X[i]|+1$, $(q)^X$ weight $N+\deg_Z(q)$, $(z_j,a)$ weight $|H_Z[j]|+1$, $(q)^Z$ weight
   $N+\deg_X(q)$, so the largest check weight is $w'=\max(w_X+1,w_Z+1,N+q_Z,N+q_X)$, with $w_X,w_Z$ the largest row
   weights and $q_X,q_Z$ the largest column weights of $H_X,H_Z$.
3. The qubit degrees (X, Z) are $(\deg_X(q)+1,\deg_Z(q)+1)$ for $(q,a)$, $(N,|H_X[i]|)$ for $s_i$ and
   $(|H_Z[j]|,N)$ for $t_j$.

*Proof.* Orthogonality, by cases: $(c_i,a)$ and $(z_j,a')$ meet in $|H_X[i]\cap H_Z[j]|$ central qubits if $a=a'$
and in none otherwise, an even number; $(c_i,a)$ and $(q)^Z$ share $(q,a)$ and $s_i$ if $q\in H_X[i]$ and nothing
otherwise; $(q)^X$ and $(z_j,a)$ likewise with $t_j$; $(q)^X$ and $(q')^Z$ meet in $N$ central qubits if $q=q'$ and
in none otherwise. $\ker\tilde H_X/\operatorname{rowsp}\tilde H_Z$ is the homology of the truncation at $T_2$, which
equals that of the full tensor product (it depends only on the two differentials at $T_2$). By the Künneth formula
over a field it is $\bigoplus_{i+j=2}H_i(\mathcal Q)\otimes H_j(\mathcal A_N)$, induced by
$[u]\otimes[\ell]\mapsto[u\otimes\ell]$ ($u\otimes\ell$ placed on the central qubits). Here
$H_0(\mathcal A_N)=H_2(\mathcal A_N)=0$ and $H_1(\mathcal A_N)=\mathbf 1^\perp/\langle\mathbf 1\rangle$ has dimension
$N-2$ because $\mathbf 1\in\mathbf 1^\perp$ for even $N$. So $k'=(N-2)k$; the rest is counting, and with independent
rows $m_X+m_Z=n-k$. $\square$

**Lemma 8 (contraction).** For $\omega\in\mathbb F_2^N$ of even weight and $z\in\mathbb F_2^{n'}$ let $z_a$ be the
restriction of $z$ to the qubits $(\cdot,a)$ and $c_\omega(z)=\sum_{a:\omega_a=1}z_a$. Then
$c_\omega(\ker\tilde H_X)\subseteq\ker H_X$, $c_\omega(\operatorname{rowsp}\tilde H_Z)\subseteq\operatorname{rowsp}H_Z$
and $c_\omega(u\otimes\ell)=\langle\omega,\ell\rangle u$.

*Proof.* If $\tilde H_Xz=0$, the check $(c_i,a)$ gives $(H_Xz_a)_i=z(s_i)$; summing over the $|\omega|$ values $a$
with $\omega_a=1$ gives $(H_Xc_\omega(z))_i=|\omega|z(s_i)=0$. On generators of $\operatorname{rowsp}\tilde H_Z$:
$c_\omega((z_j,a))$ is $H_Z[j]$ or 0, and $c_\omega((q)^Z)=|\omega|e_q=0$. The last claim is the definition. $\square$

**Proposition 9 (exact doubling; [LGCE26, Lemma 4]).** Let $(H_X,H_Z)$ be a CSS code with $k\ge1$ and let $N\ge4$
be even. Then $d_Z(\mathrm{Amp}_N(H_X,H_Z))=2d_Z$ and $d_X(\mathrm{Amp}_N(H_X,H_Z))=2d_X$.

*Proof.* Lower bound. Let $z\in\ker\tilde H_X\setminus\operatorname{rowsp}\tilde H_Z$ and let
$[\ell_1],\dots,[\ell_{N-2}]$ be a basis of $\mathbf 1^\perp/\langle\mathbf 1\rangle$. By the isomorphism in the proof
of Lemma 7, $z-\sum_ju_j\otimes\ell_j\in\operatorname{rowsp}\tilde H_Z$ for some $u_j\in\ker H_X$ with
$u_i\notin\operatorname{rowsp}H_Z$ for some $i$. The standard form induces a perfect pairing on
$\mathbf 1^\perp/\langle\mathbf 1\rangle$ (the radical of $\mathbf 1^\perp$ is $\langle\mathbf 1\rangle$), so there is
$\omega\in\mathbf 1^\perp$ with $\langle\omega,\ell_j\rangle=\delta_{ij}$; then $\omega'=\omega+\mathbf 1$ has even
weight and $\langle\omega',\ell_j\rangle=\delta_{ij}$ as well, since every $\ell_j$ has even weight. By Lemma 8,
$c_\omega(z)$ and $c_{\omega'}(z)$ are congruent to $u_i$ modulo $\operatorname{rowsp}H_Z$, so both have weight at
least $d_Z$, and as the supports of $\omega$ and $\omega'$ partition $\{0,\dots,N-1\}$,
$2d_Z\le|c_\omega(z)|+|c_{\omega'}(z)|\le\sum_a|z_a|\le|z|$. Upper bound. For a Z-logical $u$ of the base of weight
$d_Z$, $u\otimes(e_0+e_1)$ is a nontrivial logical of weight $2d_Z$ by Lemma 7, as
$e_0+e_1\in\mathbf 1^\perp\setminus\langle\mathbf 1\rangle$ for $N\ge4$. The X-statement follows by exchanging the
roles of $H_X$ and $H_Z$ (and of $s$ and $t$), which maps Definition 6 to itself. $\square$

[LGCE26, Lemma 4] gives $\lceil\alpha_Pd_P\rceil\le\tilde d_P\le d_Pd_P^A$ for amplifiers whose logical classes have
families of representatives with bounded overlap; the families $\{\omega,\omega+\mathbf 1\}$ give $\alpha_P=2$ and
$d^A=2$. The version above writes the proof out in the coordinates of our data files, so that its hypotheses reduce to
what `src/rebuild.py` checks: the code is $\mathrm{Amp}_4$ of the base entry, and $k\ge1$.

### The two codes

| slug | base entry (author on the board) | base rows | base $[[n,k,d]]$, $d_X$ / $d_Z$ | $n'=5n-k$, $k'=2k$ | $w'$ |
|---|---|---|---|---|---|
| `860-130-10` | `185-65-5` (@mathysrennela; a puncturing transform of `186-64-6`) | all 60 / 60, independent | $[[185,65,5]]$, 5 / 5 | 860, 130 | 7 |
| `463-34-12` | `96-17-6` (@vprusso; `96-16-6` with one X-check deleted) | all 39 / 40, independent | $[[96,17,6]]$, 6 / 6 | 463, 34 | 7 |

The base check matrices and witnesses in `data/bases/` are copied from the board (Apache License 2.0, see below);
`src/rebuild.py` rebuilds both amplified codes from them by Definition 6 and compares row by row, and
`src/rebuild.py --board DIR` checks that the copies coincide with the board's files. The witnesses of the entries
are $u\otimes(e_0+e_1)$ for the witnesses $u$ of the bases.

**Base distances.** At board commit `8042f9da` the two bases carry no certificate (their distances are claimed as
upper bounds; certificates made with the same two programs are part of `board-certificates-2026-11/` of this
repository and of the board's pull request #5306). exactA and exactB certify them exactly on both sides (Tanner-graph groups trivial by nauty and by
Magma, so one anchor per qubit; `logs/amplified/cert/`): for each side, the run at $W=d-1$ reports no word and the run
at $W=d$ finds one.

| base | side | exactA: $W=d-1$, nodes | exactA: $W=d$ | exactB: $W=d-1$, nodes | exactB: $W=d$ | result |
|---|---|---|---|---|---|---|
| `185-65-5` | X | 4: none, 525 | word of weight 5 | 4: none, 638 | word of weight 5 | $d_X=5$ |
| `185-65-5` | Z | 4: none, 517 | word of weight 5 | 4: none, 577 | word of weight 5 | $d_Z=5$ |
| `96-17-6` | X | 5 (parity: 4): none, 95 | word of weight 6 | 5 (parity: 4): none, 104 | word of weight 6 | $d_X=6$ |
| `96-17-6` | Z | 5: none, 389 | word of weight 6 | 5: none, 362 | word of weight 6 | $d_Z=6$ |

With Proposition 9: $d_X=d_Z=10$ for `860-130-10` and $d_X=d_Z=12$ for `463-34-12`.

**Direct computations (independent of Proposition 9).**
- tropidist at commit `b03ec2d4077f`, one thread, time limit 600 s, the board's css-sides problems
  (`src/tropidist_inputs.py`): status exact on all four sides, $d_X=d_Z=10$ for `860-130-10` (2.9 s and 1.3 s) and
  $d_X=d_Z=12$ for `463-34-12` (0.6 s and 0.6 s); the vectors it returned were rechecked with
  `src/tropidist_inputs.py check` (nontrivial logicals of the stated weight). Records in `logs/tropidist/`.
- exactA and exactB on the amplified codes, with their Tanner-graph groups of order 24 (the permutations of the four
  copies; computed with nauty and with Magma): nothing at $W=9$ resp. $W=11$ on either side (30,573 to 33,543 resp.
  16,857 to 21,387 nodes per run), words of weight 10 resp. 12 at the next $W$ (`logs/amplified/amp/`).

**Context on the board** (commit `8042f9da`). Both codes have check weight 7 and are on the frontier of the
weight-8 cell only through it: entries with check weight 8 are better in $n$, $k$ and $d$ (for example `512-132-14`
against `860-130-10`, and `300-36-13` against `463-34-12`). Among the entries with check weight at most 7, the
largest $d$ with $n\le860$, $k\ge130$ is 8 (`675-139-8` and others), and with $n\le463$, $k\ge34$ it is 10
(`382-36-10`, `384-36-10`). `463-34-12` dominates `464-32-12` (same $d$ and check weight, better $n$ and $k$), the
$[[4,2,2]]$ amplification of `96-16-6`; `96-17-6` is `96-16-6` with one X-check deleted, which raises $k$ by one and
lowers $n'=5n-k$ by one.

## Files

- `src/rebuild.py`: rebuilds and checks the three codes, the instances and generators of all searches, the
  tropidist problems, Propositions 3 and 4 (standard library only).
- `src/exactA.c`, `src/exactB.c`: the two exhaustive search programs (usage in `src/exactA_NOTES.md`,
  `src/exactB_NOTES.md`, which also contain the completeness proofs and the validation of the programs; paths named
  there refer to our working tree, not to this directory). These four files are identical to those in
  `board-certificates-2026-11/src/` of this repository. `src/same_results.py`: comparison of reruns with our logs.
  `src/tropidist_inputs.py`: the tropidist problems and the recheck of tropidist records. `src/blocks.m`: Magma,
  Proposition 5.
- `data/codes/<slug>.json`: the submitted codes: `checks` exactly as in the board entry (`X`, `Z`: lists of
  supports, qubits from 0), `witness`, `construction`. `data/bases/<base>.json`: checks and witnesses of the two base
  entries, copied from the board.
- `instances/`: inputs of the searches (format: `n m`, then one check per line): `c480.{hx,hz}` the entry's
  labelling; `c480_16.*` and `c480v.*` the two instances of Theorem 2 with their generator files (`.gens`, `.gensA` in
  exactA's format, `.gensB` in exactB's), witnesses and the permutation `c480v.perm` (`c480_16.code.json` records the
  labelling of `c480_16`); `amplified/A/`, `amplified/B/` the relabelled bases with their permutations
  (`q -> p[q]`) and groups, `amplified/amp422_<base>.*` the amplified codes in the entries' labelling.
- `logs/480-16-26/`: `first/`, `second/` the runs of Theorem 2 and the witness search, `tropidist/` the two tropidist
  records, `blocks.log` the output of `src/blocks.m`. `logs/amplified/`: `cert/` the runs on the bases, `amp/` the
  runs on the amplified codes. `logs/tropidist/`: the tropidist records of the amplified codes, their wall times
  (`runs.tsv`), and the checksums of the inputs of all six tropidist runs (`inputs.sha256`; the runs on `480-16-26`
  were named `c480.X`, `c480.Z`). `logs/rebuild.log`: output of `make check`; `logs/environment.txt`: the system it
  ran on.

The searches were run in working copies with other paths: in the logs of `logs/amplified/`, `inst/A/`, `inst/B/` and
`cand/` are `instances/amplified/A/`, `instances/amplified/B/` and `instances/amplified/`; the first runs of Theorem 2
used the instance files without a directory, the second ones `inst/`. In `logs/480-16-26/first/exactB_{X,Z}.log` the
absolute path of the program in the line `cmd:` was shortened to `bin/exactB`; otherwise the logs are unchanged.

## Credits and licenses

The design and identification of `480-16-26`, the two amplified codes, the programs, the computations and the texts
of this directory were made with Claude Opus 5.5 (Anthropic) in Claude Code. The two-block construction is that of
Lin and Pryadko (arXiv:2306.16400), the amplification that of [LGCE26]; the base codes are board entries by
@mathysrennela (`185-65-5`) and @vprusso (`96-17-6`). The files in `data/bases/` are derived from files of
github.com/unitaryfoundation/qldpc-challenge (commit `8042f9da`), which is licensed under the Apache License 2.0;
everything else is under the MIT license of this repository. tropidist is by Cheng, Guo and Pan (arXiv:2611.08247);
nauty is by McKay and Piperno.

## tropidist certificate of [[480,16,26]]

`logs/trop480/` holds an exact tropidist certification of d_X = d_Z = 26 for [[480,16,26]] (board commit b03ec2d4077f of tropidist, 16 threads per side, 5185 s and 5474 s), on a qubit relabelling of the board entry; `logs/trop480/README.md` explains the relabelling and `logs/trop480/check_perm.py` checks it.
