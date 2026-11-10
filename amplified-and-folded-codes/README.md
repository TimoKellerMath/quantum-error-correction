# Amplified and folded codes for the qLDPC challenge board

Ten quantum codes submitted to the [qLDPC challenge board](https://github.com/unitaryfoundation/qldpc-challenge)
(one pull request per code), with their check matrices, the proofs of their distances, programs and logs.

- Eight CSS codes, each one distance-amplification step (Liang, Gu, Chen, Eisert, Wang,
  [arXiv:2609.37231](https://arxiv.org/abs/2609.37231), cited as [LGCE26]) of a board entry with a certified exact
  distance. Their distance is exactly twice that of the base (Proposition 4 below, which is Lemma 4 of [LGCE26]).
- Two non-CSS stabilizer codes, folds (symplectic halvings) of two-block CSS codes along a ZX-duality
  (Burton, Durso-Sabina, Brown, [arXiv:2406.09951](https://arxiv.org/abs/2406.09951), cited as [BDB24]; Lee,
  Okada, Maskara, Kasai, Zhou, [arXiv:2609.30069](https://arxiv.org/abs/2609.30069), cited as [LOMK26]).

| board slug | code | type | construction | $w$ | qubit degree X / Z / total | $kd^2/n$ | distance |
|---|---|---|---|---|---|---|---|
| `880-40-22` | $[[880,40,22]]$ | CSS | $[[4,2,2]]$-amplified `180-20-11` | 8 | 7 / 7 / 11 | 22.00 | exact (A) |
| `866-128-12` | $[[866,128,12]]$ | CSS | $[[4,2,2]]$-amplified `186-64-6` | 7 | 6 / 6 / 10 | 21.28 | exact (A), and tropidist directly |
| `610-40-22` | $[[610,40,22]]$ | CSS | $[[4,2,2]]$-amplified `126-20-11` | 9 | 8 / 8 / 12 | 31.74 | exact (A) |
| `696-48-24` | $[[696,48,24]]$ | CSS | $[[4,2,2]]$-amplified `144-24-12` | 10 | 9 / 9 / 13 | 39.72 | exact (A) |
| `606-48-22` | $[[606,48,22]]$ | CSS | $[[4,2,2]]$-amplified `126-24-11` | 10 | 9 / 9 / 13 | 38.34 | exact (A) |
| `656-64-24` | $[[656,64,24]]$ | CSS | $[[6,4,2]]$-amplified `96-16-12` | 11 | 10 / 10 / 16 | 56.20 | exact (A) |
| `618-48-24` | $[[618,48,24]]$ | CSS | $[[6,4,2]]$-amplified `90-12-12` | 11 | 10 / 10 / 16 | 44.74 | exact (A) |
| `616-56-22` | $[[616,56,22]]$ | CSS | $[[6,4,2]]$-amplified `90-14-11` | 12 | 10 / 10 / 16 | 44.00 | exact (A) |
| `72-6-9-s` | $[[72,6,9]]$ | stabilizer | fold of the gross code `144-12-12` | 6 | total 6 | 6.75 | exact (B) |
| `255-8-22` | $[[255,8,d]]$, $16\le d\le22$ | stabilizer | fold of `510-16-26` | 6 | total 6 | 15.18 at $d=22$ | bounds (B) |

$w$ is the largest check weight. (A): Proposition 4 together with the board's exact certificates of the base
distances (tropidist, verification level `solver`); no further input. (B): see Part B; for `72-6-9-s` four
independent exhaustive computations, for `255-8-22` a proved lower bound 16 and a witness of weight 22.
The slug `72-6-9` belongs to a CSS entry with the same parameters; `-s` is the suffix the board's stabilizer
entries carry in that situation (e.g. `42-6-7-s` beside `42-6-7`). Stabilizer and CSS codes are ranked on separate
boards.

## Reproduce

Requirements: Python 3.8 or later (standard library only), a C compiler with `__builtin_popcountll` (gcc, clang),
`make`. Optional: Magma for `src/gross_k0_mindist.m`.

```
make check     # rebuild all ten codes from their definitions and check them (about 1 s); log: logs/rebuild.log
make quick     # compile src/*.c and rerun the short searches of Part B (about 3 CPU-minutes), compared with logs/
make long      # the long searches of Part B (about 1 CPU-hour)
make tropidist-inputs   # the inputs of our tropidist runs, compared with logs/tropidist/inputs.sha256
```

`python3 src/rebuild.py --board DIR`, with `DIR` a checkout of the board repository, also compares the base
entries in `data/bases/` and the fold parents with the board's files, and the ten codes with their board
entries once these exist. The searches of Part B were run in a working copy in which `instances/` was called
`codes/`, hence the paths `codes/...` in the logs; the instance files are the same. Our runs used gcc 16.2.1 with
`-O3 -march=native` on x86-64, one thread per job.

## Part A: amplified codes

### Construction

Notation: a CSS code is a pair $(H_X,H_Z)$ of $\mathbb F_2$-matrices with $n$ columns and $H_XH_Z^{\mathsf T}=0$,
whose rows need not be independent; $k=n-\operatorname{rk}H_X-\operatorname{rk}H_Z$,
$d_Z=\min\{|u|:u\in\ker H_X\setminus\operatorname{rowsp}H_Z\}$,
$d_X=\min\{|v|:v\in\ker H_Z\setminus\operatorname{rowsp}H_X\}$. For a qubit $q$, $\deg_X(q)$ and $\deg_Z(q)$ are
the numbers of rows of $H_X$, $H_Z$ containing $q$. $\mathbf 1\in\mathbb F_2^N$ is the all-one vector, $e_a$ the
unit vectors, $\langle\cdot,\cdot\rangle$ the standard bilinear form.

**Definition 1.** Let $(H_X,H_Z)$ be a CSS code on $n$ qubits with rows $H_X[0],\dots,H_X[m_X-1]$ and
$H_Z[0],\dots,H_Z[m_Z-1]$, and let $N\ge2$ be even. The code $\mathrm{Amp}_N(H_X,H_Z)$ has the qubits
- central $(q,a)$, $0\le q<n$, $0\le a<N$, with index $Nq+a$;
- $t_j$ for each Z-check $j$, index $Nn+j$;
- $s_i$ for each X-check $i$, index $Nn+m_Z+i$;

and the checks (in this order)
- X-checks $(c_i,a)=\{(q,a):q\in H_X[i]\}\cup\{s_i\}$ and $(q)^X=\{(q,a):0\le a<N\}\cup\{t_j:q\in H_Z[j]\}$;
- Z-checks $(z_j,a)=\{(q,a):q\in H_Z[j]\}\cup\{t_j\}$ and $(q)^Z=\{(q,a):0\le a<N\}\cup\{s_i:q\in H_X[i]\}$.

This is the code of [LGCE26, Eqs. (14), (15)] for the amplifier $[[N,N-2,2]]$ with $G_X=G_Z=\mathbf 1^{\mathsf T}$:
the truncation $T_3\to T_2\to T_1$ of the tensor product of the base complex
$Q_2\xrightarrow{H_Z^{\mathsf T}}Q_1\xrightarrow{H_X}Q_0$ with
$\mathcal A_N\colon\mathbb F_2\xrightarrow{\mathbf 1}\mathbb F_2^N\xrightarrow{\mathbf 1^{\mathsf T}}\mathbb F_2$,
with qubits $T_2=(Q_1\otimes A_1)\oplus(Q_2\otimes A_0)\oplus(Q_0\otimes A_2)$ in the order of [LGCE26, Eq. (13)].
$N=4$ is the $[[4,2,2]]$ amplifier, $N=6$ the $[[6,4,2]]$ amplifier. Our data files use exactly these indices and
this order of checks; `src/rebuild.py` rebuilds every amplified code from the board's base entry and the list of
base rows left out, and compares row by row.

**Lemma 2.** In the situation of Definition 1 with $N\ge4$:
1. $\mathrm{Amp}_N(H_X,H_Z)$ is a CSS code with $n'=Nn+m_X+m_Z$ and $k'=(N-2)k$. If the rows of $H_X$ and of
   $H_Z$ are linearly independent, $n'=(N+1)n-k$.
2. $(c_i,a)$ has weight $|H_X[i]|+1$, $(q)^X$ weight $N+\deg_Z(q)$, $(z_j,a)$ weight $|H_Z[j]|+1$, $(q)^Z$ weight
   $N+\deg_X(q)$. Hence the largest check weight is
   $w'=\max(w_X+1,\,w_Z+1,\,N+q_Z,\,N+q_X)$, with $w_X,w_Z$ the largest row weights and $q_X,q_Z$ the largest
   column weights of $H_X,H_Z$.
3. The qubit degrees (X, Z) are $(\deg_X(q)+1,\deg_Z(q)+1)$ for $(q,a)$, $(N,|H_X[i]|)$ for $s_i$ and
   $(|H_Z[j]|,N)$ for $t_j$.

*Proof.* Orthogonality, by cases: $(c_i,a)$ and $(z_j,a')$ meet in $|H_X[i]\cap H_Z[j]|$ central qubits if
$a=a'$ and in none otherwise, an even number as $H_XH_Z^{\mathsf T}=0$; $(c_i,a)$ and $(q)^Z$ share $(q,a)$ and
$s_i$ if $q\in H_X[i]$ and nothing otherwise; $(q)^X$ and $(z_j,a)$ likewise with $t_j$; $(q)^X$ and $(q')^Z$ meet
in $N$ central qubits if $q=q'$ and in none otherwise. $\ker\tilde H_X/\operatorname{rowsp}\tilde H_Z$ is the
homology of the truncation at $T_2$, which equals that of the full tensor product, because it only depends on the
two differentials at $T_2$. By the Künneth formula over a field it is
$\bigoplus_{i+j=2}H_i(\mathcal Q)\otimes H_j(\mathcal A_N)$, the isomorphism being induced by
$[u]\otimes[\ell]\mapsto[u\otimes\ell]$ ($u\otimes\ell$ placed on the central qubits, $(u\otimes\ell)(q,a)=u_q\ell_a$).
Here $H_0(\mathcal A_N)=\mathbb F_2/\langle\mathbf 1^{\mathsf T}\mathbb F_2^N\rangle=0$ and
$H_2(\mathcal A_N)=\ker(\mathbf 1)=0$, while $H_1(\mathcal A_N)=\mathbf 1^\perp/\langle\mathbf 1\rangle$ has
dimension $N-2$ because $\mathbf 1\in\mathbf 1^\perp$ for even $N$. So $k'=(N-2)k$, and the same holds for the
X-logicals. The rest is counting; with independent rows $m_X+m_Z=n-k$. $\square$

### Distance

**Lemma 3 (contraction).** In the situation of Definition 1, let $\omega\in\mathbb F_2^N$ have even weight. For
$z\in\mathbb F_2^{n'}$ let $z_a\in\mathbb F_2^n$ be its restriction to the qubits $(\cdot,a)$ and
$c_\omega(z)=\sum_{a:\omega_a=1}z_a$. Then $c_\omega(\ker\tilde H_X)\subseteq\ker H_X$,
$c_\omega(\operatorname{rowsp}\tilde H_Z)\subseteq\operatorname{rowsp}H_Z$ and
$c_\omega(u\otimes\ell)=\langle\omega,\ell\rangle u$.

*Proof.* Let $\tilde H_Xz=0$. The check $(c_i,a)$ gives $(H_Xz_a)_i=z(s_i)$; summing over the $|\omega|$ values
$a$ with $\omega_a=1$ gives $(H_Xc_\omega(z))_i=|\omega|\,z(s_i)=0$. On the generators of
$\operatorname{rowsp}\tilde H_Z$: $c_\omega((z_j,a))$ is $H_Z[j]$ if $\omega_a=1$ and $0$ otherwise, and
$c_\omega((q)^Z)=|\omega|e_q=0$. The last claim is the definition. $\square$

**Proposition 4 (exact doubling).** Let $(H_X,H_Z)$ be a CSS code with $k\ge1$ and arbitrary rows, and let
$N\ge4$ be even. Then $d_Z(\mathrm{Amp}_N(H_X,H_Z))=2d_Z$ and $d_X(\mathrm{Amp}_N(H_X,H_Z))=2d_X$.

*Proof.* Lower bound. Let $z\in\ker\tilde H_X\setminus\operatorname{rowsp}\tilde H_Z$ and let
$[\ell_1],\dots,[\ell_{N-2}]$ be a basis of $\mathbf 1^\perp/\langle\mathbf 1\rangle$. By the isomorphism in the
proof of Lemma 2, $z-\sum_ju_j\otimes\ell_j\in\operatorname{rowsp}\tilde H_Z$ for some $u_j\in\ker H_X$ with
$u_i\notin\operatorname{rowsp}H_Z$ for some $i$. The form $\langle\cdot,\cdot\rangle$ induces a perfect pairing on
$\mathbf 1^\perp/\langle\mathbf 1\rangle$ (its radical in $\mathbf 1^\perp$ is
$\mathbf 1^\perp\cap(\mathbf 1^\perp)^\perp=\langle\mathbf 1\rangle$), so there is $\omega\in\mathbf 1^\perp$
with $\langle\omega,\ell_j\rangle=\delta_{ij}$. Then $\omega'=\omega+\mathbf 1$ has even weight and
$\langle\omega',\ell_j\rangle=\delta_{ij}$, as every $\ell_j$ has even weight. By Lemma 3, $c_\omega(z)$ and
$c_{\omega'}(z)$ are congruent to $u_i$ modulo $\operatorname{rowsp}H_Z$, so both have weight at least $d_Z$. As the
supports of $\omega$ and $\omega'$ partition $\{0,\dots,N-1\}$,
$$2d_Z\le|c_\omega(z)|+|c_{\omega'}(z)|\le\sum_a|z_a|\le|z|.$$
Upper bound. For $u\in\ker H_X\setminus\operatorname{rowsp}H_Z$ with $|u|=d_Z$ and $\ell=e_0+e_1$, which lies in
$\mathbf 1^\perp\setminus\langle\mathbf 1\rangle$ since $N\ge4$, $u\otimes\ell$ is a nontrivial logical of weight
$2d_Z$ by Lemma 2. The X-statement follows by exchanging the roles of $H_X$ and $H_Z$ (and of $s$ and $t$), which
maps Definition 1 to itself. $\square$

Proposition 4 is [LGCE26, Lemma 4] for the amplifiers $[[N,N-2,2]]$: there, if every logical class of the
amplifier has a family of $\Sigma_P$ representatives with overlap at most $o_P$, then
$\lceil\alpha_Pd_P\rceil\le\tilde d_P\le d_Pd_P^A$ with $\alpha_X=\Sigma_Z/o_Z$, $\alpha_Z=\Sigma_X/o_X$, separately
for $P=X$ and $P=Z$, for every CSS base with $k>0$ whose rows need not be independent; the families
$\{\omega,\omega+\mathbf 1\}$ give $\Sigma_P=2$, $o_P=1$, and $d^A=2$, so $\tilde d_P=2d_P$. The proof there is the
same contraction argument; [LGCE26] attribute the overlap bound to Audoux and Couvreur,
[arXiv:1512.07081](https://arxiv.org/abs/1512.07081) (we did not check the statement numbers given there). What
the version above adds is the proof written out in the coordinates of Definition 1, which are those of our data
files, so that the hypotheses reduce to what `src/rebuild.py` checks: the code is $\mathrm{Amp}_N$ of a
presentation of the base entry, $N\in\{4,6\}$, and $k\ge1$. Lemma 5 below, which follows from the same proof,
describes all logical operators of weight $2d$; we use it to test our searches.

**Lemma 5 (minimum-weight logicals).** In the situation of Proposition 4, the Z-logicals of
$\mathrm{Amp}_N(H_X,H_Z)$ of weight $2d_Z$ are exactly the vectors $\ell\otimes(e_a+e_b)$, $a\ne b$, with $\ell$ a
Z-logical of the base of weight $d_Z$; their number is $\binom N2M_Z$, where $M_Z$ is the number of Z-logicals of
weight $d_Z$ of the base. Dually for X.

*Proof.* If $|z|=2d_Z$, all inequalities in the proof of Proposition 4 are equalities: $z$ vanishes on the qubits
$s_i$, $t_j$; $|c_\omega(z)|=|c_{\omega'}(z)|=d_Z$; and the $z_a$ with $a$ in the support of $\omega$ have pairwise
disjoint supports, as have those with $a$ in the support of $\omega'$. As $z$ vanishes on the $s_i$, the checks
$(c_i,a)$ give $z_a\in\ker H_X$ for all $a$. Some $z_{a_0}$ with $\omega_{a_0}=1$ is not in
$\operatorname{rowsp}H_Z$, since $c_\omega(z)$ is not; so
$d_Z\le|z_{a_0}|\le|c_\omega(z)|=d_Z$ and $z_a=0$ for the other $a$ with $\omega_a=1$. Likewise
$c_{\omega'}(z)=z_{b_0}$ for a single $b_0$ with $\omega'_{b_0}=1$, and $z_a=0$ for all $a\ne a_0,b_0$. As $z$
vanishes on the $t_j$, the checks $(q)^X$ give $\sum_az_a=0$, so $z_{a_0}=z_{b_0}=:\ell$ and
$z=\ell\otimes(e_{a_0}+e_{b_0})$. Conversely these vectors are nontrivial logicals of weight $2d_Z$ (proof of
Proposition 4, with $e_a+e_b$ for $e_0+e_1$), and distinct pairs $(\ell,\{a,b\})$ give distinct vectors. $\square$

### The eight codes

In each case the retained base rows are linearly independent and span the row spaces of the base, so
$n'=(N+1)n-k$ (Lemma 2), and they were chosen to minimize the check weight; the choice is recorded in
`data/codes/<slug>.json` (`drop_X`, `drop_Z`: indices of the board's rows left out). The witnesses in the board
entries are $u\otimes(e_0+e_1)$ for the witnesses $u$ of the base entries. The base distances are certified exact on
both sides by the board's certificates `certs/<base>.json` (board commit `8042f9da`), all made with tropidist
(Cheng, Guo, Pan, [arXiv:2611.08247](https://arxiv.org/abs/2611.08247), github.com/Fanerst/tropidist at commit
`b03ec2d4077f`), an exhaustive search whose completeness is proved in that paper; their verification level is
`solver` (no replayable proof object). With Proposition 4 this proves $d_X=d_Z=2d$.

| slug | base entry (authors on the board; source) | $N$ | retained X / Z rows | base $[[n,k,d]]$, certificate |
|---|---|---|---|---|
| `880-40-22` | `180-20-11`, two-block code over $C_{15}\times S_3$ (G. Lee; arXiv:2610.13880) | 4 | 80 / 80 of 90 | $[[180,20,11]]$, tropidist |
| `866-128-12` | `186-64-6` (@MathysRennela; Construction I of arXiv:2610.19451) | 4 | 61 / 61 of 62 | $[[186,64,6]]$, tropidist |
| `610-40-22` | `126-20-11`, generalized bicycle code on $\mathbb Z_{63}$ (@FarLab) | 4 | 53 / 53 of 63 | $[[126,20,11]]$, tropidist |
| `696-48-24` | `144-24-12`, two-block code over $C_6\times C_6\times C_2$ (G. Lee; arXiv:2610.13880) | 4 | 60 / 60 of 72 | $[[144,24,12]]$, tropidist |
| `606-48-22` | `126-24-11`, two-block code over $C_{63}$ (G. Lee; arXiv:2610.13880) | 4 | 51 / 51 of 63 | $[[126,24,11]]$, tropidist |
| `656-64-24` | `96-16-12`, two-block code (@FarLab) | 6 | 40 / 40 of 48 | $[[96,16,12]]$, tropidist |
| `618-48-24` | `90-12-12`, two-block code (@FarLab) | 6 | 39 / 39 of 45 | $[[90,12,12]]$, tropidist |
| `616-56-22` | `90-14-11`, two-block code (@FarLab) | 6 | 38 / 38 of 45 | $[[90,14,11]]$, tropidist |

The base check matrices and witnesses in `data/bases/` are copied from the board (Apache License 2.0, see below);
`src/rebuild.py --board DIR` checks that they coincide with the board's files.

For `616-56-22` the check weight 12 cannot be lowered by another choice of retained rows: the 7 rows of each base
matrix that are left out meet at most $7\cdot6=42$ of the 45 columns of degree 6, so a column of degree 6 remains
and $w'\ge N+6=12$ by Lemma 2.

Further evidence (computed; not used by the proof):
- tropidist, run by us at the board's commit `b03ec2d4077f` in the board's encoding (side X: $H_Z$ against a basis
  of the Z-logicals; side Z mirrored), one thread, 600 s per run (`logs/tropidist/`). `866-128-12`: $d_X=d_Z=12$
  exactly, 2.3 s per side; this proves $d=12$ without Proposition 4. `610-40-22`, side X: no result within 600 s,
  but the proved lower bound $d_X\ge20$ at the time limit. We did not try the other sides and codes.
- Information-set search (Prange with sums of two rows, our program, not included), $10^4$ trials per side with
  seeds 1000 to 1015 and $10^5$ with seeds 5000 to 5015 (`logs/amplified/isd_main.log`, `isd_deep.log`, one line per
  run: seed, lightest weight found, histogram): the lightest logical found has weight exactly $2d$ on all sixteen
  sides, and every distinct lightest word found (rechecked independently of the search) has the form of Lemma 5
  (`logs/amplified/isd_summary.json`: per side the trials, the lightest weight, the number of distinct lightest
  words found and the number $\binom N2M_{\rm base}$ predicted by Lemma 5, with $M_{\rm base}$ found by a saturated
  search on the base).
  For `866-128-12` the search found all $\binom42\cdot3813=22878$ words of weight 12 on each side. At $2d=24$ such a
  search is weak evidence by itself: a single fixed word of weight 23 would be hit less than once on average.

## Part B: folds

### Construction

**Two-block parents.** For a finite abelian group $G$ and subsets $F_A,F_B\subset G$ the two-block code has the
qubits $(1,h)$ and $(2,h)$, $h\in G$, with indices $\iota(h)$ and $|G|+\iota(h)$, where $\iota$ numbers $G$ in
lexicographic order of the coordinates (for $G=\mathbb Z/12\times\mathbb Z/6$: $\iota(i,j)=6i+j$). The X-check $g$
acts on $(1,g+F_A)\cup(2,g+F_B)$, the Z-check $g$ on $(1,g-F_B)\cup(2,g-F_A)$. This is the board's convention
$H_X=(A\mid B)$, $H_Z=(B^{\mathsf T}\mid A^{\mathsf T})$. The two parents:
- the gross code, $G=\mathbb Z/12\times\mathbb Z/6$, $F_A=\{(3,0),(0,1),(0,2)\}$, $F_B=\{(0,3),(1,0),(2,0)\}$
  ($A=x^3+y+y^2$, $B=y^3+x+x^2$; Bravyi, Cross, Gambetta, Maslov, Rall, Yoder,
  [arXiv:2308.07915](https://arxiv.org/abs/2308.07915)), $[[144,12,12]]$, the board entry `144-12-12`;
- $G=\mathbb Z/255$, $F_A=\{0,26,76\}$, $F_B=\{0,30,66\}$, $[[510,16,26]]$, the board entry `510-16-26`
  (V. Russo); its exact distance 26 is Theorem 4.1 of `weight6-exact-distances/exact_distances.pdf` in this
  repository (commit `2be467a7bd36`), a computer-assisted proof.

With these indices both matrices coincide with the board's files (`src/rebuild.py --board`).

**Definition 6 (fold).** Let $C=(H_X,H_Z)$ be a CSS code on $2m$ qubits and $\Psi$ a fixed-point-free involution of
the qubits that maps the set of supports of the X-checks onto the set of supports of the Z-checks. Let
$\sigma_0<\dots<\sigma_{m-1}$ be the qubits $q$ with $q<\Psi(q)$. For $a\in\mathbb F_2^{2m}$ let $\varphi(a)$ be
the Pauli operator on $m$ qubits with $X$ on $\{i:a_{\sigma_i}=1\}$ and $Z$ on $\{i:a_{\Psi(\sigma_i)}=1\}$ (so $Y$
on their intersection). The fold $S_\Psi$ is the stabilizer group generated by the $\varphi(a)$, $a$ a row of $H_X$;
its generators in our data are these, in the order of the X-checks.

This is the symplectic halving of [LOMK26, Def. 3, Eqs. (13) to (15)] with the pairing $\pi=\Psi$; [LOMK26] credit
the halving to [BDB24], where Thm. 3.5 shows that a CSS code on $2m$ qubits has a fixed-point-free involutory
ZX-duality exactly when it is the double of a code on $m$ qubits. The board entry `72-6-8-b` (a genon code,
arXiv:2611.09538) is the fold of the gross code along another such involution (by our computation, equal up to a
relabelling of the qubits; see the dead ends of `72-6-9-s` below).

**Lemma 7.** In the situation of Definition 6, for a Pauli operator $p=(x\mid z)$ on $m$ qubits let
$\lambda(p)\in\mathbb F_2^{2m}$ be given by $\lambda(p)_{\sigma_i}=z_i$ and $\lambda(p)_{\Psi(\sigma_i)}=x_i$. Then
$\lambda$ is a linear bijection, the symplectic product satisfies $\langle\varphi(a),p\rangle=a\cdot\lambda(p)$,
and $\lambda(\varphi(a))=a\circ\Psi$. Consequently $p$ commutes with $S_\Psi$ iff $\lambda(p)\in\ker H_X$, and
$p\in S_\Psi$ iff $\lambda(p)\in\operatorname{rowsp}H_Z$; $\lambda$ maps the nontrivial logical operators of
$S_\Psi$ bijectively onto the nontrivial Z-logicals of $C$, $k(S_\Psi)=k(C)/2$, and
$|\lambda(p)|=2y(p)+t(p)$ where $y(p)$ is the number of $Y$ factors of $p$ and $t(p)$ the number of $X$ or $Z$
factors, while the Pauli weight is $y(p)+t(p)$. In particular $d(S_\Psi)\ge\lceil d_Z(C)/2\rceil$.

*Proof.* $\langle\varphi(a),(x\mid z)\rangle=\sum_i(a_{\sigma_i}z_i+a_{\Psi(\sigma_i)}x_i)=a\cdot\lambda(p)$, as
every qubit is exactly one of $\sigma_i$, $\Psi(\sigma_i)$. By definition
$\lambda(\varphi(a))_{\sigma_i}=a_{\Psi(\sigma_i)}$ and $\lambda(\varphi(a))_{\Psi(\sigma_i)}=a_{\sigma_i}$, i.e.
$\lambda(\varphi(a))=a\circ\Psi$. The first equivalence follows as $S_\Psi$ is generated by the $\varphi(a)$ for the
rows $a$ of $H_X$. The permutation $\Psi$ maps $\operatorname{rowsp}H_X$ onto $\operatorname{rowsp}H_Z$, so
$\lambda(S_\Psi)=\{a\circ\Psi:a\in\operatorname{rowsp}H_X\}=\operatorname{rowsp}H_Z$, and the second equivalence
follows from the bijectivity of $\lambda$. Then $2k(S_\Psi)=\dim\ker H_X-\dim\operatorname{rowsp}H_Z=k(C)$; the
weights are read off the definition. $\square$

Lemma 7 is [LOMK26, Eqs. (17) to (20)] in our notation. The next lemma sharpens the last bound with the parity of
the parent and an exhaustive search over the logicals with few non-$Y$ factors.

**Lemma 8.** In the situation of Lemma 7 assume moreover that every column of $H_X$ has odd weight and
$d_Z(C)\ge D$. Then every nontrivial logical $p$ of $S_\Psi$ has an even number $t(p)$ of $X$ or $Z$ factors and
Pauli weight at least $(D+t(p))/2$. Hence, if $T$ is even and every nontrivial logical with $t(p)\le T$ has Pauli
weight at least $B$, then $d(S_\Psi)\ge\min(B,(D+T+2)/2)$.

*Proof.* $v=\lambda(p)$ is a nontrivial Z-logical of $C$ (Lemma 7), so $|v|\ge D$; summing the rows of
$H_Xv=0$ gives $\sum_q\mathrm{colwt}(q)v_q=0$, i.e. $|v|$ is even. As $|v|=2y(p)+t(p)$, $t(p)$ is even and the
Pauli weight is $y(p)+t(p)=(|v|+t(p))/2\ge(D+t(p))/2$. A nontrivial logical with $t(p)>T$ has $t(p)\ge T+2$. $\square$

**Lemma 9 (no single-qubit Clifford makes the code CSS).** Let $s_1,\dots,s_r$ generate a stabilizer group $S$
such that (i) for every $i$, the only elements of $S$ supported in $\operatorname{supp}s_i$ are $1$ and $s_i$ (up
to phase). If a product $U$ of single-qubit Clifford operators maps $S$ onto a CSS group, then $U$ maps every $s_i$
to a pure $X$- or pure $Z$-type operator.

*Proof.* Write $Us_iU^\dagger=ab$ with $a$ of $X$-type and $b$ of $Z$-type, both in $USU^\dagger$ (a CSS group is
the product of its $X$-type and its $Z$-type subgroup). Qubitwise $a$ is $I$ or $X$ and $b$ is $I$ or $Z$, so
$\operatorname{supp}(ab)=\operatorname{supp}a\cup\operatorname{supp}b$. Single-qubit Cliffords preserve supports,
so $U^\dagger aU$ and $U^\dagger bU$ are elements of $S$ supported in $\operatorname{supp}s_i$, hence in
$\{1,s_i\}$ by (i); as $ab=Us_iU^\dagger\ne1$, exactly one of $a$, $b$ is $Us_iU^\dagger$ and the other is $1$.
$\square$

The board's verifier tests the generator-level condition: whether single-qubit Cliffords exist that make every
generator pure, a parity system with one unknown per generator. For both folds `src/rebuild.py` shows that this
system has no solution and that (i) holds for every generator; by Lemma 9 no single-qubit Clifford makes the group
CSS.

### `72-6-9-s`: the $[[72,6,9]]$ fold of the gross code

The involution: $\Psi(1,h)=(1,c_1+\psi(h))$, $\Psi(2,h)=(2,c_2+\psi(h))$ with the automorphism $\psi$ of
$\mathbb Z/12\times\mathbb Z/6$ given by $\psi(1,0)=(3,4)$, $\psi(0,1)=(10,3)$ ($\psi^2=1$), $c_1=(1,2)$,
$c_2=(5,4)$; it preserves the two blocks. The fold has $n=72$, $k=6$ (rank 66) and 72 generators, all of weight 6,
none with a $Y$ factor.

$d=9$ (proved, computer-assisted; four independent computations):
1. Meet in the middle (`src/mitm.c`): every nontrivial logical of weight at most $2h$ is $ab$ with $a,b$ of weight
   at most $h$, equal syndromes and different logical labels; for $h=4$ none of the 84,965,491 Pauli operators of
   weight at most 4 pairs up, so there is no nontrivial logical of weight at most 8 (24.8 s,
   `logs/folds/mitm_gross_c4.log`). Validation (`logs/folds/mitm_validation.log`, `mitm_board_72-6-8-b.log`):
   $[[5,1,3]]$: 3; the board entries `13-1-5`, `17-1-7`, `21-1-7`: 5, 7, 7; the board entry `72-6-8-b`: 8.
2. Magma V2.30-1, `MinimumDistance` of the additive code over $\mathrm{GF}(4)$ (Brouwer-Zimmermann type): 9, in
   411 s (`src/gross_k0_mindist.m`, `logs/folds/magma_gross_k0.log`).
3. Lemma 8 with $D=12$ (the gross code has $d=12$, certified on the board, `certs/144-12-12.json`; all column
   weights of $H_X$ are 3) and $T=4$: `src/fewnony2.c` gives the least weight 12 for $t=0$ and $t=4$ and no
   nontrivial logical with $t\in\{1,2,3\}$ (0.5 s, `logs/folds/fewnony2_validation.log`), so
   $d\ge\min(12,(12+6)/2)=9$.
4. tropidist (github.com/Fanerst/tropidist at `b03ec2d4077f`), run by us in the board's encoding for stabilizer
   codes (three columns $X_i,Z_i,Y_i$ per qubit, the generators and a basis of the logicals as symplectic
   products): exact 9 in 1.8 s, one thread (`logs/tropidist/`).

Upper bound: the witness in the board entry. There are exactly 24 nontrivial logicals of weight 9, each with three
$Y$ factors (`src/count_odd.c`, exhaustive, 878 s, `logs/folds/count9_gross_k0.log`; the information-set search
`src/isd.c` with $10^5$ trials, seed 20261110, found the same 24, `logs/folds/isd_gross_k0.log`); they are in
`data/witnesses/72-6-9-s.json`.

Not CSS under single-qubit Cliffords (Lemma 9 and `src/rebuild.py`). Up to translations the gross code has five
classes of affine ZX-involutions: the block-preserving class of this fold and four block-swapping classes with
$\psi=-1$, whose folds have $d=8$ (three classes, `logs/folds/mitm_gross_c1.log` to `mitm_gross_c3.log`, one of them
the board's genon `72-6-8-b`) and $d\le6$ (one class, a witness of weight 6 in `data/witnesses/gross_c0_d6.json`).

### `255-8-22`: the $[[255,8]]$ fold of `510-16-26`

The involution swaps the blocks: $\Psi(1,h)=(2,15+154h)$, $\Psi(2,h)=(1,154(h-15))$ in $\mathbb Z/255$
($154^2\equiv1$). Fold qubit $h$ is the pair $\{(1,h),\Psi(1,h)\}$. The fold has $n=255$, $k=8$ (rank 247) and
255 generators, 204 of weight 6 and 51 of weight 5 with one $Y$ factor each. The translations $h\mapsto h+5$ and the
multiplication $h\mapsto101h$ map the generator set onto itself (`src/rebuild.py`); they generate a group of order
102.

$16\le d\le22$:
- Lower bound 16 (proved, computer-assisted): Lemma 8 with $D=26$ (the exact distance of `510-16-26`, see above;
  all column weights 3) and $T=4$. `src/fewnony2.c` enumerates every non-$Y$ support of size $t\le4$ whose smallest
  element is less than 5, which suffices because the translations by multiples of 5 are automorphisms: least
  weights 81 ($t=0$), 59 ($t=2$), 36 ($t=4$), and no nontrivial logical with $t\in\{1,3\}$ (47.4 s,
  `logs/folds/fewnony2_B510_t0_T4_A5.log`). So $d\ge\min(36,(26+6)/2)=16$.
- Upper bound 22: the witness in the board entry, one of 102 nontrivial logicals of weight 22 in
  `data/witnesses/255-8-22.json`.
- Evidence that $d=22$ (computed, heuristic as a search): three runs of `src/isd.c` with $10^6$ trials each
  (seeds 101 and 303: random qubit order in a random single-qubit Clifford frame per trial; seed 202: random
  order of the $2n$ symplectic coordinates), 31 CPU-minutes in total, found nothing lighter than 22
  (`logs/folds/isd_B510_t0.log`, `isd_B510_t0_Q2.log`). The 102 distinct words of weight 22 form two complete
  orbits of size 51 under the group of order 102; the two runs of the first kind hit them 625 times, 6.1 times per
  word (`logs/folds/c2_B510_t0.log`). Exact orbit folds (`src/orbit_fold.py` and `src/mitm.c`,
  `logs/folds/orbit_fold_B510.log`, `mitm_orbitfold_B510.log`): no nontrivial logical invariant under the rotation
  subgroup of order 3 has Pauli weight at most 24, none invariant under the subgroup of order 17 has weight at most
  238, and none is invariant under the subgroup of order 51. The reflections were not searched exactly. This does
  not prove $d=22$.

Not CSS under single-qubit Cliffords: 51 qubits carry $X$, $Y$ and $Z$ among the generators (Lemma 9 and
`src/rebuild.py`).

## Files

- `src/rebuild.py`: rebuilds and checks all ten codes (standard library only); `src/same_results.py` compares
  reruns with our logs; `src/tropidist_inputs.py` writes the tropidist inputs.
- `src/mitm.c`, `src/fewnony2.c`, `src/count_odd.c`, `src/isd.c`, `src/orbit_fold.py`: the searches of Part B
  (usage in the head of each file); `src/gross_k0_mindist.m`: the Magma computation.
- `data/codes/<slug>.json`: the submitted codes: `checks` exactly as in the board entry (`X`, `Z` lists of
  supports, or `S` a list of generators `{"X": [...], "Z": [...]}`, qubits from 0), `witness`, `construction`
  (amplified: base slug, $N$, base rows left out; folds: group, $F_A$, $F_B$, involution).
- `data/bases/<base>.json`: checks and witnesses of the eight base entries, copied from the board.
- `data/witnesses/`: all logicals of weight 9 of `72-6-9-s`, all 102 of weight 22 found for `255-8-22`, and a
  witness of weight 6 for the fold `gross_c0`.
- `instances/`: the inputs of the C programs: `n r L`, then $r$ generators and $L$ logical representatives as
  strings of $2n$ bits ($X$ part, then $Z$ part). `gross_c4` is `72-6-9-s`, `B510_c0` is `255-8-22` (its
  `.json` also holds the involution and $\Sigma$), `B510_c0_orb*` are its orbit folds, `gross_c0` to `gross_c3` the
  other folds of the gross code, `board_72-6-8-b` and `test_*` are validation instances made from board entries,
  `test513` is $[[5,1,3]]$.
- `logs/rebuild.log`: output of `make check`, `logs/environment.txt` the system it ran on; `logs/folds/`,
  `logs/amplified/`: our runs; `logs/tropidist/`: the tropidist records (`<slug>.<side>.json`, written by the
  tool), wall times (`runs.tsv`) and the checksums of their inputs (`inputs.sha256`).

## Credits and licenses

The constructions of the ten codes, the programs, the computations and the texts of this directory were made with
Claude Opus 5.5 (Anthropic) in Claude Code. The constructions are those of [LGCE26] (amplification) and of [BDB24] and
[LOMK26] (halving); the base codes are board entries by G. Lee (arXiv:2610.13880), @MathysRennela
(arXiv:2610.19451), @FarLab, S. Bravyi et al. (arXiv:2308.07915) and V. Russo, as listed above. The files in
`data/bases/` and the instances `board_72-6-8-b.code`, `test_13-1-5.code`, `test_17-1-7.code`, `test_21-1-7.code`
are derived from files of github.com/unitaryfoundation/qldpc-challenge (commit `8042f9da`), which is licensed under
the Apache License 2.0; everything else is under the MIT license of this repository. tropidist is by Cheng, Guo and
Pan (arXiv:2611.08247).
