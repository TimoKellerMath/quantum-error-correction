# Surface codes from regular maps: nine codes with exact distances

Nine quantum CSS codes on closed orientable hyperbolic surfaces: eight from orientably-regular maps of type
$\{8,8\}$ and $\{6,6\}$ in M. Conder's lists, and a two-sheeted cyclic cover of a $\{4,5\}$ map. For all
nine the minimum distance is determined exactly, on both sides, by three independent programs of ours and by
tropidist. The codes are submitted to the qLDPC challenge board
(https://github.com/unitaryfoundation/qldpc-challenge), one pull request each, as the entries named by the slugs
below; the matrices in `data/codes/` are those of the entries. Programs, computations and first drafts of the texts were made with
Claude Opus 5.5 (Anthropic) in Claude Code.

| slug | $[[n,k,d]]$ | type | map (Conder) | kind | $(V,E,F)$ | genus | $w$ | $kd^2/n$ | min.-weight logicals X / Z |
|---|---|---|---|---|---|---|---|---|---|
| `864-434-6` | $[[864,434,6]]$ | $\{8,8\}$ | C217.22 (mirror image) | chiral | (216, 864, 216) | 217 | 8 | 18.08 | 21744 / 21744 |
| `720-362-6` | $[[720,362,6]]$ | $\{8,8\}$ | C181.29 | chiral | (180, 720, 180) | 181 | 8 | 18.10 | 24600 / 24600 |
| `936-470-6` | $[[936,470,6]]$ | $\{8,8\}$ | C235.23 (mirror image) | chiral | (234, 936, 234) | 235 | 8 | 18.08 | 26364 / 26364 |
| `648-326-6` | $[[648,326,6]]$ | $\{8,8\}$ | C163.13 | chiral | (162, 648, 162) | 163 | 8 | 18.11 | 24732 / 24732 |
| `640-322-6` | $[[640,322,6]]$ | $\{8,8\}$ | C161.44 (mirror image) | chiral | (160, 640, 160) | 161 | 8 | 18.11 | 24960 / 24960 |
| `1000-502-6` | $[[1000,502,6]]$ | $\{8,8\}$ | C251.45 (mirror image) | chiral | (250, 1000, 250) | 251 | 8 | 18.07 | 22500 / 22500 |
| `294-100-6` | $[[294,100,6]]$ | $\{6,6\}$ | C50.3 (mirror image) | chiral | (98, 294, 98) | 50 | 6 | 12.24 | 3381 / 3381 |
| `486-164-6` | $[[486,164,6]]$ | $\{6,6\}$ | R82.22 | reflexible | (162, 486, 162) | 82 | 6 | 12.15 | 2835 / 2835 |
| `720-74-8` | $[[720,74,8]]$ | $\{4,5\}$ | double cover of R19.3 | cover | (288, 720, 360) | 37 | 5 | 6.58 | 100 / 6108 |

$w$ is the maximal check weight; every qubit lies in two X-checks and two Z-checks. All codes have
$d_X=d_Z=d$. The last column is the number of logical operators of weight $d$ of each type (as vectors).

## Construction

For a finite group $G=\langle r,s\rangle$ with $(rs)^2=1$ the map $M(r,s)$ has darts $G$ (right
multiplication), faces $g\langle r\rangle$, vertices $g\langle s\rangle$ and edges
$g\langle rs\rangle=\{g,grs\}$. Qubits are edges, X-checks are vertices, Z-checks are faces, incidences
counted mod 2. Thus $H_X=\partial_1$ and $H_Z=\partial_2^{\mathsf T}$ of the cellular chain complex of the
surface over $\mathbb F_2$, $k=\dim H_1(\Sigma;\mathbb F_2)=2-\chi=2g$, $d_Z$ is the length of a shortest
homologically nontrivial cycle of the graph of the map and $d_X$ that of the dual graph (Breuckmann and Terhal,
arXiv:1506.04029). If $r,s$ have orders $p,q$, the map has type $\{p,q\}$: $p$-gonal faces, $q$-valent
vertices, so X-checks of weight $q$ and Z-checks of weight $p$.

For each regular map, $G=\langle r,s\mid r^p,s^q,(rs)^2,\dots\rangle$ with the relators below has order $2E$
(coset enumeration by Magma and by `src/rebuild.py`). The relators are given in the orientation of the dart
permutations in `data/maps.json`; for an entry marked "mirror image" Conder's relators hold for
$(r^{-1},s^{-1})$. This does not change the code: $M(r,s)$ and $M(r^{-1},s^{-1})$ have the same faces and
vertices (as sets of darts, since $g\langle r^{-1}\rangle=g\langle r\rangle$), and
$\{g,grs\}\mapsto\{gr,gs^{-1}\}$ is a bijection between their edges that preserves incidences
($grs\cdot r=gs^{-1}$ by $rsr=s^{-1}$). Relators were dropped where the group without them still has order
$2E$. Conder's lists: *Regular orientable maps of genus 2 to 301* and *Chiral rotary maps of genus 2 to 301*,
https://www.math.auckland.ac.nz/~conder/ (`RegularOrientableMaps301.txt`, `ChiralMaps301.txt`, downloaded
2026-11-10); for a reflexible entry the rotation group is $\langle R,S\mid w_i,\,w_i(R^{-1},S^{-1})\rangle$,
and the mirrored words turned out to be redundant.

| code | Conder's entry | data/maps.json vs. Conder | relators besides $r^p, s^q, (rs)^2$ | dropped as redundant |
|---|---|---|---|---|
| `864-434-6` | C217.22: Type {8,8}_12, Order 1728, Self-dual | mirror image | $s^{-1}r^{3}s^{-1}r^{-2}s^{-4}r$, $s^{-1}rs^{-2}r^{-3}s^{2}r^{-1}s^{-1}r$ | none |
| `720-362-6` | C181.29: Type {8,8}_60, Order 1440, Self-dual | same | $s^{-1}r^{3}s^{-1}r^{-2}s^{4}r$, $rs^{-1}r^{-2}sr^{3}s^{-1}rs^{-2}rs^{-1}$ | $s^{-1}r^{-1}s^{2}r^{3}s^{2}r^{-1}s^{-1}r$ |
| `936-470-6` | C235.23: Type {8,8}_78, Order 1872, Self-dual | mirror image | $srs^{-2}r^{-3}s^{-2}rsr^{-1}$, $s^{-1}r^{-1}sr^{-1}sr^{-2}s^{-1}rs^{-1}rs^{-1}$ | $s^{-1}r^{3}s^{-1}r^{-2}s^{-4}r$ |
| `648-326-6` | C163.13: Type {8,8}_6, Order 1296, Self-dual | same | $s^{-1}r^{-1}s^{2}r^{3}s^{2}r^{-1}s^{-1}r$, $s^{-1}r^{-1}srs^{-1}r^{2}s^{-1}rs^{-2}r$ | $sr^{-3}sr^{2}s^{4}r^{-1}$ |
| `640-322-6` | C161.44: Type {8,8}_10, Order 1280, Self-dual | mirror image | $s^{-1}r^{2}s^{-1}r^{-3}s^{2}r^{-1}s^{-2}$, $s^{2}r^{-2}sr^{-2}sr^{-1}s^{2}r^{-1}$ | $sr^{-3}sr^{2}s^{-4}r^{-1}$ |
| `1000-502-6` | C251.45: Type {8,8}_20, Order 2000, Self-dual | mirror image | $r^{-1}s^{3}r^{-3}s^{2}r^{-1}sr^{-1}$, $s^{2}r^{-2}sr^{-2}sr^{-1}s^{2}r^{-1}$ | none |
| `294-100-6` | C50.3: Type {6,6}_14, Order 588, Self-dual | mirror image | $r^{-1}sr^{-1}sr^{-3}sr^{-1}s^{-3}$ | none |
| `486-164-6` | R82.22: Type {6,6}_6, Order 1944, Self-dual | same | $srs^{-2}r^{3}s^{3}r^{-2}$, $rs^{-2}rs^{-1}r^{2}s^{-2}rs^{-1}r$, $srs^{-1}r^{-1}sr^{2}sr^{-1}s^{-1}rs$ | none |
| base of `720-74-8` | R19.3: Type {4,5}_8, Order 1440 | same | $sr^{-2}s^{2}rs^{-1}r^{2}s^{-1}rsr^{-1}s^{-1}rs$ | none |

All eight regular maps are self-dual (Conder's flag), so $d_X=d_Z$ for them.

**The cover `720-74-8`.** Base: R19.3 above, $|G|=720$, the $\{4,5\}$ code $[[360,38,8]]$ with
$(V,E,F)=(144,360,180)$. A voltage $\omega$ on the edges with values in $\mathbb Z/2$ is 1 exactly on the 41
edges $\{w,w\,rs\}$ for the words $w$ of `data/cover_voltage.txt` (each the shortlex-least word, in the order
$r<r^{-1}<s<s^{-1}$, of one of the two darts of the edge). The cover has darts $G\times\mathbb Z/2$,
$\tilde R(g,i)=(gr,i+\omega(g))$, $\tilde S(g,i)=(gs,i)$, where $\omega(g)$ is the voltage of the edge through
$g$; faces, vertices and edges are the cycles of $\tilde R$, $\tilde S$, $\tilde R\tilde S$, and
$(V,E,F)=(288,720,360)$, genus 37, $k=74$. `data/maps.json` stores the base darts and a cohomologous voltage with 87
edges (the one drawn in the census; $\omega+\delta f$ for a 0-cochain $f$ gives an isomorphic cover); the 41-edge
voltage was obtained from it by vertex flips. No orientably-regular $\{4,5\}$ map of genus 37 exists (Conder's
lists), so the cover is not a regular map.

## How the codes were found

A census of all orientably-regular maps of type $\{p,q\}$, $3\le p,q\le8$, $1/p+1/q<1/2$, with at most 2000
edges (for $\{8,8\}$: chiral maps up to 1000 edges), as the normal subgroups $N$ of
$\langle r,s\mid r^p,s^q,(rs)^2\rangle$ of index $2E$ in which $r,s,rs$ keep the orders $p,q,2$ (Magma
`LowIndexNormalSubgroups`), of the non-orientable regular maps (extended triangle group) and of 2000 cyclic covers
of the orientable maps with at most 1000 edges (random voltage classes with recorded seeds): 2174 orientable and
340 non-orientable maps (up to duality and mirror image) and 2000 covers, 7063 codes with their X/Z duals, all with
exact $d_X,d_Z$ by the two programs below. The census
programs (Magma and Python) are not part of this package. Against the board at commit `8042f9da` (2026-11-09; CSS
codes, dominance in $(n,k,d,w)$ within the cells of check weight and locality), 254 parameter points of the census
were undominated. Of the 1005 $\{8,8\}$ maps (up to mirror image), 19 have $d=6$ and none more, the smallest
with $n=640$; among the 920 $\{6,6\}$ maps, $d=6$ is the maximum, first reached at $n=294$; all ten sampled double covers of R19.3 are
$[[720,74,8]]$. Excluded as known from the literature: $[[896,386,6]]$ ($\{7,7\}$, R193.166) and
$[[960,194,8]]$ ($\{5,5\}$, R97.56), both plotted in Fig. 3.7 of N. P. Breuckmann's thesis (arXiv:1802.01520),
and $[[588,254,6]]$, $[[273,119,5]]$ (Table 3.1 there). The nine codes here were not found in Breuckmann and
Terhal (arXiv:1506.04029), Breuckmann et al. (arXiv:1703.00590), Breuckmann's thesis, Mahmoud and Rayan
(arXiv:2610.10948), Breuckmann and Burton (arXiv:2202.06647, App. B.1), Mahmoud, Ali and Rayan (arXiv:2504.07800),
the catalog of github.com/QEC-pages/Quantum_LDPC_Codes and qecdb.org (searched 2026-11-10); Conder's lists contain
the maps, not the codes.

## Distances: what is proved, and how

- **Upper bounds:** a logical operator of weight $d$ of each type in `data/parameters.json`.
- **Lower bounds.** Every column of $H_X$ and of $H_Z$ has weight 2, so a Z-type logical operator is a cycle of the
  graph of the map and an X-type one a cycle of the dual graph, and a minimum-weight one is a circuit. Three
  programs that share no code compute the least weight of a homologically nontrivial cycle exactly:
  - `src/graphlike_A/graphlike_A.c`: breadth-first trees from every vertex with logical annotations in
    $\mathbb F_2^k$ (a Horton-type candidate set); exactness is Proposition 5 and Lemma 7 of
    `src/graphlike_A/NOTES.md`, the count of minimum-weight logicals Proposition 8.
  - `src/graphlike_B/graphlike_B.c`: shortest paths in the $k$ double covers of the check graph defined by the
    logical operators of the other type, searched bit-sliced; exactness is Corollary 4 and Propositions 8, 9 of
    `src/graphlike_B/NOTES.md`; `--ref` adds an unpruned reference search.
  - `src/circuits_check.py`: enumeration of all circuits of length at most $d$ of the check graph and a
    row-space test for each (Lemma 12 of `src/graphlike_B/NOTES.md`); it also counts the minimum-weight logicals.

  Validation (computed, 2026-11-10): graphlike_A on 33 toric, rotated and planar surface codes of known
  distance, in 14223 comparisons with exhaustive enumeration on random graph-like codes, and on the 428 board
  entries with two graph-like sides (it agrees with all 302 exact board certificates among them); graphlike_B on
  the toric, planar and rotated families, on 800 random graph-like codes against a minimum-cycle-basis check (706
  of them also against exhaustive enumeration), and on all 456 board entries with a graph-like side. On the whole census the two
  programs agree on all 7063 codes.
- **tropidist** (B. Cheng, N. Guo, F. Pan, *Certifying quantum code and circuit distances at scale*,
  arXiv:2611.08247; github.com/Fanerst/tropidist, MIT license, not included here), at commit `b03ec2d4077f`, the
  commit of the board's tropidist certificates, one thread, 600 s per run, on the problems written by
  `src/tropidist_inputs.py` (the board's "css-sides" encoding): exact on all 18 sides, the proved lower bound equals
  the weight of the vector found. Records in `logs/tropidist/`; the board certificates `certs/<slug>.json` are
  these runs. We rechecked every vector the tool returned (in the kernel, outside the other row space, weight $d$).

| code | tropidist side X: bounds, time | side Z | graphlike_A / graphlike_B (census run) |
|---|---|---|---|
| `864-434-6` | 6 = 6, 2.3 s | 6 = 6, 2.2 s | 0.008 s / 0.012 s |
| `720-362-6` | 6 = 6, 1.6 s | 6 = 6, 1.6 s | 0.006 s / 0.009 s |
| `936-470-6` | 6 = 6, 2.6 s | 6 = 6, 2.7 s | 0.009 s / 0.013 s |
| `648-326-6` | 6 = 6, 1.1 s | 6 = 6, 1.3 s | 0.005 s / 0.007 s |
| `640-322-6` | 6 = 6, 1.4 s | 6 = 6, 1.3 s | 0.004 s / 0.006 s |
| `1000-502-6` | 6 = 6, 2.9 s | 6 = 6, 2.7 s | 0.010 s / 0.016 s |
| `294-100-6` | 6 = 6, 4.2 s | 6 = 6, 0.5 s | 0.001 s / 0.002 s |
| `486-164-6` | 6 = 6, 0.9 s | 6 = 6, 1.0 s | 0.002 s / 0.004 s |
| `720-74-8` | 8 = 8, 1.9 s | 8 = 8, 1.8 s | 0.003 s / 0.009 s |

The result is a computer-assisted proof. Four independently written programs agree on every code: the three of
ours, whose methods are proved correct in their notes, and tropidist, whose completeness is proved in its paper.
This makes an undetected programming error unlikely but does not exclude it.

## Reproduce

Requirements: Python 3.8 or later (standard library only), a C compiler, `make`.

```
make check             # src/rebuild.py: matrices from the stored darts and from the presentations, Conder's
                       # relators, k, commutation, weights, witnesses (seconds)
make distances         # builds bin/graphlike_A, bin/graphlike_B; runs both and src/circuits_check.py (about a
                       # minute on one core)
make tropidist-inputs  # the tropidist problems in tropidist-inputs/
```

To rerun tropidist: clone github.com/Fanerst/tropidist, `git checkout b03ec2d4077f287f89f97b5c2c16ce5ce17c2c58`, install its
requirements (Python, numpy, stim, a C compiler), and from its root run
`python -m cdist --matrices H.mtx L.mtx --threads 1 --time-limit 600 --out rec.json` for each pair
`tropidist-inputs/<slug>.<side>.H.mtx`, `.L.mtx`. `logs/tropidist/inputs.sha256` lists the checksums of the
inputs of our runs, which `make tropidist-inputs` reproduces byte for byte.

## Layout

- `data/codes/<slug>.hx`, `.hz`: $H_X$, $H_Z$ (first line `n m`, then one line per check with its 0-based qubits).
- `data/maps.json`: dart permutations `R`, `S` (0-based, right action; dart 0 is the identity) of each map; for
  `720-74-8` the base map, the voltage on its edges (numbered by smallest dart) and `m` = 2. The matrices are the
  code of these darts with qubit $e$ = the $e$-th edge $\{x,x^{RS}\}$ in the order of its smallest dart,
  X-check $i$ = the $i$-th cycle of `S`, Z-check $j$ = the $j$-th cycle of `R`.
- `data/presentations.json`: Conder's names and verbatim relators, the presentations used above.
- `data/cover_voltage.txt`: the 41 words of the voltage of `720-74-8`.
- `data/parameters.json`: $n,k,d_X,d_Z$, witnesses, numbers of minimum-weight logical operators.
- `src/`: `rebuild.py`, `run_distances.py`, `circuits_check.py`, `tropidist_inputs.py`, `graphlike_A/`,
  `graphlike_B/` (C sources and the notes with the proofs).
- `logs/`: `rebuild.log` and `distances.jsonl` (outputs of `make check` and `make distances` in a clean copy of
  this directory, `environment.txt` the versions used), `tropidist/` (records of our tropidist runs, `runs.tsv`
  with wall times, `inputs.sha256`).
- `SHA256SUMS`: checksums of all files.

## Credits and license

Hyperbolic surface codes: Breuckmann and Terhal (arXiv:1506.04029); the maps and their names: M. Conder's lists
(relators quoted verbatim in `data/presentations.json`); the census was computed with Magma
(`LowIndexNormalSubgroups`); exact distances of codes: tropidist (Cheng, Guo and Pan, arXiv:2611.08247). MIT license
(see `../LICENSE`), except for Conder's relators, which are quoted from his lists.
