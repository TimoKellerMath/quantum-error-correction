# Exact distances of five weight-six two-block codes

The note [`exact_distances.pdf`](exact_distances.pdf) (source `exact_distances.tex`) determines the minimum
distance of five quantum CSS codes with checks of weight six and qubit degree six:

| name | group | $a$, $b$ | $[[n,k,d]]$ | $kd^2/n$ |
|---|---|---|---|---|
| C600 | $\langle x,y\mid x^{30},y^{10},yxy^{-1}=x^{19}\rangle\cong D_5\times C_{30}$ | $1+x^{24}y^8+x^{13}y^9$, $1+x^{22}y^4+x^{24}y^7$ | $[[600,24,26]]$ | 27.04 |
| C504 | $\mathbb Z/252$ | $1+x+x^{58}$, $1+x^{48}+x^{166}$ | $[[504,12,28]]$ | 18.67 |
| B510 | $\mathbb Z/255$ | $1+x^{26}+x^{76}$, $1+x^{30}+x^{66}$ | $[[510,16,26]]$ | 21.21 |
| B558 | $\mathbb Z/279$ | $1+x^{37}+x^{140}$, $1+x^{61}+x^{197}$ | $[[558,14,28]]$ | 19.67 |
| B434 | $\mathbb Z/217$ | $1+x^{29}+x^{117}$, $1+x^{19}+x^{111}$ | $[[434,10,28]]$ | 18.06 |

C600 is a two-block group algebra code $C(G;a,b)$ with $H_X=(L_a\mid R_b)$, $H_Z=(R_b^{\mathsf T}\mid L_a^{\mathsf T})$.
The other four are cyclic two-block codes in the convention of the
[qLDPC challenge board](https://github.com/unitaryfoundation/qldpc-challenge):
$H_X=(\mathrm{circ}(a)\mid\mathrm{circ}(b))$, $H_Z=(\mathrm{circ}(b)^{\mathsf T}\mid\mathrm{circ}(a)^{\mathsf T})$.
B510, B558 and B434 are the board entries `510-16-26`, `558-14-28` and `434-10-28` (submitted by V. Russo with
witness-backed upper bounds); `src/codes.py --board DIR` checks that our matrices coincide with the board's files.

## What is proved, and how

- **Upper bounds:** explicit logical operators of weight $d$ in `data/*_witnesses_{X,Z}.json`
  (for C504 the witness found at $W=28$ is in `runs/found_C504_W28.txt`).
- **Lower bounds (computer-assisted):** an exhaustive anchored search with exclusions (Algorithm 1 of the note,
  completeness: Proposition 3.1 and Corollary 3.2) reports no nontrivial logical operator of weight at most $d-2$.
  All column weights are 3, so logical operators have even weight (Lemma 2.2) and this gives $d$.
  For the cyclic codes an explicit qubit permutation exchanges $X$ and $Z$ (Lemma 2.3); for C600 both sides
  were searched, with anchors from a group of 150 code automorphisms (`data/C600_automorphisms.json`).
- **Two independent implementations** of the search, which share no code: `src/vsearch_orb.c` (A) and
  `src/cluster2.c` (B). Both were validated on the gross code $[[144,12,12]]$ and the two-gross code
  $[[288,12,18]]$ (`runs/validation.log`).

The result is a computer-assisted proof. Two independently written programs, which follow the same argument (proved in
the note), agree on every run; this makes an undetected programming error unlikely but does not exclude it.

## Reproduce

Requirements: a C compiler, Python 3 (standard library only), `make`.

```
make check      # build the codes, check n, k, weights, automorphisms and witnesses, write instances/
make validate   # both programs on the gross and two-gross codes (seconds)
make search     # all searches of Section 4 of the note, logs in rerun/ (set P=<parallel jobs>)
```

`src/jobs.txt` lists every search with its command. The logs of our runs are in `runs/`. The searches were run in a working
copy in which `instances/` was called `codes/`, hence the paths `codes/...` in their logs; the instance files there were
byte-identical to those in `instances/`, and the binaries byte-identical to those that `make` builds from `src/` with
gcc 16.2.1 (the sources differed only in comments). `SHA256SUMS` lists the checksums of the sources, data, instances
and logs. Times in the note are single-thread wall-clock seconds
(gcc 16.2.1, `-O2`, x86-64), measured with up to seven jobs running in parallel; a job run alone is up to about twice as fast.

## Layout

- `exact_distances.tex`, `exact_distances.pdf`: the note.
- `src/`: `codes.py` (constructions, checks, instance writer), `vsearch_orb.c`, `cluster2.c`, `jobs.txt`,
  `run_all.sh`, `validate.sh`.
- `data/`: check matrices (`<name>.json`), witnesses, automorphism generators of C600.
- `instances/`: the inputs of the searches, as written by `src/codes.py`.
- `runs/`: logs of the certification runs, of the validation, and of `python3 src/codes.py --board DIR` with `DIR` the
  `codes/` directory of the board repository (`codes_check.log`).

## Credits and license

Two-block group algebra codes are due to Wang, Lin and Pryadko (arXiv:2305.06890) and Lin and Pryadko
(arXiv:2306.16400), the gross and two-gross codes to Bravyi et al. (arXiv:2308.07915), and B510, B558 and B434 with
their witnesses to V. Russo (@vprusso). The search refines known complete searches: Dumer, Kovalev and Pryadko
(arXiv:1611.07164), Okada and Kasai (arXiv:2607.14091), Lee, Okada, Maskara, Kasai and Zhou (arXiv:2609.30069),
J. Y. Lee (arXiv:2610.06820), Cheng, Guo and Pan (arXiv:2611.08247), Rowshan and Devitt (arXiv:2609.39395), and the
fieldnote `fieldnotes/2026-10-13-bb360-12-exhaustive-d24.md` of C. Crocker in the qLDPC challenge repository; see
Section 3 of the note for what is taken from where. The programs in `src/`, the computations and a first draft of the note were made with
Claude Opus 5.5 (Anthropic) in Claude Code. MIT license (see `../LICENSE`), except for `data/B*_witnesses_{X,Z}.json`, which are copied from
V. Russo's board entries in the qLDPC challenge repository (Apache License 2.0,
https://www.apache.org/licenses/LICENSE-2.0).
