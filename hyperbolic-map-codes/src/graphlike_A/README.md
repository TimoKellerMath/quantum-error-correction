# graphlike_A

`graphlike_A.c` and `NOTES.md` are copied unchanged from our working directory (2026-11-10). `NOTES.md` proves that
the program computes $d_Z$ exactly when every column of $H_X$ has weight at most 2, and $d_X$ when every column of
$H_Z$ has (Lemmas 1 to 3, Proposition 5, Lemma 7), and that `--count` counts the logical operators of minimum weight
(Proposition 8). The test drivers and the results file it mentions are not part of this package; the validation is
summarized in `../../README.md`.

Build from the package root with `make bin/graphlike_A`. Usage:
`bin/graphlike_A NAME.hx NAME.hz [--side X|Z|both] [--quiet] [--count | --count-limit NODES]`; output one JSON line.
