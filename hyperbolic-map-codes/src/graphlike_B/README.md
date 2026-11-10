# graphlike_B

`graphlike_B.c` and `NOTES.md` are copied unchanged from our working directory (2026-11-10); the program was written
independently of graphlike_A and shares no code with it. `NOTES.md` proves exactness (Lemmas 1 to 3, Corollary 4,
Propositions 8 and 9, Lemma 10) and the principle of the circuit enumeration `../circuits_check.py` (Lemma 12). The
test drivers and the results file it mentions are not part of this package; the validation is summarized in
`../../README.md`.

Build from the package root with `make bin/graphlike_B`. Usage: `bin/graphlike_B [--ref] NAME.hx NAME.hz`; output
one JSON line, `--ref` adds the unpruned reference search.
