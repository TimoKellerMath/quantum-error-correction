# Quantum error correction

Work of Timo Keller on quantum error-correcting codes. Programs, computations and first drafts of the texts were made
with Claude Opus 5.5 (Anthropic) in Claude Code; the README of each directory says what.

- [`weight6-exact-distances/`](weight6-exact-distances/): exact minimum distances of five quantum LDPC
  codes with checks of weight six (a non-abelian two-block group algebra code `[[600,24,26]]` and four
  cyclic two-block codes, three of them entries of V. Russo on the qLDPC challenge board), with a short note containing the proofs, the programs of the computer-assisted
  part, their inputs and logs.
- [`hyperbolic-map-codes/`](hyperbolic-map-codes/): nine quantum surface codes, eight from regular maps of type
  {8,8} and {6,6} in M. Conder's lists and one from a double cover of a {4,5} map, with exact minimum distances
  (three programs of ours with proofs of their methods, and tropidist), submitted to the qLDPC challenge board;
  data, programs and logs.
- [`amplified-and-folded-codes/`](amplified-and-folded-codes/): eight distance-amplified CSS codes whose
  distance is exactly twice the certified distance of their bases (with the proof), and two non-CSS stabilizer
  codes obtained by halving two-block codes along ZX-dualities, submitted to the qLDPC challenge board; data, a
  rebuild script, search programs and logs.
- [`board-certificates-2026-11/`](board-certificates-2026-11/): exact minimum distances (and three intervals) of
  163 entries of other authors on the qLDPC challenge board, and four corrected side values: four programs with
  proofs of their methods, the distance argument for amplified codes, a driver that reruns everything from the
  board's files, data and logs.
- [`two-block-and-amplified-codes/`](two-block-and-amplified-codes/): a weight-six two-block group algebra code
  `[[480,16,26]]` over C12 x Dic5 with exact distance (two exhaustive search programs with the proofs of their
  completeness) and two distance-amplified codes `[[860,130,10]]` and `[[463,34,12]]` with exact distances,
  submitted to the qLDPC challenge board; data, programs, a rebuild script and logs.
- [`folded-stabilizer-codes/`](folded-stabilizer-codes/): five non-CSS stabilizer codes, folds of two-block
  entries of the qLDPC challenge board along ZX-dualities (mirror codes in the sense of Khesin and Lu), with
  exact distances (tropidist, a meet-in-the-middle search and Magma), submitted to the board; data, a rebuild
  script, the search program and logs.
- [`pair-partition-and-two-block-codes/`](pair-partition-and-two-block-codes/): five quasi-cyclic CSS codes with
  checks of weight six, lifted from sparse 6 x 14 pair-partition templates, and a weight-six two-block group
  algebra code `[[672,28,24]]` over C28 x Dic3, with exact distances (two exhaustive search programs of ours on
  both sides; dist-m4ri and, for the five pair-partition codes, tropidist agree), submitted to the qLDPC
  challenge board; definitions, data, a rebuild script, the search records and logs.

License: MIT (see `LICENSE`), except for files copied from other sources; the README of their directory names them
and their license.
