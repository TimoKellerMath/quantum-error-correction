# Check the six codes (make test) and rerun the 24 exact searches (make rerun). Run from this directory.
PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3 -march=native
# exactA.c and exactB.c: the files of two-block-and-amplified-codes/src/ in this repository (identical copies are in
# board-certificates-2026-11/src/); their sha256 must be the ones in src/programs.sha256
PROGS ?= ../two-block-and-amplified-codes/src

.PHONY: test sums check related bin rerun clean

# the checksums of every file, then src/rebuild.py (about 10 s)
test: sums check

sums:
	@sha256sum -c --quiet SHA256SUMS && echo "SHA256SUMS: every listed file matches"

check:
	$(PYTHON) src/rebuild.py

# 672-28-24/related/k_covers.py (README, Part B, Lemma 8), compared with our output (a few seconds)
related:
	$(PYTHON) 672-28-24/related/k_covers.py | cmp - 672-28-24/related/k_covers.out && echo "k_covers: output equals k_covers.out"

bin: bin/exactA bin/exactB

bin/exactA: $(PROGS)/exactA.c src/programs.sha256
	cd $(PROGS) && sha256sum -c $(CURDIR)/src/programs.sha256
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $(PROGS)/exactA.c -lpthread

bin/exactB: $(PROGS)/exactB.c src/programs.sha256
	cd $(PROGS) && sha256sum -c $(CURDIR)/src/programs.sha256
	@mkdir -p bin
	$(CC) $(CFLAGS) -pthread -o $@ $(PROGS)/exactB.c

# the 24 searches of records/ (about 45 CPU-seconds, most of it on 672-28-24), compared with our logs
rerun: bin
	$(PYTHON) src/rerun.py

clean:
	rm -rf bin rerun
