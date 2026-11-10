# Surface codes from regular maps: checks and exact distances (see README.md)
CC ?= cc
CFLAGS ?= -O3 -march=native
PYTHON ?= python3

all: check distances

bin/graphlike_A: src/graphlike_A/graphlike_A.c
	mkdir -p bin
	$(CC) $(CFLAGS) -std=gnu11 -Wall -Wextra -o $@ $<

bin/graphlike_B: src/graphlike_B/graphlike_B.c
	mkdir -p bin
	$(CC) $(CFLAGS) -std=c11 -Wall -Wextra -o $@ $<

check:
	$(PYTHON) src/rebuild.py

distances: bin/graphlike_A bin/graphlike_B
	$(PYTHON) src/run_distances.py

tropidist-inputs:
	$(PYTHON) src/tropidist_inputs.py tropidist-inputs

sums:
	sha256sum -c SHA256SUMS

clean:
	rm -rf bin tropidist-inputs

.PHONY: all check distances tropidist-inputs sums clean
