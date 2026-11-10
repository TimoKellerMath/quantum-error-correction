# Rebuild and check the five codes, and rerun the exhaustive searches of README.md. Run from this directory.
PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3 -march=native
MAGMA ?= magma
GAP ?= gap
SAME := $(PYTHON) src/same_results.py
CODES := 60-9-8 56-6-8-s 32-9-6 48-6-8-s 33-2-7
# meet-in-the-middle runs SLUG:a:b (lower bound a + b = d - 1, and the positive controls a + b = d)
MITM := 60-9-8:4:3 56-6-8-s:4:3 32-9-6:3:2 32-9-6:3:3 48-6-8-s:4:3 48-6-8-s:4:4 33-2-7:3:3 33-2-7:4:3

.PHONY: check quick magma gap tropidist-inputs clean

check:
	$(PYTHON) src/rebuild.py

bin/mitm: src/mitm.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $<

# the eight meet-in-the-middle runs (about 10 CPU-seconds, 0.4 GB of memory for 48-6-8-s:4:4)
quick: bin/mitm
	@mkdir -p rerun/mitm
	@for j in $(MITM); do s=$${j%%:*}; r=$${j#*:}; a=$${r%%:*}; b=$${r#*:}; \
	  ./bin/mitm instances/$$s.txt $$a $$b > rerun/mitm/$${s}_a$${a}_b$${b}.log 2>&1 \
	  && $(SAME) rerun/mitm/$${s}_a$${a}_b$${b}.log logs/mitm/$${s}_a$${a}_b$${b}.log || exit 1; done
	@echo "quick: all 8 reruns agree with the logs"

# Magma's MinimumDistance on the five codes (about 90 CPU-seconds, most of it for 60-9-8)
magma:
	@mkdir -p rerun/magma
	@for s in $(CODES); do $(MAGMA) -b magma/$$s.m > rerun/magma/$$s.log 2>&1 \
	  && $(SAME) rerun/magma/$$s.log logs/magma/$$s.log || exit 1; done
	@echo "magma: all 5 reruns agree with the logs"

# GAP's numbering of the groups of the four parents in GAP's SmallGroups library
gap:
	@mkdir -p rerun
	$(GAP) -q src/gap_elements.g | tr -d '\r' > rerun/gap_elements.log
	cmp rerun/gap_elements.log logs/gap_elements.log && echo "gap: same output as logs/gap_elements.log"

tropidist-inputs:
	$(PYTHON) src/tropidist_inputs.py write tropidist-inputs $(addsuffix .P,$(CODES))
	cd tropidist-inputs && sha256sum -c ../logs/tropidist/inputs.sha256

clean:
	rm -rf bin rerun tropidist-inputs
