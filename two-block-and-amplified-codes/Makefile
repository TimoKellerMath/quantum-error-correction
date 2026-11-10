# Rebuild and check the three codes, and rerun the exact searches of README.md. Run from this directory.
PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3 -march=native
MAGMA ?= magma
BIN := bin/exactA bin/exactB
SAME := $(PYTHON) src/same_results.py
I := instances/amplified
# base runs: BASE:SIDE:W (logs/amplified/cert/), direct runs on the amplified codes: BASE:SIDE:W (logs/amplified/amp/)
CERT := 185-65-5:X:4 185-65-5:X:5 185-65-5:Z:4 185-65-5:Z:5 96-17-6:X:5 96-17-6:X:6 96-17-6:Z:5 96-17-6:Z:6
AMP := 185-65-5:X:9 185-65-5:X:10 185-65-5:Z:9 185-65-5:Z:10 96-17-6:X:11 96-17-6:X:12 96-17-6:Z:11 96-17-6:Z:12

.PHONY: check bin quick long long-first tropidist-inputs blocks clean

check:
	$(PYTHON) src/rebuild.py

bin: $(BIN)

bin/exactA: src/exactA.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $< -lpthread

bin/exactB: src/exactB.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -pthread -o $@ $<

# the 32 searches on the bases 185-65-5, 96-17-6 and on 860-130-10, 463-34-12 (a few seconds)
quick: $(BIN)
	@mkdir -p rerun/amplified/cert rerun/amplified/amp
	@for j in $(CERT); do b=$${j%%:*}; r=$${j#*:}; s=$${r%%:*}; w=$${r#*:}; \
	  ./bin/exactA $(I)/A/$${b}_r1.hx $(I)/A/$${b}_r1.hz --side $$s --W $$w --group $(I)/A/$${b}_r1.gens \
	    --mode anchor > rerun/amplified/cert/A_$${b}_$${s}_W$$w.log 2>&1 \
	  && $(SAME) rerun/amplified/cert/A_$${b}_$${s}_W$$w.log logs/amplified/cert/A_$${b}_$${s}_W$$w.log \
	  && ./bin/exactB $(I)/B/$$b.hx $(I)/B/$$b.hz --side $$s --W $$w --group $(I)/B/$$b.gens \
	    > rerun/amplified/cert/B_$${b}_$${s}_W$$w.log 2>&1 \
	  && $(SAME) rerun/amplified/cert/B_$${b}_$${s}_W$$w.log logs/amplified/cert/B_$${b}_$${s}_W$$w.log \
	  || exit 1; done
	@for j in $(AMP); do b=$${j%%:*}; r=$${j#*:}; s=$${r%%:*}; w=$${r#*:}; t=amp422_$$b; \
	  ./bin/exactA $(I)/$$t.hx $(I)/$$t.hz --side $$s --W $$w --group $(I)/A/$$t.gens --mode anchor \
	    > rerun/amplified/amp/A_$${t}_$${s}_W$$w.log 2>&1 \
	  && $(SAME) rerun/amplified/amp/A_$${t}_$${s}_W$$w.log logs/amplified/amp/A_$${t}_$${s}_W$$w.log \
	  && ./bin/exactB $(I)/$$t.hx $(I)/$$t.hz --side $$s --W $$w --group $(I)/B/$$t.gens \
	    > rerun/amplified/amp/B_$${t}_$${s}_W$$w.log 2>&1 \
	  && $(SAME) rerun/amplified/amp/B_$${t}_$${s}_W$$w.log logs/amplified/amp/B_$${t}_$${s}_W$$w.log \
	  || exit 1; done
	@echo "quick: all 32 reruns agree with the logs"

# the second set of searches on 480-16-26 (about 11 CPU-minutes)
long: $(BIN)
	@mkdir -p rerun/480-16-26/second
	@for s in Z X; do \
	  ./bin/exactA instances/c480v.hx instances/c480v.hz --side $$s --W 24 --group instances/c480v.gensA \
	    --mode pattern > rerun/480-16-26/second/exactA_pattern_$$s.log 2>&1 \
	  && $(SAME) rerun/480-16-26/second/exactA_pattern_$$s.log logs/480-16-26/second/exactA_pattern_$$s.log \
	  || exit 1; done
	@for s in Z X; do \
	  ./bin/exactB instances/c480v.hx instances/c480v.hz --side $$s --W 24 --group instances/c480v.gensB \
	    > rerun/480-16-26/second/exactB_orbits_$$s.log 2>&1 \
	  && $(SAME) rerun/480-16-26/second/exactB_orbits_$$s.log logs/480-16-26/second/exactB_orbits_$$s.log \
	  || exit 1; done
	@echo "long: all 4 reruns agree with the logs"

# the first set of searches on 480-16-26 (about 25 CPU-minutes)
long-first: $(BIN)
	@mkdir -p rerun/480-16-26/first
	@for s in Z X; do \
	  ./bin/exactA instances/c480_16.hx instances/c480_16.hz --side $$s --W 25 --group instances/c480_16.gens \
	    > rerun/480-16-26/first/exactA_$$s.log 2>&1 \
	  && $(SAME) rerun/480-16-26/first/exactA_$$s.log logs/480-16-26/first/exactA_$$s.log \
	  || exit 1; done
	@for s in Z X; do \
	  ./bin/exactB instances/c480_16.hx instances/c480_16.hz --side $$s --W 25 --group instances/c480_16.gensB \
	    > rerun/480-16-26/first/exactB_$$s.log 2>&1 \
	  && $(SAME) rerun/480-16-26/first/exactB_$$s.log logs/480-16-26/first/exactB_$$s.log \
	  || exit 1; done
	@echo "long-first: all 4 reruns agree with the logs"

tropidist-inputs:
	$(PYTHON) src/tropidist_inputs.py write tropidist-inputs 480-16-26.X 480-16-26.Z 860-130-10.X 860-130-10.Z \
	    463-34-12.X 463-34-12.Z
	cd tropidist-inputs && sha256sum -c ../logs/tropidist/inputs.sha256

# optional, needs Magma
blocks:
	$(MAGMA) -b src/blocks.m

clean:
	rm -rf bin rerun tropidist-inputs
