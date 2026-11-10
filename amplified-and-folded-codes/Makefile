# Rebuild and check the ten codes, and rerun the searches of Part B of README.md. Run from this directory.
PYTHON ?= python3
CC ?= cc
CFLAGS ?= -O3
BIN := bin/mitm bin/fewnony2 bin/count_odd bin/isd
SAME := $(PYTHON) src/same_results.py

.PHONY: check bin quick long tropidist-inputs clean

check:
	$(PYTHON) src/rebuild.py

bin: $(BIN)

bin/%: src/%.c
	@mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $<

# about 3 CPU-minutes
quick: $(BIN)
	@mkdir -p rerun
	./bin/mitm instances/gross_c4.code 4 16 > rerun/mitm_gross_c4.log 2> rerun/mitm_gross_c4.err
	$(SAME) rerun/mitm_gross_c4.log logs/folds/mitm_gross_c4.log
	./bin/mitm instances/test513.code 2 4 > rerun/mitm_test513.log 2> /dev/null
	$(SAME) rerun/mitm_test513.log logs/folds/mitm_validation.log
	./bin/mitm instances/test_13-1-5.code 3 4 > rerun/mitm_13-1-5.log 2> /dev/null
	$(SAME) rerun/mitm_13-1-5.log logs/folds/mitm_validation.log
	./bin/mitm instances/test_17-1-7.code 4 4 > rerun/mitm_17-1-7.log 2> /dev/null
	$(SAME) rerun/mitm_17-1-7.log logs/folds/mitm_validation.log
	./bin/mitm instances/test_21-1-7.code 4 4 > rerun/mitm_21-1-7.log 2> /dev/null
	$(SAME) rerun/mitm_21-1-7.log logs/folds/mitm_validation.log
	./bin/fewnony2 instances/gross_c4.code 4 > rerun/fewnony2_gross_c4.log
	$(SAME) rerun/fewnony2_gross_c4.log logs/folds/fewnony2_validation.log
	./bin/fewnony2 instances/B510_c0.code 4 5 > rerun/fewnony2_B510_c0.log
	$(SAME) rerun/fewnony2_B510_c0.log logs/folds/fewnony2_B510_t0_T4_A5.log
	for s in 85 15 5; do $(PYTHON) src/orbit_fold.py instances/B510_c0.json $$s rerun/B510_c0_orb$$s.code \
	    > rerun/orbit_fold_$$s.log && cmp rerun/B510_c0_orb$$s.code instances/B510_c0_orb$$s.code \
	    && $(SAME) rerun/orbit_fold_$$s.log logs/folds/orbit_fold_B510.log || exit 1; done
	./bin/mitm instances/B510_c0_orb5.code 3 4 > rerun/mitm_orb5.log 2> /dev/null
	$(SAME) rerun/mitm_orb5.log logs/folds/mitm_orbitfold_B510.log
	./bin/mitm instances/B510_c0_orb15.code 7 16 > rerun/mitm_orb15.log 2> /dev/null
	$(SAME) rerun/mitm_orb15.log logs/folds/mitm_orbitfold_B510.log
	./bin/isd instances/gross_c4.code 100000 20261110 0.3 rerun/gross_k0 > rerun/isd_gross_k0.log
	$(SAME) rerun/isd_gross_k0.log logs/folds/isd_gross_k0.log
	@echo "quick: all reruns agree with the logs"

# about 1 CPU-hour
long: $(BIN)
	@mkdir -p rerun
	./bin/mitm instances/board_72-6-8-b.code 4 16 > rerun/mitm_board_72-6-8-b.log 2> /dev/null
	$(SAME) rerun/mitm_board_72-6-8-b.log logs/folds/mitm_board_72-6-8-b.log
	for i in 1 2 3; do ./bin/mitm instances/gross_c$$i.code 4 16 > rerun/mitm_gross_c$$i.log 2> /dev/null \
	    && $(SAME) rerun/mitm_gross_c$$i.log logs/folds/mitm_gross_c$$i.log || exit 1; done
	./bin/mitm instances/B510_c0_orb85.code 4 32 > rerun/mitm_orb85.log 2> /dev/null
	$(SAME) rerun/mitm_orb85.log logs/folds/mitm_orbitfold_B510.log
	./bin/count_odd instances/gross_c4.code 4 16 100 > rerun/count9_gross_k0.log 2> /dev/null
	$(SAME) rerun/count9_gross_k0.log logs/folds/count9_gross_k0.log
	./bin/isd instances/B510_c0.code 1000000 101 0.0 rerun/B510_t0_Q > rerun/isd_B510_Q.log
	$(SAME) rerun/isd_B510_Q.log logs/folds/isd_B510_t0.log
	./bin/isd instances/B510_c0.code 1000000 202 1.0 rerun/B510_t0_B > rerun/isd_B510_B.log
	$(SAME) rerun/isd_B510_B.log logs/folds/isd_B510_t0.log
	./bin/isd instances/B510_c0.code 1000000 303 0.0 rerun/B510_t0_Q2 > rerun/isd_B510_Q2.log
	$(SAME) rerun/isd_B510_Q2.log logs/folds/isd_B510_t0_Q2.log
	@echo "long: all reruns agree with the logs"

tropidist-inputs:
	$(PYTHON) src/tropidist_inputs.py write tropidist-inputs 72-6-9-s.P 866-128-12.X 866-128-12.Z 610-40-22.X
	cd tropidist-inputs && sha256sum -c ../logs/tropidist/inputs.sha256

clean:
	rm -rf bin rerun tropidist-inputs
