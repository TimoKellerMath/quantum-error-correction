ifeq ($(origin CC),default)
CC = gcc
endif
CFLAGS ?= -O2
P ?= 4

.PHONY: all check validate search pdf clean

all: check validate search

bin/vsearch_orb: src/vsearch_orb.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $<

bin/cluster2: src/cluster2.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $<

# builds the codes, checks n, k, weights, automorphisms and witnesses, writes instances/
check:
	python3 src/codes.py

validate: bin/vsearch_orb bin/cluster2 check
	sh src/validate.sh

# all lower-bound searches of exact_distances.tex, Section 4 (logs in rerun/)
search: bin/vsearch_orb bin/cluster2 check
	sh src/run_all.sh $(P)

pdf:
	pdflatex exact_distances.tex && pdflatex exact_distances.tex

clean:
	rm -rf bin rerun *.aux *.log *.out
