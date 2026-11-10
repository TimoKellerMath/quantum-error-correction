# Build the four programs; BOARD is a checkout of github.com/unitaryfoundation/qldpc-challenge.
CC ?= gcc
CFLAGS ?= -O3 -march=native -Wall -Wextra
BOARD ?= ../qldpc-challenge

.PHONY: all check check-long clean
all: bin/graphlike_A bin/graphlike_B bin/exactA bin/exactB

bin/graphlike_A: src/graphlike_A.c
	mkdir -p bin
	$(CC) $(CFLAGS) -std=gnu11 -o $@ $<

bin/graphlike_B: src/graphlike_B.c
	mkdir -p bin
	$(CC) $(CFLAGS) -std=c11 -o $@ $<

bin/exactA: src/exactA.c
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $< -lpthread

bin/exactB: src/exactB.c
	mkdir -p bin
	$(CC) $(CFLAGS) -pthread -o $@ $<

# every entry except the four long partial searches (a few minutes on one core)
check: all
	python3 src/certify.py --board $(BOARD)

# including the long partial searches (about one CPU-hour more)
check-long: all
	python3 src/certify.py --board $(BOARD) --long

clean:
	rm -rf bin
