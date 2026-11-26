# Reproduce every result file. Tested on macOS 26 (Apple clang 21, Apple M4 Pro, 48 GB RAM).
SHELL := /bin/bash
CC ?= cc
CFLAGS ?= -O3 -march=native -pthread
THREADS ?= 14

all: kmec kmec2 kmec_unl kmec_iso tests/test_bayesball

kmec: kmec.c
	$(CC) $(CFLAGS) -o $@ $<
kmec2: kmec2.c
	$(CC) $(CFLAGS) -o $@ $<
kmec_unl: kmec_unl.c
	$(CC) $(CFLAGS) -o $@ $<
kmec_iso: kmec_iso.c
	$(CC) $(CFLAGS) -o $@ $<
tests/test_bayesball: tests/test_bayesball.c kmec2.c
	$(CC) -O2 -pthread -o $@ $<

# Labeled table, n = 2..6, exact signatures (about 1 minute).
labeled: kmec
	for n in 2 3 4 5 6; do ./kmec $$n; done > results_labeled.txt

# Unlabeled table, n = 2..6, exact signatures (about 3 minutes; n = 6 needs about 2 GB RAM).
# Output columns: n k labeled unlabeled.
unlabeled: kmec_unl
	for n in 2 3 4 5 6; do ./kmec_unl $$n; done > results_unlabeled.txt

# Labeled row n = 7 with hashed signatures. Each run takes about 35-50 minutes on 14 cores
# and about 30 GB RAM. The two seeds must agree.
n7: kmec2
	./kmec2 7 5 $(THREADS) 0 > results_n7_seed0.txt
	./kmec2 7 5 $(THREADS) 7 > results_n7_seed7.txt
	diff results_n7_seed0.txt results_n7_seed7.txt && echo "n = 7: both seeds agree"

# Quick checks (about 2 minutes).
check: kmec kmec2 tests/test_bayesball
	./tests/test_bayesball
	for n in 2 3 4 5 6; do ./kmec2 $$n $$((n-2)) $(THREADS) 3 2>/dev/null; done | diff - results_labeled.txt && echo "kmec2 (seed 3) matches results_labeled.txt"
	python3 kmec_ref.py 5 2>/dev/null | diff - <(grep '^5 ' results_labeled.txt) && echo "networkx reference matches row n = 5"

# Exact counts via isomorphism classes (needs nauty's geng and directg on PATH).
# n = 2..7 take seconds; n = 8 takes about 28 minutes and about 31 GB (see README).
exact: kmec_iso
	for n in 2 3 4 5 6 7; do ./run_iso.sh $$n 0 $$((n-2)) $(THREADS) > /dev/null; done
	./run_iso.sh 8 0 6 $(THREADS) > /dev/null

.PHONY: all labeled unlabeled n7 exact check
