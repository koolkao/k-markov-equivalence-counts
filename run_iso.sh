#!/usr/bin/env bash
# Exact counts via isomorphism classes: ./run_iso.sh n [kmin] [kmax] [threads]
# Needs nauty's geng and directg on PATH (tested with nauty 2.9.3).
# Fails if the pipe breaks, if the number of unlabeled DAGs read differs from A003087(n),
# or if kmec_iso reports any warning.
set -euo pipefail
n=$1; kmin=${2:-0}; kmax=${3:-$((n-2))}; threads=${4:-14}
A003087=(1 1 2 6 31 302 5984 243668 20286025)
out=results_iso_n${n}.txt; log=iso_n${n}.log
geng -q "$n" | directg -q -a -T | ./kmec_iso "$n" "$kmin" "$kmax" "$threads" > "$out" 2> "$log"
read_count=$(sed -n 's/^read \([0-9]*\) unlabeled DAGs.*/\1/p' "$log")
if [ "$read_count" != "${A003087[$n]}" ]; then echo "ERROR: read $read_count DAGs, expected A003087($n) = ${A003087[$n]}" >&2; exit 1; fi
if grep -q WARNING "$log"; then echo "ERROR: warnings in $log" >&2; exit 1; fi
echo "OK: $read_count unlabeled DAGs; results in $out"
cat "$out"
