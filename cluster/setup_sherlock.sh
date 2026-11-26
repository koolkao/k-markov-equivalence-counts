#!/usr/bin/env bash
# One-time setup on a Sherlock login node: build nauty (geng, directg) and kmec_cluster, then run
# quick checks. Usage, from the repository root:
#   bash cluster/setup_sherlock.sh
# Puts binaries in $HOME/kmec-bin and creates logs/ (Slurm does not create log directories).
# Compute nodes mix Intel and AMD CPUs, so everything is built for x86-64-v2, not -march=native.
# x86-64-v2 needs GCC >= 11. If the system gcc is older, run `module load gcc` first.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$HOME/kmec-bin; SRC=$HOME/kmec-src
mkdir -p "$BIN" "$SRC" "$REPO/logs"
gcc --version | head -1
if [ ! -x "$BIN/directg" ]; then
  ( cd "$SRC"
    curl -fL -o nauty2_9_3.tar.gz https://pallini.di.uniroma1.it/nauty2_9_3.tar.gz
    tar xzf nauty2_9_3.tar.gz
    cd nauty2_9_3
    ./configure --enable-generic CFLAGS="-O3 -march=x86-64-v2" > configure.log 2>&1
    make -j8 geng directg > make.log 2>&1
    cp geng directg "$BIN"/ )
fi
gcc -O3 -march=x86-64-v2 -pthread -Wall -o "$BIN/kmec_cluster" "$REPO/kmec_cluster.c"
export PATH=$BIN:$PATH
# Checks: nauty's unlabeled DAG counts (A003087) and exact counts for n = 6, 7.
for n in 5 6 7; do printf "n=%s: " "$n"; geng -q "$n" | directg -a -u 2>&1 | grep -o '[0-9]* digraphs generated'; done
for n in 6 7; do
  kmec_cluster direct "$n" 0 $((n-2)) 4 16 | grep -v '^#' | awk '{print $1,$2,$4,$3}' \
    | diff - <(awk -v n="$n" '!/^#/ && $1==n' "$REPO/results_exact.txt") && echo "kmec_cluster n=$n matches results_exact.txt"
done
echo "ulimit -n (open files) = $(ulimit -n); the partition step opens P bucket files per task"
echo "Setup done."
