#!/usr/bin/env bash
# One-time setup on a Sherlock login node: build nauty (geng, directg) and kmec_cluster, then run
# quick checks. Usage, from the repository root:
#   bash cluster/setup_sherlock.sh
# Puts binaries in $HOME/kmec-bin and creates logs/ (Slurm does not create log directories).
# Compute nodes mix Intel and AMD CPUs, so everything is built for x86-64-v2, not -march=native.
# x86-64-v2 needs GCC >= 11. If the system gcc is older, run `module load gcc` first.
# Exits non-zero if any check fails.
set -euo pipefail
REPO=$(cd "$(dirname "$0")/.." && pwd)
BIN=$HOME/kmec-bin; SRC=$HOME/kmec-src
mkdir -p "$BIN" "$SRC" "$REPO/logs"
gcc --version | sed -n 1p
[ -x /usr/bin/time ] || { echo "/usr/bin/time not found (the sbatch scripts use it)"; exit 1; }
if [ ! -x "$BIN/directg" ] || [ ! -x "$BIN/geng" ]; then
  ( cd "$SRC"
    curl -fL -o nauty2_9_3.tar.gz https://pallini.di.uniroma1.it/nauty2_9_3.tar.gz
    tar xzf nauty2_9_3.tar.gz
    cd nauty2_9_3
    ./configure --enable-generic CFLAGS="-O3 -march=x86-64-v2" > configure.log 2>&1 || { echo "nauty configure failed; see $SRC/nauty2_9_3/configure.log"; exit 1; }
    make -j8 geng directg > make.log 2>&1 || { echo "nauty build failed; see $SRC/nauty2_9_3/make.log"; exit 1; }
    cp geng directg "$BIN"/ )
fi
gcc -O3 -march=x86-64-v2 -pthread -Wall -o "$BIN/kmec_cluster" "$REPO/kmec_cluster.c"
export PATH=$BIN:$PATH
fail() { echo "SETUP CHECK FAILED: $*" >&2; exit 1; }
# nauty reproduces A003087 (unlabeled DAG counts)
expect=( [5]=302 [6]=5984 [7]=243668 )
for n in 5 6 7; do
  got=$(geng -q "$n" | directg -a -u 2>&1 | sed -n 's/.* \([0-9]*\) digraphs generated.*/\1/p')
  [ "$got" = "${expect[$n]}" ] || fail "nauty n=$n gave '$got', expected ${expect[$n]}"
  echo "nauty n=$n: $got unlabeled DAGs (A003087) OK"
done
# kmec_cluster reproduces the exact counts
for n in 6 7; do
  out=$(kmec_cluster direct "$n" 0 $((n-2)) 4 64) || fail "kmec_cluster direct n=$n exited with an error"
  diff <(echo "$out" | awk '!/^#/{print $1,$2,$4,$3}') \
       <(awk -v n="$n" '!/^#/ && $1==n' "$REPO/results_exact.txt") > /dev/null || fail "kmec_cluster n=$n differs from results_exact.txt"
  echo "kmec_cluster n=$n matches results_exact.txt"
done
echo "ulimit -n (open files) = $(ulimit -n); partition tasks open P bucket files each"
echo "Setup done."
