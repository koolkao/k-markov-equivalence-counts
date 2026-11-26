#!/usr/bin/env bash
# List the largest buckets after the partition, for piloting the reduce on the hardest cases.
#   bash cluster/rank_buckets.sh N P [count]
# Prints "bucket DAGs" for the largest buckets, then a comma-separated list for sbatch --array=.
set -euo pipefail
N=$1; P=$2; C=${3:-8}
BASE=${BASE:-$SCRATCH/kmec/n${N}_P${P}}
# Sizes come from file metadata (no reading of the ~62 GB of data). awk reads all input, so no SIGPIPE.
find "$BASE/parts" -name 'b*.bin' -printf '%f %s\n' \
  | awk '{b=$1; sub(/^b/,"",b); sub(/\.bin$/,"",b); s[b]+=$2} END{for (b in s) print b, s[b]/18}' \
  | sort -k2nr \
  | awk -v c="$C" 'NR<=c {print > "/dev/stderr"; l = l (NR>1?",":"") $1} END{print l}'
