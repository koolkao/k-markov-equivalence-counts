#!/bin/bash
# Phase 1: write every unlabeled DAG on N nodes to its bucket file. One array task per part m of M.
#   sbatch --array=0-63 --export=ALL,N=9,P=512,M=64,S=64 cluster/partition.sbatch
# Output: $SCRATCH/kmec/n$N_P$P/parts/part<m>/b<b>.bin and part<m>/DONE
#SBATCH --job-name=kmec-part
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=logs/part_%A_%a.log
set -euo pipefail
export PATH=$HOME/kmec-bin:$PATH
ulimit -n 8192 2>/dev/null || ulimit -n "$(ulimit -Hn)"
BASE=$SCRATCH/kmec/n${N}_P${P}
mkdir -p "$BASE/parts"
m=$SLURM_ARRAY_TASK_ID
echo "host $(hostname) start $(date -Is)"
/usr/bin/time -v kmec_cluster partition "$N" "$P" "$M" "$m" "$S" "$SLURM_CPUS_PER_TASK" "$BASE/parts"
echo "end $(date -Is)"
