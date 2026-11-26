#!/bin/bash
# Phase 2: count classes in each bucket. One array task per bucket b of P.
#   sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=2,KMAX=5 cluster/reduce.sbatch
#   sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=6,KMAX=6 cluster/reduce.sbatch
#   sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=7,KMAX=7 cluster/reduce.sbatch
# If a bucket fails (e.g. out of memory), rerun it as MSUB sub-buckets, one task per J:
#   sbatch --array=0-3 --mem=64G --export=ALL,N=9,P=512,KMIN=7,KMAX=7,B=137,MSUB=4 cluster/reduce.sbatch
# (with B set, the array index is the sub-bucket J). combine.sh accepts either form per bucket.
#SBATCH --job-name=kmec-red
#SBATCH --partition=normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=logs/red_%A_%a.log
set -euo pipefail
export PATH=$HOME/kmec-bin:$PATH
BASE=$SCRATCH/kmec/n${N}_P${P}
RES=$BASE/res_k${KMIN}-${KMAX}
mkdir -p "$RES"
if [ -n "${B:-}" ]; then b=$B; msub=$MSUB; j=$SLURM_ARRAY_TASK_ID; else b=$SLURM_ARRAY_TASK_ID; msub=1; j=0; fi
out=$RES/b${b}_s${j}of${msub}
echo "host $(hostname) bucket $b sub $j/$msub start $(date -Is)"
/usr/bin/time -v kmec_cluster reduce "$N" "$KMIN" "$KMAX" "$P" "$b" "$msub" "$j" "$SLURM_CPUS_PER_TASK" "$BASE/parts" > "$out.tmp" 2> "$out.err"
mv "$out.tmp" "$out.txt"
echo "end $(date -Is)"
grep -E "Maximum resident|Elapsed" "$out.err" || true
