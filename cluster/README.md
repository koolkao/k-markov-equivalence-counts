# Running n = 9 on Stanford's Sherlock cluster

This runs the two-phase pipeline in `kmec_cluster.c`:
1. **Partition:** every unlabeled DAG goes to a bucket file on `$SCRATCH`. The bucket is decided by a
   class invariant, so each class lands in exactly one bucket.
2. **Reduce:** each bucket counts its own classes exactly.
3. **Combine:** sums the buckets and checks everything.

Expected cost for n = 9: about 20 core-hours for the partition and a few hundred core-hours for the
reduce. Sherlock's shared partitions are free.

Run every command from the repository root, because Slurm writes the logs to `logs/`.

## 0. One-time setup (login node)

```
git clone https://github.com/koolkao/k-markov-equivalence-counts.git
cd k-markov-equivalence-counts
bash cluster/setup_sherlock.sh
```

This builds nauty and kmec_cluster into `~/kmec-bin` and runs checks. It must end with
"Setup done."; it exits with an error if any check fails. If gcc is too old for
`-march=x86-64-v2`, run `module load gcc` first and rerun the script.

## 1. Validate on n = 8 (about 15 minutes in total)

```
sbatch --array=0-3  --export=ALL,N=8,P=64,M=4,S=64 cluster/partition.sbatch
# after they finish:
sbatch --array=0-63 --export=ALL,N=8,P=64,KMIN=6,KMAX=6 cluster/reduce.sbatch
# after they finish:
bash cluster/combine.sh 8 64 6 6   # must print ALL CHECKS PASSED, and k=6 must equal results_exact.txt
```

## 2. n = 9: partition

```
sbatch --array=0-63 --export=ALL,N=9,P=512,M=64,S=64 cluster/partition.sbatch     # about 62 GB on $SCRATCH
```

The reduce step refuses to run until all 64 parts are done. To submit the reduce in advance, add
`--dependency=afterok:<partition job id>`.

## 3. n = 9: reduce, one job per k range

Memory grows with k, so each large k gets its own job:

| k range | `--mem` to start with |
|---|---|
| `KMIN=2,KMAX=3` | 16G |
| `KMIN=4,KMAX=4` | 32G |
| `KMIN=5,KMAX=5` | 48G |
| `KMIN=6,KMAX=6` | 48G |
| `KMIN=7,KMAX=7` | 48G |

**Pilot first, on the largest buckets.** The typical bucket is far from the worst case:

```
big=$(bash cluster/rank_buckets.sh 9 512 8)       # the 8 largest buckets, e.g. "229,292,72,..."
sbatch --array=$big --mem=48G --export=ALL,N=9,P=512,KMIN=7,KMAX=7 cluster/reduce.sbatch
# check memory in the logs ("Maximum resident set size") or with: sacct -j <jobid> --format=JobID,MaxRSS,Elapsed
```

Then submit all buckets for each k range, adjusting `--mem` from the pilot. Sherlock limits how many
jobs one user can have queued (see `sacctmgr show qos normal format=MaxSubmitPU`); if `sbatch`
refuses, submit one or two k ranges at a time and the next ones when they finish:

```
sbatch --array=0-511 --mem=16G --export=ALL,N=9,P=512,KMIN=2,KMAX=3 cluster/reduce.sbatch
sbatch --array=0-511 --mem=32G --export=ALL,N=9,P=512,KMIN=4,KMAX=4 cluster/reduce.sbatch
sbatch --array=0-511 --mem=48G --export=ALL,N=9,P=512,KMIN=5,KMAX=5 cluster/reduce.sbatch
sbatch --array=0-511 --mem=48G --export=ALL,N=9,P=512,KMIN=6,KMAX=6 cluster/reduce.sbatch
sbatch --array=0-511 --mem=48G --export=ALL,N=9,P=512,KMIN=7,KMAX=7 cluster/reduce.sbatch
```

**When a bucket fails** (for example, out of memory), rerun just that bucket as sub-buckets.
Sub-buckets split by the KMIN-signature, so even the largest bucket divides evenly. For example,
bucket 229 of the k = 7 job, in 8 parts:

```
sbatch --array=0-7 --mem=48G --export=ALL,N=9,P=512,KMIN=7,KMAX=7,B=229,MSUB=8 cluster/reduce.sbatch
```

## 4. Combine and check

```
bash cluster/combine.sh 9 512 2 3
bash cluster/combine.sh 9 512 4 4
bash cluster/combine.sh 9 512 5 5
bash cluster/combine.sh 9 512 6 6
bash cluster/combine.sh 9 512 7 7
```

Each run checks:
- that the partition is complete and its total is A003087(9);
- that every bucket has a complete result;
- that the result headers match and no ERROR lines occur;
- that the kept total is A003087(9).

The k = 7 run also checks T(9,7) = A007984(9) = 326266056291213.

k = 0..2 is computed on a single machine (`kmec_cluster direct 9 0 2 <threads> 512`). The cluster's
k = 2 result is an independent cross-check of it.

## Notes

- `$SCRATCH` files not modified for 90 days are purged. Copy `res_k*/combined.txt` and the logs
  back to the repository when done.
- Each partition task opens P bucket files; the script raises `ulimit -n`.
