# Running n = 9 on Stanford's Sherlock cluster

This runs the two-phase pipeline in `kmec_cluster.c`:
1. **Partition:** every unlabeled DAG goes to a bucket file on `$SCRATCH`. The bucket is decided by a
   class invariant, so each class lands in exactly one bucket.
2. **Reduce:** each bucket counts its own classes exactly.
3. **Combine:** sums the buckets and checks the totals.

Expected cost for n = 9: about 20 core-hours for the partition and about 300 core-hours for the
reduce. Sherlock's shared partitions are free.

## 0. One-time setup (login node)

```
git clone https://github.com/koolkao/k-markov-equivalence-counts.git
cd k-markov-equivalence-counts
bash cluster/setup_sherlock.sh        # builds nauty + kmec_cluster into ~/kmec-bin, runs checks
```

The script should end with "kmec_cluster n=7 matches results_exact.txt" and "Setup done.". If gcc
is too old for `-march=x86-64-v2`, run `module load gcc` first and rerun the script.

## 1. Validate on n = 8 (about 15 minutes in total)

```
sbatch --array=0-3   --export=ALL,N=8,P=16,M=4,S=32 cluster/partition.sbatch
# when those finish:
sbatch --array=0-15  --export=ALL,N=8,P=16,KMIN=2,KMAX=6 cluster/reduce.sbatch
# when those finish:
bash cluster/combine.sh 8 16 2 6      # must print ALL CHECKS PASSED and match results_exact.txt
```

## 2. n = 9

```
sbatch --array=0-63  --export=ALL,N=9,P=512,M=64,S=64 cluster/partition.sbatch     # about 62 GB on $SCRATCH
# when done:
sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=2,KMAX=5 cluster/reduce.sbatch
sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=6,KMAX=6 cluster/reduce.sbatch
sbatch --array=0-511 --export=ALL,N=9,P=512,KMIN=7,KMAX=7 --mem=48G cluster/reduce.sbatch
```

To pilot first, submit `--array=0-7` of the k = 7 reduce, check memory with `sacct` or the
`Maximum resident` line in the logs, then submit the rest.

When a reduce task runs out of memory, rerun that bucket as sub-buckets with more memory. For
example, bucket 137 of the k = 7 run:

```
sbatch --array=0-3 --mem=64G --export=ALL,N=9,P=512,KMIN=7,KMAX=7,B=137,MSUB=4 cluster/reduce.sbatch
```

Finally:

```
bash cluster/combine.sh 9 512 2 5
bash cluster/combine.sh 9 512 6 6
bash cluster/combine.sh 9 512 7 7     # checks T(9,7) = A007984(9) = 326266056291213
```

k = 0..2 for n = 9 runs on a single machine with `kmec_cluster direct 9 0 2 <threads> 512`. The
k = 2 result from the cluster reduce is an independent cross-check of it.

## Notes

- `$SCRATCH` files not modified for 90 days are purged. Copy `res_k*/combined.txt` and the logs
  back to the repository when done.
- Each partition task opens P bucket files. The script raises `ulimit -n`.
