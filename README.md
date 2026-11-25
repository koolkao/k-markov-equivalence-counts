# Counts of k-Markov equivalence classes of DAGs

T(n,k) is the number of equivalence classes of directed acyclic graphs (DAGs) on n nodes, where two
DAGs are equivalent if they entail the same d-separation statements u ⊥ v | S for every conditioning
set S ⊆ V \ {u, v} with |S| <= k (0 <= k <= n-2).

- **Names for this relation.** Kocaoglu (NeurIPS 2023, Definition 3.1) calls it *k-Markov
  equivalence*; Kook and Mogensen (2026) call it *k-weak equivalence*. It is also the special case
  C = {S : |S| <= k} of the *C-Markov equivalence* of Lee, Ribeiro and Kocaoglu (UAI 2025).
- **Why pairwise statements are enough.** Set statements X ⊥ Y | S hold iff they hold for every pair.
- **The ends of the range.** k = n-2 is ordinary Markov equivalence. k = 0 compares only marginal
  independences.

The counts are computed by brute force from the definition. Every DAG is enumerated, its
d-separation statements are computed, and distinct statement sets are counted. No characterization
theorem is used.

## Results

### Labeled DAGs (`results_labeled.txt`; n = 7 in `results_n7_seed*.txt`)

| n \ k | 0 | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|---|
| 2 | 2 | | | | | |
| 3 | 8 | 11 | | | | |
| 4 | 49 | 137 | 185 | | | |
| 5 | 462 | 2962 | 7452 | 8782 | | |
| 6 | 6424 | 99640 | 554950 | 982070 | 1067825 | |
| 7 | 129425 | 4882180 | 66029917 | 213365427 | 299380027 | 312510571 |

### Up to isomorphism (`results_unlabeled.txt`; columns are n, k, labeled, unlabeled)

| n \ k | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| 2 | 2 | | | | |
| 3 | 4 | 5 | | | |
| 4 | 9 | 16 | 20 | | |
| 5 | 21 | 64 | 122 | 142 | |
| 6 | 56 | 329 | 1222 | 2015 | 2201 |

The right edge of this table (2, 5, 20, 142, 2201) is the number of ordinary Markov equivalence
classes (essential graphs) up to isomorphism.

## Programs

| file | what it does | range |
|---|---|---|
| `kmec.c` | Exact signatures; d-separation on the moralized ancestral graph; single thread | n <= 6 |
| `kmec2.c` | Node-by-node DAG enumeration; bitmask Bayes-ball d-separation; 128-bit streamed fingerprints; sharded hash set; pthreads. Usage: `./kmec2 n [kmax] [threads] [seed]` (the built-in usage message omits `seed`) | n <= 7 |
| `kmec_unl.c` | Unlabeled counts: exact labeled classes, each canonicalized over all n! relabelings; also counts unlabeled DAGs as a check | n <= 6 |
| `kmec_ref.py` | Independent reference using networkx `is_d_separator` (needs networkx >= 3.3) | n <= 5 |
| `tests/test_bayesball.c` | Checks kmec2's Bayes-ball routine against moralization d-separation on random DAGs, n = 3..8 | |
| `tests/recount_ref.py` | Independent Python recount with its own d-separation and canonicalization, labeled and unlabeled | n <= 5 |

## Reproducing

```
make            # build everything (cc -O3 -march=native -pthread)
make check      # Bayes-ball test, kmec2 vs results_labeled.txt, networkx reference for n = 5
make labeled    # regenerate results_labeled.txt (n = 2..6), about 1 minute
make unlabeled  # regenerate results_unlabeled.txt (n = 2..6), about 3 minutes, about 2 GB RAM
make n7         # n = 7 with seeds 0 and 7, then diff them: 35-50 minutes and about 30 GB RAM per run
python3 tests/recount_ref.py 5
```

## Validation

- **Agreement between implementations.** `kmec.c` and `kmec_ref.py` agree for n = 2..5. `kmec2.c`
  agrees exactly with `kmec.c` for n = 2..6, with 1 to 256 threads and several seeds.
  `tests/recount_ref.py` reproduces both tables for n <= 5.
- **Bayes-ball test.** `tests/test_bayesball.c` runs 22,512,000 random (DAG, u, v, S) cases with zero
  mismatches. The routine follows Shachter (1998), Bayes-Ball, Algorithm 2 and Theorem 2.
- **Known sequences reproduced.**
  - Labeled DAG counts (3, 25, 543, 29281, 3781503, 1138779265) match OEIS A003024.
  - Unlabeled DAG counts (2, 6, 31, 302, 5984) match A003087.
  - The labeled right edge k = n-2 matches A007984 (2, 11, 185, 8782, 1067825, 312510571).
  - The labeled column k = 0 matches A046165 (2, 8, 49, 462, 6424, 129425), as the bijection with
    minimal covers in Deligeorgaki, Markham, Misra and Solus (2023) predicts.
  - The unlabeled column k = 0 matches A048194 (2, 4, 9, 21, 56). This is consistent with Textor,
    Idelberger and Liskiewicz (2015).
- **Unlabeled right edge.** For n <= 4 it agrees with the list of unlabeled essential graphs in
  Andersson, Madigan and Perlman (1997), Table 1. An independent reviewer's count based on skeletons
  and v-structures, using no d-separation code, gives 2, 5, 20, 142, 2201.
- **Row n = 7 uses hashed signatures.** The run has to use 128-bit fingerprints of the signatures, so
  a collision could merge two classes, but nothing can split one. Three runs with different hash
  initializations give identical counts:
  - seed A: an earlier build whose hash function equals the current code's seed 0; that source is not
    archived;
  - seed 7;
  - seed 0: rerun with the archived source; `n7_seed0_source.sha256` holds the SHA-256 of `kmec2.c`
    used for the run.

  The two ends of row 7 and the DAG count match OEIS. The run logs are `n7_seed*.log`; the
  environment is in `n7_seed0_env.txt` (Apple M4 Pro, 14 threads, Apple clang 21, macOS 26). The runs
  took 19-47 minutes each (19 minutes on an otherwise idle machine), with peak memory of 29-31 GB.

## Novelty and literature

`oeis_search_2026-11-25.txt` records OEIS searches made on 2026-11-25 (search terms URL-encoded,
with a control query that finds A007984). None of these were found:
- the labeled or unlabeled triangle read by rows;
- the interior columns;
- the k = n-3 diagonal;
- any interior n = 7 term;
- the unlabeled Markov equivalence class counts 2, 5, 20, 142, 2201.

A full-text review of the literature found no published counts of these classes for 1 <= k <= n-3.
The k = 0 columns are known (see above). For unlabeled essential graphs, only n <= 4 is implicit in
Andersson, Madigan and Perlman (1997).

## References

- S. A. Andersson, D. Madigan, M. D. Perlman. A characterization of Markov equivalence classes for
  acyclic digraphs. Ann. Statist. 25 (1997) 505-541. doi:10.1214/aos/1031833662
- D. Deligeorgaki, A. Markham, P. Misra, L. Solus. Combinatorial and algebraic perspectives on the
  marginal independence structure of Bayesian networks. Algebraic Statistics 14 (2023) 233-286.
  doi:10.2140/astat.2023.14.233
- M. Kocaoglu. Characterization and learning of causal graphs with small conditioning sets.
  NeurIPS 36 (2023). arXiv:2301.09028
- L. Kook, S. W. Mogensen. Exact graph learning via integer programming. arXiv:2601.20589 (2026)
- K. Lee, B. Ribeiro, M. Kocaoglu. Constraint-based causal discovery from a collection of
  conditioning sets. UAI 2025, PMLR 286, 2486-2516
- R. D. Shachter. Bayes-Ball: The rational pastime. UAI 1998. arXiv:1301.7412
- J. Textor, A. Idelberger, M. Liskiewicz. Learning from pairwise marginal independencies.
  UAI 2015. arXiv:1508.00280
- M. Wienoebst, M. Liskiewicz. Recovering causal structures from low-order conditional
  independencies. AAAI 2020, 10302-10309. arXiv:2010.02675

## License and citation

MIT (see `LICENSE`). To cite, see `CITATION.cff`.
