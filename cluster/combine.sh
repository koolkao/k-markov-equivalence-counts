#!/usr/bin/env bash
# Check the partition and sum the reduce results for one k range.
#   bash cluster/combine.sh N P KMIN KMAX
# Reads $SCRATCH/kmec/n<N>_P<P>/ (override with BASE=...). Exits non-zero on any problem:
#  * missing partition parts, parts that disagree on n, P or M, or a generated total != A003087(N);
#  * a bucket without a complete result: either b<b>_s0of1.txt, or all of b<b>_s<j>of<m>.txt for
#    some m (if several sub-bucket splits exist, the first complete one is used);
#  * a result whose header does not match (N, KMIN, KMAX, P, b, msub, j) or that lacks a k row;
#  * ERROR lines in the .err file of any result that is used;
#  * kept total != A003087(N) (every DAG must be counted in exactly one bucket);
#  * a known value that does not match (T(N,0) = A046165, T(N,N-2) = A007984, U(N,0) = A048194).
set -euo pipefail
shopt -s nullglob
N=$1; P=$2; KMIN=$3; KMAX=$4
BASE=${BASE:-$SCRATCH/kmec/n${N}_P${P}}
RES=$BASE/res_k${KMIN}-${KMAX}
A003087=(1 1 2 6 31 302 5984 243668 20286025 3424938010)
fail() { echo "FAIL: $*" >&2; exit 1; }

# --- partition
dones=("$BASE"/parts/part*/DONE)
[ ${#dones[@]} -gt 0 ] || fail "no partition DONE files in $BASE/parts"
M=$(awk 'NR==1{print $6}' "${dones[0]}")
[ "${#dones[@]}" = "$M" ] || fail "partition: ${#dones[@]} of $M parts done"
awk -v n="$N" -v p="$P" -v m="$M" '$2!=n || $4!=p || $6!=m {bad=1} END{exit bad}' "${dones[@]}" || fail "DONE files disagree on n, P or M"
gen=$(awk '{s+=$NF} END{printf "%.0f", s}' "${dones[@]}")
[ "$gen" = "${A003087[$N]}" ] || fail "partition generated $gen DAGs, expected ${A003087[$N]}"
echo "partition OK: $M parts, $gen DAGs"

# --- one complete result set per bucket
files=()
for b in $(seq 0 $((P-1))); do
  if [ -s "$RES/b${b}_s0of1.txt" ]; then files+=("$RES/b${b}_s0of1.txt"); continue; fi
  chosen=""
  for m in $(for f in "$RES"/b${b}_s*of*.txt; do f=${f##*of}; echo "${f%.txt}"; done | sort -un); do
    [ "$m" -gt 1 ] || continue
    ok=1; for j in $(seq 0 $((m-1))); do [ -s "$RES/b${b}_s${j}of${m}.txt" ] || ok=0; done
    if [ $ok = 1 ]; then chosen=$m; break; fi
  done
  [ -n "$chosen" ] || fail "bucket $b has no complete result"
  for j in $(seq 0 $((chosen-1))); do files+=("$RES/b${b}_s${j}of${chosen}.txt"); done
done

# --- validate the chosen files
nk=$((KMAX-KMIN+1))
for f in "${files[@]}"; do
  base=$(basename "$f" .txt); b=${base#b}; b=${b%%_s*}; j=${base#*_s}; j=${j%%of*}; m=${base##*of}
  hdr=$(sed -n 1p "$f")
  [[ "$hdr" == "# n $N kmin $KMIN kmax $KMAX P $P b $b msub $m j $j "* ]] || fail "bad header in $f: $hdr"
  rows=$(awk '!/^#/' "$f" | wc -l | tr -d ' ')
  [ "$rows" = "$nk" ] || fail "$f has $rows k rows, expected $nk"
  if [ -f "${f%.txt}.err" ] && grep -q ERROR "${f%.txt}.err"; then fail "ERROR in ${f%.txt}.err"; fi
done
kept=$(for f in "${files[@]}"; do sed -n 1p "$f"; done | awk '{s+=$NF} END{printf "%.0f", s}')
[ "$kept" = "${A003087[$N]}" ] || fail "kept total $kept != A003087($N) = ${A003087[$N]}"
echo "reduce OK: ${#files[@]} result files, kept total $kept"

# --- sum and check known values
for f in "${files[@]}"; do awk '!/^#/' "$f"; done \
  | awk '{U[$2]+=$3; T[$2]+=$4; n=$1} END{for (k in U) printf "%d %d %.0f %.0f\n", n, k, T[k], U[k]}' | sort -k2n > "$RES/combined.txt"
echo "# n k labeled unlabeled"; cat "$RES/combined.txt"
awk 'BEGIN{ A046165[9]=152424420; A007984[9]=326266056291213; A048194[9]=2223;
            A046165[8]=3731508;   A007984[8]=212133402500;    A048194[8]=557;
            A046165[7]=129425;    A007984[7]=312510571;       A048194[7]=164; bad=0 }
  { n=$1; k=$2;
    if (k==0   && (n in A046165)) { ok=($3==A046165[n]); print (ok?"OK  ":"FAIL"), "T("n",0) = A046165"; bad+=!ok }
    if (k==n-2 && (n in A007984)) { ok=($3==A007984[n]); print (ok?"OK  ":"FAIL"), "T("n","k") = A007984"; bad+=!ok }
    if (k==0   && (n in A048194)) { ok=($4==A048194[n]); print (ok?"OK  ":"FAIL"), "U("n",0) = A048194"; bad+=!ok } }
  END { exit bad>0 }' "$RES/combined.txt" || fail "a known-value check failed"
echo "ALL CHECKS PASSED"
