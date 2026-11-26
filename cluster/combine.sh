#!/usr/bin/env bash
# Check the partition and sum the reduce results.
#   ./cluster/combine.sh N P KMIN KMAX
# Reads $SCRATCH/kmec/n<N>_P<P>/ (override with BASE=...). Exits non-zero on any problem:
#  * missing or incomplete partition parts, or generated total != A003087(N);
#  * a bucket with no complete result (either b<b>_s0of1.txt, or all of b<b>_s<j>of<m>.txt, j < m);
#  * a result whose header does not match (N, KMIN, KMAX, P, b, msub, j) or that lacks a k row;
#  * ERROR lines in any .err file;
#  * kept total != A003087(N) (every DAG must be counted in exactly one bucket);
#  * a known value that does not match (T(N,0) = A046165, T(N,N-2) = A007984, U(N,0) = A048194).
set -euo pipefail
N=$1; P=$2; KMIN=$3; KMAX=$4
BASE=${BASE:-$SCRATCH/kmec/n${N}_P${P}}
RES=$BASE/res_k${KMIN}-${KMAX}
A003087=(1 1 2 6 31 302 5984 243668 20286025 3424938010)
fail() { echo "FAIL: $*" >&2; exit 1; }

# partition
nd=$(ls "$BASE"/parts/part*/DONE 2>/dev/null | wc -l | tr -d " ")
M=$(awk '{print $6}' "$(ls "$BASE"/parts/part*/DONE | head -1)")
[ "$nd" = "$M" ] || fail "partition: $nd of $M parts done"
gen=$(cat "$BASE"/parts/part*/DONE | awk '{s+=$NF} END{printf "%.0f", s}')
[ "$gen" = "${A003087[$N]}" ] || fail "partition generated $gen DAGs, expected ${A003087[$N]}"
echo "partition OK: $M parts, $gen DAGs"

# collect exactly one complete result set per bucket
files=()
for b in $(seq 0 $((P-1))); do
  if [ -s "$RES/b${b}_s0of1.txt" ]; then files+=("$RES/b${b}_s0of1.txt"); continue; fi
  subs=( $(ls "$RES"/b${b}_s*of*.txt 2>/dev/null | grep -v "_s0of1.txt" || true) )
  [ ${#subs[@]} -gt 0 ] || fail "bucket $b has no result"
  m=$(basename "${subs[0]}" .txt | sed 's/.*of//')
  for j in $(seq 0 $((m-1))); do [ -s "$RES/b${b}_s${j}of${m}.txt" ] || fail "bucket $b sub-bucket $j of $m missing"; files+=("$RES/b${b}_s${j}of${m}.txt"); done
done
if grep -l ERROR "$RES"/*.err >/dev/null 2>&1; then fail "ERROR lines in: $(grep -l ERROR "$RES"/*.err | tr '\n' ' ')"; fi

# validate headers and k rows, then sum
nk=$((KMAX-KMIN+1))
for f in "${files[@]}"; do
  base=$(basename "$f" .txt); b=${base#b}; b=${b%%_s*}; j=${base#*_s}; j=${j%%of*}; m=${base##*of}
  hdr=$(head -1 "$f")
  echo "$hdr" | grep -q "^# n $N kmin $KMIN kmax $KMAX P $P b $b msub $m j $j " || fail "bad header in $f: $hdr"
  [ "$(grep -vc '^#' "$f")" = "$nk" ] || fail "wrong number of k rows in $f"
done
kept=$(for f in "${files[@]}"; do head -1 "$f"; done | awk '{s+=$NF} END{printf "%.0f", s}')
[ "$kept" = "${A003087[$N]}" ] || fail "kept total $kept != A003087($N) = ${A003087[$N]}"
echo "reduce OK: ${#files[@]} result files, kept total $kept"
for f in "${files[@]}"; do grep -v '^#' "$f"; done \
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
