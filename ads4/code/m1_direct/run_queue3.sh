#!/bin/bash
# queue of k = 3 (C3H) direct D' runs: "Delta J grid" triples separated by ';'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MD_KIND=C3
cd "$(dirname "$0")"
IFS=';' read -ra JOBS <<< "$1"
for j in "${JOBS[@]}"; do
  set -- $j
  f=out/dq3_D$1_J$2_${3//,/-}.npz
  [ -f $f ] || python3 md_dquadJ.py $1 $2 $3 > out/log_dq3_$1_$2_$3.txt 2>&1
done
