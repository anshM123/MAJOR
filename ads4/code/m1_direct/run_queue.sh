#!/bin/bash
# queue of direct D' runs: args = "Delta kind grid" triples separated by ';'
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
IFS=';' read -ra JOBS <<< "$1"
for j in "${JOBS[@]}"; do
  set -- $j
  f=out/dq_D$1_$2_${3//,/-}.npz
  [ -f $f ] || python3 md_dquad.py $1 $2 $3 > out/log_dq_$1_$2_$3.txt 2>&1
done
