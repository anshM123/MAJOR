#!/bin/bash
# per-dimension convergence scan of the direct D' integral (E2, J = 0) at Delta = $1, base grid 10,10,6,10,10,16,10
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
D=$1
for g in 14,10,6,10,10,16,10 10,14,6,10,10,16,10 10,10,9,10,10,16,10 10,10,6,14,10,16,10 10,10,6,10,14,16,10 10,10,6,10,10,24,10 10,10,6,10,10,16,14; do
  f=out/dq_D${D}_E2_${g//,/-}.npz
  [ -f $f ] || python3 md_dquad.py $D E2 $g > out/log_dq_${D}_E2_${g}.txt 2>&1
done
