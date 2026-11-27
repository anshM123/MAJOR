#!/bin/bash
# grid convergence of the direct D' integral, kind E2 (EE k=2, J=0) and S (scalar), Delta = 10
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
D=${1:-10}
for g in 8,8,6,8,8,12,8 10,10,6,10,10,16,10 12,12,8,12,12,20,12; do
  python3 md_dquad.py $D E2 $g > out/log_dq_${D}_E2_${g}.txt 2>&1
done
