#!/bin/bash
# per-dimension convergence scan of the D' rule for the EE integrand
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
D=${1:-30}
for g in 12,10,6,12,10,20,10 12,10,6,12,10,10,16 12,14,6,12,10,10,10 12,10,6,12,14,10,10 12,10,6,12,10,28,10; do
  [ -f cache/dq_D${D}_E_${g//,/-}.npz ] || python3 pcyc_dquad.py $D E $g > cache/log_dq_${D}_E_$g.txt 2>&1
done
