#!/bin/bash
# production D' quadratures (pcyc_dquad, nu = 1..16), one chain per Delta
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
D=$1
for job in "S 12,12,6,12,12,12,12" "S 14,16,8,14,16,16,16" "E 12,12,6,12,12,20,12" "E 12,14,8,12,14,24,16"; do
  set -- $job
  [ -f cache/dq_D${D}_$1_${2//,/-}.npz ] || python3 pcyc_dquad.py $D $1 $2 > cache/log_dq_${D}_$1_$2.txt 2>&1
done
