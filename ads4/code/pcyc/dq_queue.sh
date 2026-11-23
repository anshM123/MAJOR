#!/bin/bash
# production D' quadratures (pcyc_dquad), one chain
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
for job in "30 E 12,10,6,12,10,10,10" "30 E 12,14,8,12,14,14,14" "30 S 14,16,8,14,16,16,16" \
           "100 S 12,10,6,12,10,10,10" "100 E 12,10,6,12,10,10,10" "100 S 14,16,8,14,16,16,16" "100 E 12,14,8,12,14,14,14"; do
  set -- $job
  [ -f cache/dq_D$1_$2_${3//,/-}.npz ] || python3 pcyc_dquad.py $1 $2 $3 > cache/log_dq_$1_$2_$3.txt 2>&1
done
