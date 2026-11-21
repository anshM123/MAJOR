#!/bin/sh
# Stage 3 background queue: exact forward moments, 2 workers (one core each).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for w in 0 1; do
  ( for Dg in 100 300 1000; do python3 stage3_fwexact.py $Dg $w 2; done ) > cache/stage3_fw_w$w.log 2>&1 &
done
wait
