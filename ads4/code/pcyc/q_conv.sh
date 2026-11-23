#!/bin/bash
# convergence scans of the D' product rule (pcyc_dquad) -- 2 chains
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
( for g in 12,10,6,12,10,10,10 16,14,8,16,14,14,14 20,18,10,20,18,16,16; do python3 pcyc_dquad.py 30 S $g > cache/log_dq_30_S_$g.txt 2>&1; done ) &
( for g in 12,10,6,12,10,10,10 16,14,8,16,14,14,14; do python3 pcyc_dquad.py 30 E $g > cache/log_dq_30_E_$g.txt 2>&1; done ) &
wait
