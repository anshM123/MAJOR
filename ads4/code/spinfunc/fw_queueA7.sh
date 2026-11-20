#!/bin/bash
# fw queue A7 (core 1): larger-x nodes at (1100, 100), x = 0.2, 0.1 (nu = 218.8, 109.4), needed to fix the per-row normalisation a_h
# (with only first-lobe fold nodes, (2J+1)x <= 6, the a_h fit is ill-conditioned: §2.6).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
for x in 0.2 0.1; do python3 fw_points.py 100 $(python3 -c "print(round($x*(1100-102)**0.5*(1100+99)**0.5, 6))") 1100; done
echo FW_QUEUE_A7_DONE
