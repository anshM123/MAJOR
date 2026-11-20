#!/bin/bash
# fw queue A (FOLD_WINDOW_LARGEJ.md): direct fold-window nodes, Delta_g = 300, J = 100 then 150 (n3 contour, nu_f = 13 below 13).
# Resumable: finished nodes are skipped (re-run this script after a restart).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 100 $nu 400; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 150 $nu 450; done
echo FW_QUEUE_A_DONE
