#!/bin/bash
# fw queue A2 (FOLD_WINDOW_LARGEJ.md; supersedes fw_queueA.sh after its first node): direct fold-window nodes, reordered by
# decisiveness: (400, 100) [Delta_g 300], (1064, 64), (1100, 100) [Delta_g 1000], (364, 64), then (450, 150), (1200, 200) (2 nodes each:
# precision quantification).  n3 contour, nu_f = 13 below 13; FW_M = 10, FW_R = 8.  Resumable (finished nodes skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "fw_points.py 100 16.7 400" > /dev/null; do sleep 20; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 100 $nu 400; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 64 $nu 1064; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 100 $nu 1100; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 64 $nu 364; done
for nu in 16.7 25.0; do python3 fw_points.py 150 $nu 450; done
for nu in 16.7 25.0; do python3 fw_points.py 200 $nu 1200; done
echo FW_QUEUE_A2_DONE
