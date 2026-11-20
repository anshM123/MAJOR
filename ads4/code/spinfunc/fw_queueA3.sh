#!/bin/bash
# fw queue A3 (FOLD_WINDOW_LARGEJ.md; supersedes fw_queueA2.sh): direct nodes in priority order.  Fold window nu = 10 ... 33.4 (n3,
# nu_f = 13 below 13), bulk points x = 0.2 ... 0.6 (n3, nu_f = nu).  FW_M = 10, FW_R = 8.  Resumable (finished nodes skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "^python3 fw_points.py" > /dev/null; do sleep 20; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 100 $nu 400; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 64 $nu 1064; done
# bulk (1064, 64): x = 0.2, 0.4, 0.6  (m = 1063.97)
for x in 0.2 0.4 0.6; do python3 fw_points.py 64 $(python3 -c "print(round($x*(1064-66)**0.5*(1064+63)**0.5, 6))") 1064; done
for nu in 16.7 25.0 13.3 33.4 20.8; do python3 fw_points.py 100 $nu 1100; done
# bulk (1100, 100): x = 0.3, 0.5  (m = 1094.4)
for x in 0.3 0.5; do python3 fw_points.py 100 $(python3 -c "print(round($x*(1100-102)**0.5*(1100+99)**0.5, 6))") 1100; done
for nu in 16.7 25.0 13.3 33.4 20.8; do python3 fw_points.py 64 $nu 364; done
for nu in 16.7 25.0; do python3 fw_points.py 150 $nu 450; done
for nu in 16.7 25.0; do python3 fw_points.py 200 $nu 1200; done
echo FW_QUEUE_A3_DONE
