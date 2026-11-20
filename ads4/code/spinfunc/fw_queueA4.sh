#!/bin/bash
# fw queue A4 (supersedes fw_queueA3.sh after its J = 100 nodes): core 1 -- (1064, 64) fold + bulk, (364, 64) fold, (450, 150) probe.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "^python3 fw_points.py" > /dev/null; do sleep 20; done
for nu in 16.7 25.0 13.3 33.4 20.8 10.0; do python3 fw_points.py 64 $nu 1064; done
for x in 0.2 0.4 0.6; do python3 fw_points.py 64 $(python3 -c "print(round($x*(1064-66)**0.5*(1064+63)**0.5, 6))") 1064; done
for nu in 16.7 25.0 13.3 33.4 20.8; do python3 fw_points.py 64 $nu 364; done
for nu in 16.7 25.0; do python3 fw_points.py 150 $nu 450; done
echo FW_QUEUE_A4_DONE
