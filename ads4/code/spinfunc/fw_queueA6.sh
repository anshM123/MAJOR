#!/bin/bash
# fw queue A6 (core 1, after the (450, 150) probe node; supersedes the rest of fw_queueA4.sh): extra fold-window nodes at (364, 64) to reduce the fit-class ambiguity of the smooth part.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "^python3 fw_points.py 150" > /dev/null; do sleep 30; done
for nu in 18.5 11.7 29.0 15.0 23.0; do python3 fw_points.py 64 $nu 364; done
echo FW_QUEUE_A6_DONE
