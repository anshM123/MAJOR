#!/bin/bash
# fw queue B6 (core 2; supersedes fw_queueB5.sh after its first (1100, 100) node): remaining (1100, 100) fold nodes, then one
# (1200, 200) precision probe (the (1100, 100) bulk points are dropped: §4 shows 3 bulk points cannot test a bulk smooth term).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "^python3 fw_points.py 100 16.7 1100" > /dev/null; do sleep 20; done
for nu in 25.0 13.3 33.4 20.8; do python3 fw_points.py 100 $nu 1100; done
python3 fw_points.py 200 16.7 1200
echo FW_QUEUE_B6_DONE
