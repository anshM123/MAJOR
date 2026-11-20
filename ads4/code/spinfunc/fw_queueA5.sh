#!/bin/bash
# fw queue A5: core 2 (after the basis queue B4) -- (1100, 100) fold + bulk, (1200, 200) probe.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 FW_M=10 FW_R=8 Z2_NUSPLIT=9
while pgrep -f "fw_queueB4.sh" > /dev/null; do sleep 30; done
for nu in 16.7 25.0 13.3 33.4 20.8; do python3 fw_points.py 100 $nu 1100; done
for x in 0.3 0.5; do python3 fw_points.py 100 $(python3 -c "print(round($x*(1100-102)**0.5*(1100+99)**0.5, 6))") 1100; done
for nu in 16.7 25.0; do python3 fw_points.py 200 $nu 1200; done
echo FW_QUEUE_A5_DONE
