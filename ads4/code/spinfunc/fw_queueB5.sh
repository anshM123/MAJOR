#!/bin/bash
# fw queue B5 (core 2, after B4; supersedes fw_queueA5.sh): fixed-engine K2 row (1064, 64) X35 (EE + C3H, state + reference) for the bulk
# smooth test at J = 64, then the direct nodes of (1100, 100) (fold + bulk) and the (1200, 200) probe.  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "fw_queueB4.sh" > /dev/null; do sleep 30; done
python3 fw_rows.py row EE,C3H 64 1064 X35
export FW_M=10 FW_R=8 Z2_NUSPLIT=9
for nu in 16.7 25.0 13.3 33.4 20.8; do python3 fw_points.py 100 $nu 1100; done
for x in 0.3 0.5; do python3 fw_points.py 100 $(python3 -c "print(round($x*(1100-102)**0.5*(1100+99)**0.5, 6))") 1100; done
for nu in 16.7 25.0; do python3 fw_points.py 200 $nu 1200; done
echo FW_QUEUE_B5_DONE
