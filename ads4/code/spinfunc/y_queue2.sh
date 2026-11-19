#!/bin/bash
# EE_DG100.md: worker 2 after y1_queue.sh -- direct EE points on the Delta = 100 row, then even extra rows (K2)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f y1_queue.sh > /dev/null; do sleep 20; done
python3 y4_direct.py n3even 4 6 8 10
python3 y4_direct.py n1 4 6 8 10
python3 y4_direct.py n3odd 5 7 9
python3 y4_direct.py p8 5 7 9
python3 y5_n4.py run all >> cache/y5_w1.log 2>&1
python3 y1_rows.py evenx 4 1.05,1.1,1.2,1.3
python3 y1_rows.py evenx 10 1.1,1.2,1.3
