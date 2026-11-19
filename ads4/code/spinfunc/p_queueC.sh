#!/bin/bash
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1
while pgrep -f "p3_run.py C3H 5" > /dev/null; do sleep 15; done
L1_DEG=8 P3_DS=100 python3 p3_run.py C2 7 >> cache/p3_run_k3_C2.log 2>&1
L1_DEG=8 P3_DS=100 python3 p3_run.py C3H 7 >> cache/p3_run_k3_C3H.log 2>&1
