#!/bin/bash
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1
L1_DEG=8 P3_DS=100 python3 p3_run.py C2 5 >> cache/p3_run_k3_C2.log 2>&1
L1_DEG=8 P3_DS=100 python3 p3_run.py C3H 5 >> cache/p3_run_k3_C3H.log 2>&1
L1_DEG=6 P3_ANTI=1 python3 -c "
import p3_odd_action as P3
for D in (400.37,200.37,800.37,100.37):
    m=P3.mphys(D,19); P3.get('C2',19,D,[x*m for x in (0.2,0.4)])
" >> cache/p5_J19.log 2>&1
