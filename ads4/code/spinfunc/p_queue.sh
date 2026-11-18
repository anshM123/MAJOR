#!/bin/bash
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=6
if [ "$1" == A ]; then
  while pgrep -f p5_run.sh > /dev/null; do sleep 20; done
  python3 p3_run.py C3H 5 7 9 > cache/p3_run_C3H.log 2>&1
else
  while pgrep -f "p3_run.py C2" > /dev/null; do sleep 20; done
  P3_ANTI=1 python3 -c "
import p5_N as P, p3_odd_action as P3
for D in (200.37,400.37,800.37,100.37):
    m=P3.mphys(D,19); P3.get('C2',19,D,[x*m for x in (0.2,0.4)])
" > cache/p5_J19.log 2>&1
fi
