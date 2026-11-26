#!/bin/bash
# DV runs of the Regge-contour decomposition (ct_check.py) + exact real-axis values (realint.py) at adversarial states
ulimit -v 4000000; export OMP_NUM_THREADS=1
cd /home/user/MAJOR/ads4/code/analytic
for st in "$@"; do
  set -- $st
  python3 ct_check.py $1 $2 $3 $4 ${5:-150} ${6:-3000} > logs/ct_$1_$2_$3_Y$4.log 2>&1
  python3 realint.py $1 $2 $3 ${7:-400} ${6:-3000} >> logs/ct_$1_$2_$3_Y$4.log 2>&1
done
