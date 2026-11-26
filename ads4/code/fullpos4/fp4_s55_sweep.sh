#!/bin/bash
# FULLPOS4 §1: author's §55 contour (fp4_s55.py) vs exact ladder sum; p = 3, 4; b = 0, 4; non-integer N; nu <= 1 and > 1
ulimit -v 4000000; export OMP_NUM_THREADS=1
for st in "2 7.5 0 3 1.2" "2 7.5 0 3 0.8" "2 5.5 0 3 0.8" "10 30.5 0 3 1.2" "7 23.25 0 4 2.7" "5 26.7 4 3 0.6" "6 34.3 4 4 4.1" "6 17.9 0 3 5.5" "4 19.75 4 4 0.3" "10 43.5 0 3 12.0" "8 40.2 4 3 7.7" "9 31.1 0 4 0.95"; do
  timeout 1800 python3 fp4_s55.py $st 32 40 2>&1 | grep -v Warn
done
