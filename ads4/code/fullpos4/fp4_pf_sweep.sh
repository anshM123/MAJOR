#!/bin/bash
# FULLPOS4 §1: verification sweep of the outgoing contour (fp4_pf.py) against the exact terminating kernel
ulimit -v 4000000; export OMP_NUM_THREADS=1
for st in "3 0.8 2 0 7.5" "3 1.2 10 0 32.5" "4 2.7 7 0 25.25" "3 0.6 9 4 30.7" "4 4.1 10 4 40.3" "3 5.5 6 0 19.9" "4 0.3 8 4 23.75" "3 1.7 3 0 9.3" "3 12.0 10 0 45.5"; do
  python3 fp4_pf.py $st def + 2>&1 | grep -v Warn
done
