#!/bin/bash
# w1 queue A5 (replaces the rest of A3): Delta_g = 1000 near-gap rows (Delta = 1040) with 35 x nodes (x-resolution), J = 32, 40; (336, 22)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "w1_rows.py even EE,C3H 32 1040" > /dev/null; do sleep 20; done
MID=0.0625,0.0875,0.1125,0.1375,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775
W1_TAG=_a python3 w1_rows.py even EE,C3H 32 1040 $MID
W1_TAG=_a python3 w1_rows.py even EE,C3H 40 1040
W1_TAG=_a python3 w1_rows.py even EE,C3H 40 1040 $MID
W1_TAG=_a python3 w1_rows.py even EE,C3H 22 336
W1_TAG=_a python3 w1_rows.py even EE,C3H 48 1040
