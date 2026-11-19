#!/bin/bash
# w1 queue A3 (replaces the rest of A2): near-gap-mass rows at Delta_g = 1000 (Delta = 1040, J^2/m ~ 1-2) and Delta_g = 300 (Delta = 336)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "w1_rows.py even EE,C3H 24 1600" > /dev/null; do sleep 20; done
W1_TAG=_a python3 w1_rows.py refb 36,40,44
W1_TAG=_a python3 w1_rows.py even EE,C3H 32 1040
W1_TAG=_a python3 w1_rows.py even EE,C3H 40 1040
W1_TAG=_a python3 w1_rows.py even EE,C3H 22 336
W1_TAG=_a python3 w1_rows.py even EE,C3H 36 1040
W1_TAG=_a python3 w1_rows.py even EE,C3H 44 1040
W1_TAG=_a python3 w1_rows.py even EE,C3H 24,28 800
W1_TAG=_a python3 w1_rows.py even EE,C3H 22,26,30 1600
