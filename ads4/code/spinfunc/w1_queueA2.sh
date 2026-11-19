#!/bin/bash
# w1 queue A2 (after A): EE reference-basis points (EE-own helicity basis), x-resolution test (J = 24 mid nodes), new even J (EE + C3H)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f w1_queueA.sh > /dev/null; do sleep 20; done
W1_TAG=_a python3 w1_rows.py refb 18,20,16,14,12,10,8,6,4,22,26,30
W1_TAG=_a python3 w1_rows.py even EE,C3H 24 1600 0.0625,0.0875,0.1125,0.1375,0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775
W1_TAG=_a python3 w1_rows.py even EE,C3H 22,26,30 1600
W1_TAG=_a python3 w1_rows.py even EE,C3H 24,28 800
W1_TAG=_a python3 w1_rows.py refb 34,36,38,40
W1_TAG=_a python3 w1_rows.py even EE,C3H 36,40 1600
W1_TAG=_a python3 w1_rows.py even EE,C3H 22,26 800
W1_TAG=_a python3 w1_rows.py even EE,C3H 34,38 1600
