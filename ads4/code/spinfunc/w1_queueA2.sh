#!/bin/bash
# w1 queue A2 (after A): new even J on Delta = 1600 / 800 (EE + C3H), reference-basis points
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f w1_queueA.sh > /dev/null; do sleep 20; done
W1_TAG=_a python3 w1_rows.py refb 22,26,30,34,36,38,40
W1_TAG=_a python3 w1_rows.py even EE,C3H 22,26,30 1600
W1_TAG=_a python3 w1_rows.py even EE,C3H 24,28 800
W1_TAG=_a python3 w1_rows.py even EE,C3H 36,40 1600
W1_TAG=_a python3 w1_rows.py even EE,C3H 22,26 800
W1_TAG=_a python3 w1_rows.py even EE,C3H 34,38 1600
