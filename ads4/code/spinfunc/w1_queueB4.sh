#!/bin/bash
# w1 queue B4 (replaces the rest of B2/B3): near-gap-mass rows at Delta_g = 300 for J >= 24 (where J = 24 fails), intermediate masses,
# odd J = 15, 19; then Delta = 3200 rows
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f "w3_direct.py run 20 1600" > /dev/null; do sleep 20; done
W1_TAG=_b python3 w1_rows.py even EE,C3H 28 336
W1_TAG=_b python3 w1_rows.py even EE,C3H 24,28 500
W1_TAG=_b python3 w1_rows.py odd EE,C3H 15,19 336,400
W1_TAG=_b python3 w1_rows.py even EE,C3H 32 336
W1_TAG=_b python3 w1_rows.py refb 48
W1_TAG=_b python3 w1_rows.py even EE,C3H 40,48 3200
