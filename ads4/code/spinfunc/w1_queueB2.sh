#!/bin/bash
# w1 queue B2 (after B): direct EE RQMC cross-checks; J = 24 at the near-gap mass (Delta = 336: K2 EE + C3H, then direct EE + C3H);
# odd J = 15, 17, 19 (EE + C3H); Delta = 3200 rows for large J
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
while pgrep -f w1_queueB.sh > /dev/null; do sleep 20; done
python3 w3_direct.py run 18 336 0.075,0.1,0.125,0.15,0.175,0.2
W1_TAG=_b python3 w1_rows.py even EE,C3H 24 336
python3 w3_direct.py run 24 336 0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4
python3 w3_direct.py run 20 1600 0.05,0.075,0.1,0.15,0.2
W1_TAG=_b python3 w1_rows.py odd EE,C3H 15,17,19 400,800,1600
python3 w3_direct.py run 16 400 0.075,0.1,0.125,0.15,0.175,0.2
W1_TAG=_b python3 w1_rows.py refb 48,56
W1_TAG=_b python3 w1_rows.py even EE,C3H 40,48,32,56 3200
