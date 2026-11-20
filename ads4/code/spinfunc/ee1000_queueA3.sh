#!/bin/bash
# ee1000 queue A3 (replaces the rest of A2): even-J K2 rows at Delta = 1000 + J, 35 x nodes (0.05 ... 0.8; the LO/HI nodes of A2 were
# dropped: K2 is not reliable at nu <~ 50 nor at x > 0.8 for J >= 40, EE_DG1000_LARGEJ §3). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
while pgrep -f "ee1000_rows.py even EE,C3H 40 1040" > /dev/null; do sleep 20; done
python3 ee1000_rows.py refb 48,56,64
python3 ee1000_rows.py even EE,C3H 36 1036 X35
python3 ee1000_rows.py even EE,C3H 48 1048 X35
python3 ee1000_rows.py even EE,C3H 64 1064 X35
python3 ee1000_rows.py even EE,C3H 56 1056 X35
echo QUEUE_A3_DONE
