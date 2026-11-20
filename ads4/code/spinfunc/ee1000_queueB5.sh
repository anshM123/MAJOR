#!/bin/bash
# ee1000 queue B5 (replaces the rest of B4): odd-J K2 rows at Delta = 1000 + J (35 nodes; J^2/m <= 0.93, K2 converged). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=B
while pgrep -f "ee1000_direct.py 64 1064" > /dev/null; do sleep 20; done
for J in 31 23 27 19; do python3 ee1000_rows.py odd EE,C3H $J $((1000 + J)) X35; done
echo QUEUE_B5_DONE
