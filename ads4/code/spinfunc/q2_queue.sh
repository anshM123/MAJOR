#!/bin/bash
# q2 production queues (S1_S4_FIXES.md); two processes max, one thread each.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
if [ "$1" = A ]; then
  Q2_TAG=_a python3 q2_run.py even 12,14,16,18,20 400,800,1600
elif [ "$1" = B ]; then
  Q2_TAG=_b python3 q2_run.py odd 11,13 400,800,1600
  Q2_TAG=_b python3 q2_run.py even 12,14,16,18,20 200
  Q2_TAG=_b python3 q2_run.py odd 11,13 200
  Q2_TAG=_b python3 q2_run.py odd 15,17,19 400,800,1600
fi
