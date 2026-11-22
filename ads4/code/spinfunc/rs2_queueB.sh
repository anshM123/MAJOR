#!/bin/bash
# HRS_CLOSE: bookkeeping for the J = 2 state (Delta_g = 300, 1.05 Lambda) and its reference
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 rs2_bookkeep.py 314.935 2 4 > cache/rs2_queueB_J2.log 2>&1
python3 rs2_bookkeep.py 20000 2 2 > cache/rs2_queueB_J2ref.log 2>&1
