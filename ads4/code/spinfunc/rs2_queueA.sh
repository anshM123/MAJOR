#!/bin/bash
# HRS_CLOSE: bookkeeping queue (J = 0 states and references; J = 2 state)
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 rs2_bookkeep.py 1000 0 4 > cache/rs2_queueA_D1000.log 2>&1
python3 rs2_bookkeep.py 20000 0 4 > cache/rs2_queueA_D20000.log 2>&1
python3 rs2_bookkeep.py 100 0 4 > cache/rs2_queueA_D100.log 2>&1
