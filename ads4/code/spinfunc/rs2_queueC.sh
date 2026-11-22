#!/bin/bash
# HRS_CLOSE: K = 5 engine images (a^5 term) for the J = 0 gap states
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 rs2_engine5.py 300 0 > cache/rs2_queueC_D300.log 2>&1
python3 rs2_engine5.py 1000 0 > cache/rs2_queueC_D1000.log 2>&1
