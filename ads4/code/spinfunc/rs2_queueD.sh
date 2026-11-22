#!/bin/bash
# HRS_CLOSE: parity bound (one alpha-arc) and K = 5 images with v'-u', R kept to relative a^4 ('full')
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
while pgrep -f "rs2_bookkeep.py 20000 2" > /dev/null; do sleep 20; done
python3 rs2_parity.py 300 0 > cache/rs2_queueD_parity.log 2>&1
python3 rs2_engine5.py 300 0 full > cache/rs2_queueD_D300full.log 2>&1
python3 rs2_engine5.py 1000 0 full > cache/rs2_queueD_D1000full.log 2>&1
