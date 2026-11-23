#!/bin/bash
ulimit -v 3000000; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
while pgrep -f "cr_seed2.sh" > /dev/null; do sleep 10; done
python3 pcyc_farw.py > cache/log_farw.txt 2>&1
