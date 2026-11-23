#!/bin/bash
# P2_CONTRACT: per-state exact contraction (2 processes)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 pc_run.py 0 24 > cache/run_a.out 2>&1 &
python3 pc_run.py 24 47 > cache/run_b.out 2>&1 &
wait
