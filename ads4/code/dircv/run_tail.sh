#!/bin/bash
# DIRCV.md §3.7: per-point samples at (1100, 100), nu = 13.3 (two seeds, two processes)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 dc_tail.py 100 1100 13.3 9 4 7000 > cache/tail_7000.log 2>&1 &
python3 dc_tail.py 100 1100 13.3 9 4 8000 > cache/tail_8000.log 2>&1 &
wait
