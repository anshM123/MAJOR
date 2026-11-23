#!/bin/bash
# DIRCV.md §7: 4x the P_J-twin statistics at nu = 13.3 (+168 replicas in two processes)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
python3 dc_ptwin.py 100 1100 13.3 84 5001 _y1 > cache/run_y1.log 2>&1 &
python3 dc_ptwin.py 100 1100 13.3 84 6001 _y2 > cache/run_y2.log 2>&1 &
wait
