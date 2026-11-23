#!/bin/bash
# DIRCV.md §4: new independent-seed nodes at (1100, 100), two processes (resumable; finished nodes are skipped).
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
(PW_TAG=_x1 PW_SEED=1001 PW_R=24 python3 dc_points.py 100 1100 13.3,16.3 > cache/run_x1.log 2>&1;
 PW_TAG=_x3 PW_SEED=3003 PW_R=16 python3 dc_points.py 100 1100 20.8,25 > cache/run_x3.log 2>&1) &
(PW_TAG=_x2 PW_SEED=2002 PW_R=24 python3 dc_points.py 100 1100 16.3,13.3 > cache/run_x2.log 2>&1;
 PW_TAG=_x4 PW_SEED=4004 PW_R=16 python3 dc_points.py 100 1100 25,20.8 > cache/run_x4.log 2>&1) &
wait
