#!/bin/bash
# DIRCV.md §8: new window nodes at (1100,100) in nu = 11-15 (avoiding the D_h zero at nu = 13.09) + the PonC twin test at 13.3
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
(PW_TAG=_z1 PW_SEED=9101 PW_R=16 python3 dc_points.py 100 1100 11.0,12.3,15.0 > cache/run_z1.log 2>&1) &
(python3 dc_ponc.py 100 1100 13.3 8 4242 _c0 noref > cache/run_c0.log 2>&1;
 PW_TAG=_z2 PW_SEED=9201 PW_R=16 python3 dc_points.py 100 1100 11.7,14.2 > cache/run_z2.log 2>&1) &
wait
