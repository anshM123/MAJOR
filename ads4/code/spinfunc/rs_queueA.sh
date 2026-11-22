#!/bin/bash
# HRS_CHECK queue A: n3 overlap nodes (Delta 1000, 300), then n4 below-fold nodes at Delta 300, then Delta = 100 (engine + nodes).
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 8.02,10.02,13.02,16.02,20.02 1000
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 8.02,10.02,13.02,16.02,20.02 300
RS_TAG=_n4 RS_NUSPLIT=9.5 Z2_M=9 Z2_R=8 python3 rs_points.py 0 8.02,7.02,6.02,5.02,4.02,3.02,2.02 300
