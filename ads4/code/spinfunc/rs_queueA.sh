#!/bin/bash
# HRS_CHECK queue A (restart after the shared-reference fix): n4 below-fold nodes at Delta = 300 (8.02 kept from the first run).
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
RS_TAG=_n4 RS_NUSPLIT=9.5 Z2_M=9 Z2_R=8 python3 rs_points.py 0 3.02,2.02,4.02,5.02,6.02,7.02 300
