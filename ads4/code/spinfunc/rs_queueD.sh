#!/bin/bash
# HRS_CHECK queue D: diagnostics of the n3-node deviation at nu ~ 11-14 (Delta = 300): k = 2 vs k = 4, contour (nu_f) dependence,
# independent Sobol seed.
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
RS_KINDS=EE4,EE RS_TAG=_n3k2 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,14.02,20.02 300
N3_NUF_LOW=14.5 RS_TAG=_n3f145 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,11.02,20.02 300
N3_NUF_LOW=12.4 RS_TAG=_n3f124 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 300
N3_NUF=17 RS_TAG=_n3f17 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 14.02,16.02,20.02 300
Z2_SEED=9191 RS_TAG=_n3s2 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 300
