#!/bin/bash
# HRS_CHECK queue F: sample-size dependence (heavy-tail test) of the k = 4 deviation at nu = 12.02, Delta = 300: 2^12 x 2 x 8 points,
# new Sobol seed (and the 20.02 kappa anchor); after queue B3 (Delta = 100 n4) has finished.
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
while pgrep -f rs_queueB3.sh > /dev/null; do sleep 30; done
Z2_SEED=5151 RS_KINDS=EE4,EE RS_TAG=_n3m12 RS_NUSPLIT=0 Z2_M=12 Z2_R=8 python3 rs_points.py 0 12.02,20.02 300
