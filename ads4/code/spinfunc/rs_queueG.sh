#!/bin/bash
# HRS_CHECK queue G: independent Sobol seeds at nu = 12.02 (with the 20.02 kappa anchor) for Delta = 1000 and 100: is the
# seed-4242 deviation at nu ~ 11-14 (dG proportional to Delta, i.e. log-ratio ~ 1/Delta) a sampling fluctuation or an a^1-type error?
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
if [ "$1" = "1000" ]; then
  Z2_SEED=9191 RS_TAG=_n3s2 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 1000
  Z2_SEED=5151 RS_TAG=_n3s3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 1000
else
  Z2_SEED=9191 RS_TAG=_n3s2 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 100
  Z2_SEED=5151 RS_TAG=_n3s3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 12.02,20.02 300
fi
