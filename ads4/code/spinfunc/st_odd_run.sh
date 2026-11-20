#!/bin/bash
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 W2_DMIN=1000 EE1000_JODD=17,19,21,23,25,27,29,31 EE1000_JEVEN=
B="EE1000_MODEL=shape EE1000_XMIN=0.05"
for v in row none only; do
  env EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=$v SM_TAG=fodd_${v}_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
done
# ensemble variants with the fix
env EE1000_MODEL=shape EE1000_K=5 EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=row SM_TAG=fodd_K5_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=shape EE1000_K=6 EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=row SM_TAG=fodd_K6_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.9 EE1000_XMIN=0.05 SM_CURV=row SM_TAG=fodd_x09_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=spline EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=row SM_TAG=fodd_spline_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0 SM_CURV=row SM_TAG=fodd_withlo_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 EE1000_NPN=64 SM_CURV=row SM_TAG=fodd_npn64_row_kappa nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
env EE1000_MODEL=shape EE1000_K=4 EE1000_XMAX=0.8 EE1000_XMIN=0.05 SM_CURV=row SM_CURVIMG=bessel SM_TAG=fodd_row_bessel nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
echo ODD_DONE
