#!/bin/bash
# SMOOTH_BOUND.md §4: all corrected-chain runs (curvature double count), kappa image (default), stage post.  Sequential, 1 core.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
IMG=${SM_CURVIMG:-kappa}
( export Q3_JEVEN=0,2,4,6,8,10,12,14,16,18,20,22,24,28,32 Q3_JODD=5,7,9,11,13,15 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
  for v in row only; do SM_CURV=$v SM_TAG=${v}_$IMG nice python3 sm_chainfix.py w2 1000 300 > /dev/null 2>&1; done )
( export Z4_EST=best Z4_JS=24,28 Q3_JEVEN=24,28 Q3_JODD=5 Q3_EXTRA_M=1.02,1.05,1.1,1.15,1.2,1.3,1.4,1.5,1.75,2,2.5,3,4,5
  for v in row only; do SM_CURV=$v SM_JACFIX=0 SM_TAG=${v}_jac0_$IMG nice python3 sm_chainfix.py z4 300 > /dev/null 2>&1; done
  SM_CURV=row SM_JACFIX=1 SM_TAG=row_jac1_$IMG nice python3 sm_chainfix.py z4 300 > /dev/null 2>&1 )
( export EE1000_XMAX=0.8 EE1000_K=4 EE1000_MODEL=shape EE1000_XMIN=0.05 W2_DMIN=1000
  for v in row only; do
    EE1000_JEVEN=32,36,40 EE1000_JODD=17 SM_CURV=$v SM_TAG=f1_${v}_$IMG nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
    EE1000_JEVEN=40 EE1000_DSTATE=1 SM_CURV=$v SM_TAG=fd40_${v}_$IMG nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
    EE1000_JEVEN=48 EE1000_DSTATE=1 SM_CURV=$v SM_TAG=d48_${v}_$IMG nice python3 sm_chainfix.py ee1000 > /dev/null 2>&1
  done )
echo "# chainfix done ($IMG)"
