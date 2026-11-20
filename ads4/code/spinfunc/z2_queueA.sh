#!/bin/bash
# z2 production queue A: n3 (above the fold) nodes on the Delta = 336 rows, J = 24 then 28, with Delta_ref = 2e4 (CRN)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 Z2_M=11
python3 z2_points.py 24 13.344062350,16.680077937,11.676054556,20.850097422,25.020116906,33.360155875 336
python3 z2_points.py 28 13.331346519,16.664183148,11.664928204,20.830228935,24.996274722,33.328366297 336
