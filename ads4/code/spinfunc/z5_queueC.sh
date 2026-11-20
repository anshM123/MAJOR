#!/bin/bash
# z5 queue C (J = 32 at Delta_g = 300): K2 EE + C3H row on Delta = 336 (19 nodes, w1_rows.py, tag _z5 -> cache/w1_pts_*_J32_k2_z5.pkl),
# then the z2 direct small-nu nodes (n3 above the fold, n4 below) on that row.  Started when a core is free (<= 2 processes).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
W1_TAG=_z5 python3 w1_rows.py even EE,C3H 32 336
Z2_M=11 python3 z2_points.py 32 13.316696287,16.645870359,11.652109251,20.807337948,24.968805538,33.291740717 336
Z2_M=9 python3 z2_points.py 32 9.987522215,8.322935179,6.658348143,4.993761108,3.995008886 336
