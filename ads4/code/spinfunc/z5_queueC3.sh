#!/bin/bash
# z5 queue C3 (after the second container restart): J = 32 direct nodes on the Delta = 336 row (K2 row done: w1_*_J32_k2_z5.pkl),
# lowest-nu n3 nodes first, then n3 (nu_f = 13) at nu ~ 10, 8.3 (tag x).  z2 skips finished nodes, so re-running resumes.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
Z2_M=11 python3 z2_points.py 32 11.652109251,13.316696287,16.645870359,20.807337948,24.968805538,33.291740717 336
Z2_TAG=x Z2_NUSPLIT=8 Z2_M=11 python3 z2_points.py 32 9.987522215,8.322935179 336
