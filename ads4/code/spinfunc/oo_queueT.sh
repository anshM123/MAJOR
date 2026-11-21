#!/bin/bash
# nu < 8.3 window at Delta_g = 300 (HYPOTHESES_CLOSE2.md §3): extra n4 thimble nodes (A+ + A-, CRN state/reference, z2_points unchanged),
# tagged caches (z2<tag>_pts_*): J = 32 nu = 4.994 (x = 0.015; first thimble node at J = 32), J = 28 nu = 4.999 two new seeds.
# Waits for oo_queueK.sh.  Resumable.
cd "$(dirname "$0")"
while pgrep -f "oo_queueK[.]sh" > /dev/null; do sleep 60; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Z2_M=10
Z2_TAG=oo Z2_SEED=4242 python3 z2_points.py 32 4.9937611076 336
Z2_TAG=oo3 Z2_SEED=6161 python3 z2_points.py 28 4.999254944 336
Z2_TAG=oo4 Z2_SEED=7171 python3 z2_points.py 28 4.999254944 336
echo "QUEUE OO-T DONE"
