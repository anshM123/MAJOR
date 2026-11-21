#!/bin/bash
# nu < 8.3 window (HYPOTHESES_CLOSE2.md §3), continuation: an INDEPENDENT second seed of the J = 32 thimble node nu = 4.994 (tag oo2),
# then the J = 28 extra seeds (oo3, oo4).  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Z2_M=10
Z2_TAG=oo2 Z2_SEED=6161 python3 z2_points.py 32 4.9937611076 336
Z2_TAG=oo3 Z2_SEED=6161 python3 z2_points.py 28 4.999254944 336
Z2_TAG=oo4 Z2_SEED=7171 python3 z2_points.py 28 4.999254944 336
echo "QUEUE OO-T2 DONE"
