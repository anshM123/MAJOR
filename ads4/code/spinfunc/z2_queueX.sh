#!/bin/bash
# z2 queue X: n3 contour (nu_f = 13) at the below-fold nodes nu ~ 10, 8.3 (cross-check of / replacement for n4; tagged cache z2x_*)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 Z2_TAG=x Z2_NUSPLIT=8 Z2_M=11
while pgrep -f "z2_points.py 28 9.998509889 336" > /dev/null; do sleep 20; done
python3 z2_points.py 28 8.332091574 336
