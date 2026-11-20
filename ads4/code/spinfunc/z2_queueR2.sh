#!/bin/bash
# z2 resume queue R2: J = 28 n3 nodes, then J = 28 n4 nodes
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
Z2_M=11 python3 z2_points.py 28 13.331346519,16.664183148,11.664928204,20.830228935,24.996274722,33.328366297 336
Z2_M=9 python3 z2_points.py 28 9.998509889,8.332091574,6.665673259,4.999254944 336
