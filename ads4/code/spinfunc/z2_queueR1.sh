#!/bin/bash
# z2 resume queue R1 (after the container restart): J = 24 remaining n3 nodes, then J = 24 n4 nodes (z2 skips finished nodes)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
Z2_M=11 python3 z2_points.py 24 11.676054556,20.850097422,25.020116906,33.360155875 336
Z2_M=9 python3 z2_points.py 24 10.008046762,8.340038969,6.672031175,5.004023381 336
