#!/bin/bash
# z2 queue B2: J = 24 n4 nodes (queue B's J = 24 run failed on the nu_c rule; fixed), after queue B finishes
cd "$(dirname "$0")"
while kill -0 32612 2>/dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 Z2_M=9
python3 z2_points.py 24 10.008046762,8.340038969,6.672031175,5.004023381,4.003218705 336
