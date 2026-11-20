#!/bin/bash
# z2 production queue B: n4 (A+ + A- thimbles, below the fold) nodes on the Delta = 336 rows, J = 24 then 28, Delta_ref = 2e4 (CRN)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 Z2_M=9
python3 z2_points.py 24 10.008046762,8.340038969,6.672031175,5.004023381,4.003218705 336
python3 z2_points.py 28 9.998509889,8.332091574,6.665673259,4.999254944,3.999403956 336
