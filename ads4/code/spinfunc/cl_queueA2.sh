#!/bin/bash
# cl_queueA2 (RESULT_D_CLEANUP.md, task 1): odd J = 9, 5, 7 -- the near-gap states at Delta_g = 300 (m ~ 1.02-1.1 Lam) lie between the
# Delta = 200 and 400 rows, so the 200 row's small-x nodes (K2 at x = 0.05, nu ~ 10) also enter; direct n3 nodes nu ~ 8.3-20 there.
# Starts after cl_queueA.sh.  Resumable.
cd "$(dirname "$0")"
while pgrep -f "cl_queueA[.]sh" > /dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Z2_NUSPLIT=5
ulimit -v 3000000
X200=0.042,0.05,0.059,0.067,0.084,0.1
for J in 9 5 7; do CL_MODE=hc HC_METHOD=n3 python3 cl_points.py $J 200 $X200 > cache/cl_qA_J${J}_200.log 2>&1; done
echo "QUEUE CL-A2 DONE"
