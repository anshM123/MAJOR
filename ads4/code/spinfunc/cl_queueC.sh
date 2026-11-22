#!/bin/bash
# cl_queueC (RESULT_D_CLEANUP.md, task 3): direct n1 nodes on (1048, 48) at K2-grid x (consistency check of the fixed-engine K2 row
# used by FOLD_WINDOW_LARGEJ §4).  Starts after cl_queueB.sh.   -> cache/cl_pts_<kind>_J48_D1048.pkl
cd "$(dirname "$0")"
while pgrep -f "cl_queueB[.]sh" > /dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_MREF=11 HC_M=12 CL_FWPATCH=1
ulimit -v 3000000
python3 cl_points.py 48 1048 0.2,0.4,0.6,0.75 > cache/cl_qC_J48_1048.log 2>&1
echo "QUEUE CL-C DONE"
