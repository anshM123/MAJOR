#!/bin/bash
# cl_queueC3 (RESULT_D_CLEANUP.md, task 3): second worker for the (1048, 48) direct sub-row, x from the top down (cl_queueC2 goes up;
# finished nodes are skipped by both, so they meet in the middle).  Started after the small-nu queue had finished.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_MREF=11 HC_M=10 CL_FWPATCH=1
ulimit -v 3000000
python3 cl_points.py 48 1048 0.8,0.775,0.75,0.725,0.7,0.675,0.65,0.625,0.6,0.575,0.55,0.525,0.5,0.475,0.45,0.425,0.4 > cache/cl_qC3_J48_1048.log 2>&1
echo "QUEUE CL-C3 DONE"
