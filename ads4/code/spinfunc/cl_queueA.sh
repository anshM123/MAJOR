#!/bin/bash
# cl_queueA (RESULT_D_CLEANUP.md, task 1): direct small-nu n3 nodes at Delta_g = 300 on the Delta = 400 rows that carry the near-gap
# states, lowest-margin states first.  Even J: 3x3 (EE block + OO block, oo_points method); odd J: 1x1 odd integrand (hc_points n3).
# x = 0.021 ... 0.05 (nu ~ 8.3-20), 2^10*2*8 points, state + reference 2e4 on CRN.  Resumable (finished nodes skipped).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 Z2_NUSPLIT=5
ulimit -v 3000000
X400=0.021,0.025,0.0293,0.0335,0.042,0.05
CL_MODE=oo python3 cl_points.py 8 400 $X400  > cache/cl_qA_J8_400.log 2>&1
CL_MODE=hc HC_METHOD=n3 python3 cl_points.py 9 400 $X400  > cache/cl_qA_J9_400.log 2>&1
CL_MODE=oo python3 cl_points.py 10 400 $X400 > cache/cl_qA_J10_400.log 2>&1
CL_MODE=hc HC_METHOD=n3 python3 cl_points.py 11 400 $X400 > cache/cl_qA_J11_400.log 2>&1
CL_MODE=hc HC_METHOD=n3 python3 cl_points.py 7 400 $X400  > cache/cl_qA_J7_400.log 2>&1
CL_MODE=hc HC_METHOD=n3 python3 cl_points.py 13 400 $X400 > cache/cl_qA_J13_400.log 2>&1
echo "QUEUE CL-A DONE"
