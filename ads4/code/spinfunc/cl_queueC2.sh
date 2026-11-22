#!/bin/bash
# cl_queueC2 (RESULT_D_CLEANUP.md, task 3): the 4-node check on (1048, 48) found the fixed-engine K2 row off by 300-3100 c2 at x >= 0.4,
# so a complete direct n1 sub-row on the K2 grid x = 0.125 ... 0.8 (step 0.025) is computed, to inject into the J = 48 chain
# (HC_M=10 from the 4th node on, for time; partial replacement of an oscillating row is meaningless; HYPOTHESES_CLOSE §5).  fw structure basis (CL_FWPATCH=1).  Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CL_MODE=hc HC_METHOD=n1 HC_ADAPT=1 HC_MREF=11 HC_M=10 CL_FWPATCH=1
ulimit -v 3000000
XS=$(python3 -c "print(','.join('%g' % (0.125 + 0.025*i) for i in range(28)))")
python3 cl_points.py 48 1048 $XS > cache/cl_qC2_J48_1048.log 2>&1
echo "QUEUE CL-C2 DONE"
