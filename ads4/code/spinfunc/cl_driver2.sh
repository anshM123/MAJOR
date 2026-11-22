#!/bin/bash
# cl_driver2 (RESULT_D_CLEANUP.md): Delta_g = 300 H_edge analyses on rows whose small-nu nodes already exist (HYPOTHESES_CLOSE/2,
# PARITY_ODD_CLOSE): J = 12, 14 (3x3 oo points), 15 (odd, hc_pts), 18, 20, 22 (hc_pts).
cd "$(dirname "$0")"
wB() { until grep -q "^# done" cache/cl_qB_J$1_$2.log 2>/dev/null; do sleep 30; done; }
wB 12 400; ./cl_teven.sh 12 400 edge
wB 15 400; wB 15 336; ./cl_todd.sh 15 336,400 336,400 small edge
wB 14 400; ./cl_teven.sh 14 400 edge
wB 18 336; wB 18 400; CL_Z3SRC=hc ./cl_teven.sh 18 336,400 edge
wB 20 400; CL_Z3SRC=hc ./cl_teven.sh 20 400 edge
wB 22 336; CL_Z3SRC=hc ./cl_teven.sh 22 336 edge
echo "CL-DRIVER2 300 DONE"
