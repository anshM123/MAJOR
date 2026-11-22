#!/bin/bash
# cl_driver (RESULT_D_CLEANUP.md): runs the chain analyses as soon as the direct points they need are finished (cl_queueA: small-nu,
# cl_queueB: edge).  Each analysis runs with the edge worker paused (cl_pause.sh), so at most 2 of this agent's processes are busy.
cd "$(dirname "$0")"
wA() { until grep -q "^# done" cache/cl_qA_J$1_$2.log 2>/dev/null; do sleep 30; done; }
wB() { until grep -q "^# done" cache/cl_qB_J$1_$2.log 2>/dev/null; do sleep 30; done; }
# Delta_g = 300
wB 4 400; ./cl_teven.sh 4 400 edge
wB 6 400; ./cl_teven.sh 6 400 edge
wA 8 400; wB 8 400; ./cl_teven.sh 8 400 oo edge
wA 9 400; wB 9 400; wB 9 200; ./cl_todd.sh 9 400 400,200 small edge
wA 10 400; wB 10 400; ./cl_teven.sh 10 400 oo edge
wA 11 400; wB 11 400; ./cl_todd.sh 11 400 400 small edge
wB 12 400; ./cl_teven.sh 12 400 edge
wB 15 400; wB 15 336; ./cl_todd.sh 15 336,400 336,400 small edge
wB 14 400; ./cl_teven.sh 14 400 edge
wA 7 400; wB 7 400; wB 7 200; ./cl_todd.sh 7 400 400,200 small edge
wA 13 400; wB 13 400; ./cl_todd.sh 13 400 400 small edge
echo "CL-DRIVER 300 DONE"
