#!/bin/bash
# cl_driver (RESULT_D_CLEANUP.md): Delta_g = 300 analyses that need the small-nu nodes of cl_queueA (and the edge rows of cl_queueB);
# each analysis runs through cl_pause.sh (serialised, one worker paused).  J = 4, 5, 6 were run before (cl_teven.sh 4/6 edge, cl_todd.sh 5).
cd "$(dirname "$0")"
wA() { until grep -q "^# done" cache/cl_qA_J$1_$2.log 2>/dev/null; do sleep 30; done; }
wB() { until grep -q "^# done" cache/cl_qB_J$1_$2.log 2>/dev/null; do sleep 30; done; }
wA 8 400; wB 8 400; ./cl_teven.sh 8 400 oo edge
wA 9 400; wB 9 400; wB 9 200; ./cl_todd.sh 9 400 400,200 small edge
wA 10 400; wB 10 400; ./cl_teven.sh 10 400 oo edge
wA 11 400; wB 11 400; ./cl_todd.sh 11 400 400 small edge
wA 7 400; wB 7 400; wB 7 200; ./cl_todd.sh 7 400 400,200 small edge
wA 13 400; wB 13 400; ./cl_todd.sh 13 400 400 small edge
echo "CL-DRIVER 300 DONE"
