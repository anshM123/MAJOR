#!/bin/bash
# cl_driver1000 (RESULT_D_CLEANUP.md): H_edge at Delta_g = 1000 for low J -- runs cl_t1000.sh once the edge rows of a J are finished
# (all computed edge rows of that J are injected; the gap states at m ~ 1000 sit between the 400/800 and 1600 rows).
cd "$(dirname "$0")"
wB() { until grep -q "^# done" cache/cl_qB_J$1_$2.log 2>/dev/null; do sleep 30; done; }
wB 4 1600; ./cl_t1000.sh 4 400,1600
wB 6 1600; ./cl_t1000.sh 6 400,1600
wB 7 1600; wB 7 800; ./cl_t1000.sh 7 200,400,800,1600
wB 8 1600; ./cl_t1000.sh 8 400,1600
wB 9 1600; wB 9 800; ./cl_t1000.sh 9 200,400,800,1600
wB 10 1600; ./cl_t1000.sh 10 400,1600
wB 11 800; ./cl_t1000.sh 11 400,800
wB 13 800; ./cl_t1000.sh 13 400,800
wB 12 1600; wB 12 800; ./cl_t1000.sh 12 400,800,1600
wB 14 1600; wB 14 800; ./cl_t1000.sh 14 400,800,1600
wB 16 1600; wB 16 800; ./cl_t1000.sh 16 400,800,1600
echo "CL-DRIVER 1000 DONE"
