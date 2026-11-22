#!/bin/bash
# cl_driver4 (RESULT_D_CLEANUP.md): re-run of J = 18 (Delta_g = 300) edge after the failed-node rule (CL_MAXREL) was added (the (336, 18)
# x = 0.9 node has 121 % RQMC error).  Starts after cl_driver3.sh.
cd "$(dirname "$0")"
while pgrep -f "cl_driver3[.]sh" > /dev/null; do sleep 30; done
CL_Z3SRC=hc ./cl_teven.sh 18 336,400 edge
echo "CL-DRIVER4 DONE"
