#!/bin/bash
# cl_driver3 (RESULT_D_CLEANUP.md): re-run of J = 12 and 15 (Delta_g = 300) with the dense near-gap mass grid (the first runs used the
# sparse grid; J >= 11 has no exact-skeleton states).  Starts after cl_driver2.sh.
cd "$(dirname "$0")"
while pgrep -f "cl_driver2[.]sh" > /dev/null; do sleep 30; done
./cl_teven.sh 12 400 edge
./cl_todd.sh 15 336,400 336,400 small edge
echo "CL-DRIVER3 DONE"
