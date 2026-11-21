#!/bin/bash
# after hc_queueA2.sh: complete direct sub-rows x = 0.125 ... 0.8 on (336, 28), (336, 32) (K2 beyond its criterion; the K2 error
# oscillates with the P_J lobes, so a partial replacement of the row is meaningless -- the whole [0.125, 0.8] is replaced).  Resumable.
cd "$(dirname "$0")"
while pgrep -f "hc_queueA2[.]sh" > /dev/null; do sleep 60; done
export OMP_NUM_THREADS=1 HC_METHOD=n1 HC_ADAPT=1 HC_M=12 HC_MREF=11
python3 hc_points.py 28 336 0.25,0.35,0.45,0.55,0.65,0.7,0.75,0.8,0.125,0.15,0.175
python3 hc_points.py 32 336 0.25,0.35,0.45,0.55,0.65,0.7,0.75,0.8,0.125,0.15,0.175
echo "QUEUE A3 DONE"
