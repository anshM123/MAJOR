#!/bin/bash
# H_OO / item 3 (HYPOTHESES_CLOSE2.md): J = 16 row 400, direct nodes x = 0.075, 0.1 (the K2 node at x = 0.075 is off the direct trend).
cd "$(dirname "$0")"
while pgrep -f "oo_queueD2[.]sh" > /dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OO_M=10 OO_R=8
python3 oo_points.py 16 400 0.075,0.1
echo "QUEUE OO-D3 DONE"
