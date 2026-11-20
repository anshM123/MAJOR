#!/bin/bash
# Delta_g = 300 small-nu direct nodes (z4 protocol; n3 hybrid contour, nu_f = 13 below 13, CRN state/reference, M = 10) for the
# smallest-margin J <= 22 rows; waits for hc_queueOO.sh.  Resumable.
cd "$(dirname "$0")"
while pgrep -f hc_queueOO.sh > /dev/null; do sleep 60; done
export OMP_NUM_THREADS=1 HC_METHOD=n3 HC_M=10 Z2_NUSPLIT=5
python3 hc_points.py 22 336 0.025,0.03,0.035,0.04,0.05,0.0625,0.075,0.1
python3 hc_points.py 20 400 0.021,0.025,0.0293,0.0335,0.042,0.05,0.075,0.1
python3 hc_points.py 18 336 0.025,0.03,0.035,0.04,0.05,0.0625,0.075,0.1
echo "QUEUE B DONE"
