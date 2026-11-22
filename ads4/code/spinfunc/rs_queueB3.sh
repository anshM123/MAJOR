#!/bin/bash
# HRS_CHECK queue B3: after the Delta = 1000 n4 job, the Delta = 100 sensitivity state (engine K = 4 + crystal nodes; the dG error
# scales like m^2/Delta, so Delta = 100 is ~3x more sensitive than 300 at the same cost).
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
while pgrep -f "rs_points.py 0 3.02,2.02,4.02,5.02,6.02 1000" > /dev/null; do sleep 30; done
python3 rs_engine.py 100 0
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 10.02,16.02,20.02 100
RS_TAG=_n4 RS_NUSPLIT=9.5 Z2_M=9 Z2_R=8 python3 rs_points.py 0 3.02,5.02,2.02 100
