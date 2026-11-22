#!/bin/bash
# HRS_CHECK queue E (after D): the nu ~ 11-14 n3 nodes at Delta = 100 and 1000 (Delta scaling of the deviation seen at 300).
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
while pgrep -f rs_queueD.sh > /dev/null; do sleep 30; done
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 11.02,12.02,14.02,25.02 100
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 11.02,12.02,14.02,25.02 1000
