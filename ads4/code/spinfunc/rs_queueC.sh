#!/bin/bash
# HRS_CHECK queue C: J = 2 (++ , m = 1.05 Lambda at Delta_g = 300: Delta = 314.935, the bg_eeq state) after queue A; then extra
# cheap n3 overlap nodes at J = 0.
cd "$(dirname "$0")"
ulimit -v 3000000
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 L1_DEG=6
while pgrep -f "rs_points.py 0 3.02,2.02,4.02,5.02,6.02,7.02 300" > /dev/null; do sleep 30; done
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 2 10.02,16.02,20.02 314.935
RS_TAG=_n4 RS_NUSPLIT=9.5 Z2_M=9 Z2_R=8 python3 rs_points.py 2 3.02,5.02 314.935
RS_TAG=_n3 RS_NUSPLIT=0 Z2_M=10 Z2_R=8 python3 rs_points.py 0 11.02,12.02,14.02,18.02,25.02 300
