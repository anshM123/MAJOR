#!/bin/bash
# z5 queue C2 (replaces the tail of z5_queueC.sh): after the J = 32 K2 row (w1_rows, tag _z5) finishes, the J = 32 direct nodes on the
# Delta = 336 row: n3 above the fold, then n3 (nu_f = 13) at nu ~ 10, 8.3 (tag x, as for J = 24, 28); thimble nodes nu < 8 not run.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
while pgrep -f "w1_rows.py even EE,C3H 32 336" > /dev/null; do sleep 30; done
Z2_M=11 python3 z2_points.py 32 13.316696287,16.645870359,11.652109251,20.807337948,24.968805538,33.291740717 336
Z2_TAG=x Z2_NUSPLIT=8 Z2_M=11 python3 z2_points.py 32 9.987522215,8.322935179 336
