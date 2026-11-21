#!/bin/bash
# HWIN_THIMBLE queue C (after A): scalar-only exact-scalar validation nodes, independent seed (resumable)
cd "$(dirname "$0")"
while pgrep -f pw_queueA.sh > /dev/null; do sleep 30; done
export OMP_NUM_THREADS=1 PW_ADAPT=1 PW_TAG=_v PW_M=10 PW_R=8 PW_SEED=5151 PW_MODES="1,-1"
python3 pw_points.py 64 1064 13.3,33.4,50
python3 pw_points.py 40 1040 13.3,20.8,30
python3 pw_points.py 100 1100 25,40
