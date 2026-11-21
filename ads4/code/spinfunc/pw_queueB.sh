#!/bin/bash
# HWIN_THIMBLE production queue B (resumable)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 PW_KINDS=EE,C3H PW_ADAPT=1 PW_TAG=_w PW_M=9 PW_R=8
python3 pw_points.py 100 1100 16.3,13.3,20.8,25,10
python3 pw_points.py 100 400 10,13.3
