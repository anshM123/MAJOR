#!/bin/bash
# HWIN_THIMBLE production queue A (resumable: finished nodes are skipped)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 PW_KINDS=EE,C3H PW_ADAPT=1 PW_TAG=_w PW_M=9 PW_R=8
python3 pw_points.py 64 1064 13.3,20.8,10,25
python3 pw_points.py 64 364 10,11.7,13.3
