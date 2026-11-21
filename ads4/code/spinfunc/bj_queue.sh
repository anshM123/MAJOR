#!/bin/bash
# bj queue (BELOW_FOLD_SPINJ.md): depth-2 engine runs, one at a time (the other core runs the depth-3 J=28 state).
cd "$(dirname "$0")"
ulimit -v 3000000
export BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 OMP_NUM_THREADS=1
run() { python3 bj_eeq.py "$@" > cache/bj_eeq_D$1_J$2_d$3.log 2>&1; }
run 336 32 2 8 win
run 20000 32 2 4 win
run 336 28 2 8 win
run 20000 28 2 4 win
BJ_PATCH=1 run 1100 100 2 12 fw
BJ_PATCH=1 run 20000 100 2 6 fw
BJ_PATCH=1 run 364 64 2 10 fw
BJ_PATCH=1 run 20000 64 2 6 fw
