#!/bin/bash
# bj queue 2 (BELOW_FOLD_SPINJ.md §4): item-3 states, depth 2, fw (patched) structure basis.
cd "$(dirname "$0")"
ulimit -v 3000000
export BJ_ADAPT=1 BJ_PRUNE=1e-12 BJ_KTOL=1e-6 OMP_NUM_THREADS=1 BJ_PATCH=1
python3 bj_eeq.py 1100 100 2 5,6 fw > cache/bj_eeq_D1100_J100_d2.log 2>&1
python3 bj_eeq.py 20000 100 2 4 fw > cache/bj_eeq_D20000_J100_d2.log 2>&1
python3 bj_eeq.py 364 64 2 5,6 fw > cache/bj_eeq_D364_J64_d2.log 2>&1
python3 bj_eeq.py 20000 64 2 4 fw > cache/bj_eeq_D20000_J64_d2.log 2>&1
