#!/bin/bash
# EDGE_R2.md §2: K2 basis records + shadows (two lanes; <= 2 processes)
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
( python3 e2_k2.py refb 34,38,42,46,50 > cache/k2_refb.out 2>&1 ) &
( for JD in 26:1026 30:1030 34:1034 38:1038 42:1042 46:1046 50:1050; do
    python3 e2_k2.py state ${JD%%:*} ${JD##*:} 0.2,0.3,0.4,0.5,0.6,0.7,0.8
  done > cache/k2_state.out 2>&1
  python3 e2_k2.py oshadow 33 1033 >> cache/k2_state.out 2>&1; python3 e2_k2.py oshadow 33 20000 >> cache/k2_state.out 2>&1 ) &
wait
