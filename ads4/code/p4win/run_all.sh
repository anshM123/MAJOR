#!/bin/bash
# P4_WINDOW.md: reproduce everything (~0.6 CPU-h; at most 2 processes).  $1 = path to CHMRS ancillary regge_moments_order8.m (2106.10274).
set -e
cd "$(dirname "$0")"
ulimit -v 3000000; export OMP_NUM_THREADS=1
printf '*\n!.gitignore\n' > cache/.gitignore
python3 p4w_series.py "$1"                                   # §2.2  CHMRS series orders 1..4, Legendre-reduced; flat parts g_1..g_4
python3 p4w_exact.py > cache/exact.out                       # §3.1  exact scalar: P4 remainder, series vs exact, window splits
python3 p4w_check.py > cache/check.out & python3 p4w_316.py > cache/p316.out; wait   # §3.2-3.3 envelope + Lemma P4W on the exact scalar; §3.4
P4W_T1=15 python3 p4w_states.py > cache/states15.out & P4W_T1=10 python3 p4w_states.py > cache/states10.out; wait  # §4 (47 states)
P4W_T1=40 python3 p4w_states.py > cache/states40.out & python3 p4w_L.py > cache/L.out; wait                          # §4 sensitivity, §5
