#!/bin/bash
# u4 queue 2 (EE_FUNCTIONAL.md §4): small-|lam| LPs, 1 thread
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
for dg in 300 1000; do U4_DG=$dg U4_LAM=-0.002,-0.001,0.001,0.002 python3 u4_combo_lp.py lp; done
