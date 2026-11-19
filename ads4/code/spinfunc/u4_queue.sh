#!/bin/bash
# u4 queue (EE_FUNCTIONAL.md §4): sequential, 2 threads
export OMP_NUM_THREADS=2
cd "$(dirname "$0")"
for dg in 100 1000; do U4_DG=$dg U4_LAM=0 python3 u4_combo_lp.py check; done
for dg in 100 300 1000; do U4_DG=$dg U4_LAM=0 python3 u4_combo_lp.py lp; done
for dg in 300 1000; do U4_DG=$dg U4_LAM=-0.05,-0.02,-0.01,-0.005,0.005,0.01 python3 u4_combo_lp.py check; done
for dg in 300 1000; do U4_DG=$dg U4_LAM=-0.02,-0.01,-0.005,0.005,0.01 python3 u4_combo_lp.py lp; done
