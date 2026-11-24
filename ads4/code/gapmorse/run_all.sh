#!/bin/sh
# GAP_MORSE reproduction (<= 2 processes; outputs in cache/).  Taylor models of degree 8 on Hermitian balls.
ulimit -v 3000000
export OMP_NUM_THREADS=1
cd "$(dirname "$0")"
# 1. certified Taylor-model data, scalar k=2 integrand with Q_J (phase data are J-independent)
GC_DEG=8 python3 gm_run.py 0,2,4,5 0.5,0.3,0.7,0.85 0.03,0.04,0.05,0.06,0.07 > cache/run_a.log 2>&1 &
GC_DEG=8 python3 gm_run.py 6,7,8 0.5,0.3,0.7,0.85 0.03,0.04,0.05,0.06,0.07 > cache/run_b.log 2>&1
wait
GC_DEG=8 python3 gm_run.py 0,2,4,5 0.1 0.02,0.025,0.03,0.035,0.04 > cache/run_c.log 2>&1 &
GC_DEG=8 python3 gm_run.py 6,7,8 0.1 0.02,0.025,0.03,0.035,0.04 > cache/run_d.log 2>&1
wait
# 2. Delta_1 table (Lemma M)
python3 gm_table.py > cache/run_table.log 2>&1
