#!/bin/bash
# t1 queue (RESULT_D_CONDITIONAL.md): extended band maps (fine quadrature, beta <= 2), one process at a time
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1
T1_DG=300 T1_BMAX=2.0 T1_NQ=10 T1_NL=40 T1_NB=70 T1_TAG=_ext python3 t1_band.py r5_lp_300_b0.12_dw0.1.pkl r5_lp_300_b0.2_dw0.1.pkl r5_lp_300_b0.25_dw0.0.pkl r5_lp_300_b0.3_dw0.0.pkl > cache/t1_queue_300.log 2>&1
T1_DG=1000 T1_BMAX=2.0 T1_NQ=10 T1_NL=40 T1_NB=70 T1_TAG=_ext python3 t1_band.py r5_lp_1000_b0.066_dw0.1.pkl r5_lp_1000_b0.12_dw0.1.pkl r5_lp_1000_b0.16_dw0.0.pkl > cache/t1_queue_1000.log 2>&1
T1_DG=100 T1_BMAX=2.0 T1_NQ=10 T1_NL=40 T1_NB=70 T1_TAG=_ext python3 t1_band.py > cache/t1_queue_100.log 2>&1
