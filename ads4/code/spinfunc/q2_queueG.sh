#!/bin/bash
# queue G: J = 32, 36 on the Delta = 1600 row (Delta_g = 1000 crossing, §2.6)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 q2_refbasis.py 32
Q2_TAG=_f python3 q2_run.py even 32 1600 0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8
python3 q2_refbasis.py 36
Q2_TAG=_f python3 q2_run.py even 36 1600 0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8
