#!/bin/bash
# queue H: a row at the binding states' own mass for J = 18 (Delta = 336: m = 334 = 1.12 Lam(300)), K2 full grid, then direct RQMC C2
# at x = 0.075 ... 0.2 (removes the Delta = 400 -> m ~ 320-390 transfer for J = 18; §2.4 / §4 item 1)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Q2_TAG=_h python3 q2_run.py even 18 336 0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8
python3 q4_direct.py run C2 18 336 0.075,0.1,0.125,0.15,0.175,0.2
