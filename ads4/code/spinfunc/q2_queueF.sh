#!/bin/bash
# queue F: J = 24, 28 on the Delta = 1600 row only (c_off and Gamma_hh are m-independent at fixed x, §2.2(b)), for the Delta_g = 1000
# crossing estimate (§2.6); full node set; reference-basis points.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 q2_refbasis.py 24,28
Q2_TAG=_f python3 q2_run.py even 24,28 1600 0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8
