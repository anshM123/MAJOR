#!/bin/bash
# queue E: finer small-x nodes at Delta = 400 (K2) for J = 16, 18, 20, then direct RQMC (C2) at the small-x nodes, J = 16, 18, 20
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Q2_TAG=_e python3 q2_run.py even 20,18,16 400 0.075,0.125,0.175
python3 q4_direct.py run C2 20 400 0.075,0.125,0.175
python3 q4_direct.py run C2 18 400 0.1,0.125,0.15,0.075,0.175,0.2,0.25
python3 q4_direct.py run C2 16 400 0.075,0.125,0.175
