#!/bin/bash
# q2 queue C: mid x-nodes (0.35 ... 0.75) for J = 12-20 -- the 0.1 spacing of the l9 grid under-resolves the P_J oscillation
# (period in x ~ 0.18 at J = 20) -- Delta = 400 first (Delta_g = 300), then 800, 1600 (Delta_g = 1000).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Q2_TAG=_d python3 q2_run.py even 20,18,16,14,12 400 0.35,0.45,0.55,0.65,0.75
Q2_TAG=_d python3 q2_run.py even 20,18,16,14,12 800,1600 0.35,0.45,0.55,0.65,0.75
