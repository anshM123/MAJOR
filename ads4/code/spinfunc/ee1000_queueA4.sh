#!/bin/bash
# ee1000 queue A4 (replaces the rest of A3): K2 is not converged at J = 48 near the gap (|M4| ~ |M2| ~ 2 % of the envelope, J^2/m = 2.2;
# EE_DG1000_LARGEJ §3), so the even J >= 56 K2 rows were dropped; instead a direct (RQMC) state sub-row at (1048, 48), EE and C3H,
# x = 0.1 ... 0.7 (K2 reference). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
for xs in 0.3,0.4,0.5 0.6,0.2,0.35 0.45,0.55,0.25 0.15,0.65,0.7 0.1,0.125,0.175 0.225,0.275,0.325 0.375,0.425,0.475 0.525,0.575,0.625 0.675; do
  python3 ee1000_direct.py 48 1048 $xs EE
  python3 ee1000_direct.py 48 1048 $xs C3H
done
echo QUEUE_A4_DONE
