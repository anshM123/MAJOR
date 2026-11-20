#!/bin/bash
# ee1000 queue R1 (after the 2nd container restart): direct (RQMC) state sub-row at (1048, 48), kind $1 (EE or C3H), x = 0.1 ... 0.7
# step 0.025 (K2 reference), then the odd-J K2 rows listed in $2 (35 nodes). Resumable; run one instance per kind (2 cores).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=R$1
for xs in 0.3,0.4,0.5 0.6,0.2,0.35 0.45,0.55,0.25 0.15,0.65,0.7 0.1,0.125,0.175 0.225,0.275,0.325 0.375,0.425,0.475 0.525,0.575,0.625 0.675; do
  python3 ee1000_direct.py 48 1048 $xs $1
done
for J in ${2//,/ }; do python3 ee1000_rows.py odd EE,C3H $J $((1000 + J)) X35; done
echo QUEUE_R1_$1_DONE
