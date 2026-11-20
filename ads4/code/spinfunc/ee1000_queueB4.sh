#!/bin/bash
# ee1000 queue B4 (replaces the rest of B3): a direct-RQMC sub-row at (1040, 40) on x = 0.2 ... 0.8 (state direct, K2 reference, as
# W2_DSTATE) to measure the K2 systematic in the integrated correction; direct points at (1064, 64); odd-J K2 rows (35 nodes). Resumable.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=B
while pgrep -f "ee1000_rows.py odd EE,C3H 17 1017" > /dev/null; do sleep 20; done
python3 ee1000_direct.py 40 1040 0.4,0.3,0.45,0.35,0.55,0.65,0.25,0.2,0.75,0.8,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.325,0.375,0.225,0.275
python3 ee1000_direct.py 64 1064 0.3,0.6,0.45
python3 ee1000_rows.py odd EE,C3H 31 1031 X35
python3 ee1000_rows.py odd EE,C3H 23 1023 X35
python3 ee1000_rows.py odd EE,C3H 27 1027 X35
python3 ee1000_rows.py odd EE,C3H 19 1019 X35
python3 ee1000_rows.py odd EE,C3H 21 1021 X35
python3 ee1000_rows.py odd EE,C3H 25 1025 X35
python3 ee1000_rows.py odd EE,C3H 29 1029 X35
echo QUEUE_B4_DONE
