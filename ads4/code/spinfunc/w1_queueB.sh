#!/bin/bash
# w1 queue B: EE odd rows (p3 engine) at the q^2 keys -- J = 5, 7, 9 (rows 200 ... 1600), 11, 13 (400, 800); then even J <= 10 (400, 1600)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
W1_TAG=_b python3 w1_rows.py odd EE 5 400,200,800,1600 keys
W1_TAG=_b python3 w1_rows.py odd EE 7,9 400,200,800,1600 keys
W1_TAG=_b python3 w1_rows.py odd EE 11,13 400,800 keys
W1_TAG=_b python3 w1_rows.py even EE 10,8,6,4,2,0 400,1600 keys
