#!/bin/bash
# w1 queue A (EE_NEARGAP.md): EE even rows at the q^2 keys -- J = 24, 28, 32 on Delta = 1600 (C3H exists), then J = 12, 14, 16
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
W1_TAG=_a python3 w1_rows.py even EE 24,28,32 1600 keys
W1_TAG=_a python3 w1_rows.py refb 24,28,32
W1_TAG=_a python3 w1_rows.py even EE 16,14,12 400,800,1600 keys
