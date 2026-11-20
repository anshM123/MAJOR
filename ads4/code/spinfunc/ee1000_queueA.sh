#!/bin/bash
# ee1000 queue A (EE_DG1000_LARGEJ.md): even-J K2 rows at the state's own near-gap mass (Delta = 1000 + J), 39 x nodes (0.05 ... 0.9).
# Resumable: re-run the script; existing points are skipped.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 EE1000_Q=A
python3 ee1000_rows.py even EE,C3H 40 1040 X39
python3 ee1000_rows.py refb 48,56,64
python3 ee1000_rows.py even EE,C3H 36 1036 X39
python3 ee1000_rows.py even EE,C3H 48 1048 X39
python3 ee1000_rows.py even EE,C3H 64 1064 X39
python3 ee1000_rows.py even EE,C3H 56 1056 X39
python3 ee1000_rows.py even EE,C3H 32 1040 HI
echo QUEUE_A_DONE
