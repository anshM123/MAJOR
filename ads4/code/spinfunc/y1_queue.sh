#!/bin/bash
# EE_DG100.md: K2 EE rows at Delta_g = 100 (worker 2)
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 y1_rows.py even 4,6,8,10 100
python3 y1_rows.py odd 5,7,9 100
python3 y1_rows.py even 10 200
python3 y1_rows.py even 0,2 100
python3 y1_rows.py odd 5,7,9 xrows
python3 y1_rows.py even 4,6,8,10,0,2 400
