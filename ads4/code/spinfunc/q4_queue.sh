#!/bin/bash
# q4 direct (RQMC) validation points on the Delta = 400 row (S1_S4_FIXES.md §2.2); one process, one thread.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
python3 q4_direct.py run C2 20 400 0.05,0.1,0.15,0.2,0.25,0.3
python3 q4_direct.py run C2 16 400 0.05,0.1,0.15,0.2,0.25,0.3
python3 q4_direct.py run C2 20 400 0.4,0.5,0.6,0.7,0.8
python3 q4_direct.py run C3H 20 400 0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8
python3 q4_direct.py run C2 16 400 0.4,0.5,0.6,0.7,0.8
python3 q4_direct.py run C3H 16 400 0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8
