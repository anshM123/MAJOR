#!/bin/bash
# qf queue C (H_fwd^Q, parity-odd +- at odd J = 7, 9, 11): references; J = 7 at K = 4 (both gaps), others K = 3 (depth 4):
# their failure thresholds are s_fail >= 8.9 (qf_tol), so the K = 3 truncation (~0.1-0.5) is ample.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
run odd 20000 7 2 3
run odd 1007.0089775159122 7 4 5
run odd 307.0297506031432 7 4 5
run odd 20000 9 2 3
run odd 1009.0089775159122 9 3 4
run odd 441.8507564998612 7 3 4
run odd 309.0297506031432 9 3 4
run odd 441.88936039600236 9 3 4
run odd 20000 11 2 3
run odd 1011.0089775159122 11 3 4
