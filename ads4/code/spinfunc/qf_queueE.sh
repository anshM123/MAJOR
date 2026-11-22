#!/bin/bash
# qf queue E (H_fwd^Q, even J = 8): reference and Dg = 1000 gap (K = 3), helicity bases.
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
run even 20000 8 2 3
run even 1008.0089775159122 8 3 4
L1_DEG=6 python3 qf_helbasis.py 8 20000 1008.0089775159122 >> cache/qf_helbasis.log 2>&1
