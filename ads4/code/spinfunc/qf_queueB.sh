#!/bin/bash
# qf queue B (H_fwd^Q, even J = 6 2x2 images -> ++ and +-): reference, Dg = 300 gap, 1.475 Lam (+- J = 6 delta = 0.05 binding mass).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
run even 20000 6 2 3
run even 306.0297506031432 6 4 5
run even 441.8348597942253 6 4 5
