#!/bin/bash
# qf queue F: J = 8 at 1.475 Lam (Dg = 300), K = 3, in parallel with queue E (two J = 8 K = 3 images, ~1.2 GB each).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
run even 441.86892346077286 8 3 4
while pgrep -f qf_queueE.sh > /dev/null; do sleep 30; done
L1_DEG=6 python3 qf_helbasis.py 8 441.86892346077286 >> cache/qf_helbasis.log 2>&1
