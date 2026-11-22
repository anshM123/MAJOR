#!/bin/bash
# qf queue D (after the J = 6 gap image of queue B): J = 6 at 1.475 Lam (K = 3), then even J = 8 (reference, Dg = 1000 gap, 1.475 Lam @ Dg = 300; K = 3).
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ulimit -v 3000000
run() { python3 qf_eeq.py "$@" >> cache/qf_eeq_$1_D$2_J$3.log 2>&1; }
# chain validation on BELOW_FOLD3's J = 4 images (must reproduce s_d(+-) = -0.516, s_d(++) = +0.517 at Dg = 300), then the J = 6 gap chain
python3 qf_qan.py even 304.0297506031432 4 4 300 > cache/qf_qan_even_J4val.log 2>&1
python3 qf_qan.py even 306.0297506031432 6 4 300 > cache/qf_qan_even_J6gap.log 2>&1
QF_XSMIN=0.25 python3 qf_qan.py even 306.0297506031432 6 4 300 > cache/qf_qan_even_J6gap_xs.log 2>&1
run even 441.8348597942253 6 3 4
run even 1006.0089775159122 6 4 5
L1_DEG=6 python3 qf_helbasis.py 6 1006.0089775159122 >> cache/qf_helbasis.log 2>&1
# (J = 8 reference and Dg = 1000 gap moved to qf_queueE.sh)
while pgrep -f qf_queueE.sh > /dev/null; do sleep 30; done      # keep total RSS < 3 GB: J = 8 images one at a time
run even 441.86892346077286 8 3 4
L1_DEG=6 python3 qf_helbasis.py 8 441.86892346077286 >> cache/qf_helbasis.log 2>&1
