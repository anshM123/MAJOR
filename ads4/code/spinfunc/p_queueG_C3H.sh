#!/bin/bash
# PARITY_ODD_ROUTEA.md sect. 11: K3 (L1_DEG=8) Laplace values at x <= 0.25 on the Delta_g = 100 rows (J = 9 at Delta = 100 into p3's
# own _k3 cache; extra rows m = 1.05 (J = 5) / 1.1 (J = 7, 9) Lam into cache/p8_k3_*), kind C3H; after p_queueE.sh. Resumable.
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=8
while pgrep -f "p_queueE.sh" > /dev/null; do sleep 30; done
P8_MFL=1.05 P8K_XS=0.05,0.1,0.15 python3 p8_k2rows.py C3H 5 >> cache/p8_k3rows_C3H.log 2>&1
P8_MFL=1.1 P8K_XS=0.05,0.1,0.15 python3 p8_k2rows.py C3H 7 9 >> cache/p8_k3rows_C3H.log 2>&1
P3_DS=100 P3_XS=0.05,0.1,0.15 python3 p3_run.py C3H 9 >> cache/p3_run_k3_C3H.log 2>&1
