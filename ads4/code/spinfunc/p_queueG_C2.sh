#!/bin/bash
# PARITY_ODD_ROUTEA.md sect. 11: K3 (L1_DEG=8) Laplace values at x <= 0.25 on the Delta_g = 100 rows (J = 9 at Delta = 100 into p3's
# own _k3 cache; extra rows m = 1.05 (J = 5) / 1.1 (J = 7, 9) Lam into cache/p8_k3_*), kind C2; after p_queueD.sh. Resumable.
# Then the p_queueF.sh K3 validation points (Delta = 400; run p_queueF.sh body here so that <= 2 cores are used).
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=8
while pgrep -f "p_queueD.sh" > /dev/null; do sleep 30; done
P8_MFL=1.05 P8K_XS=0.05,0.1,0.15,0.2,0.25 python3 p8_k2rows.py C2 5 >> cache/p8_k3rows_C2.log 2>&1
P8_MFL=1.1 P8K_XS=0.05,0.1,0.15,0.2,0.25 python3 p8_k2rows.py C2 7 9 >> cache/p8_k3rows_C2.log 2>&1
P3_DS=100 P3_XS=0.05,0.1,0.15,0.2,0.25 python3 p3_run.py C2 9 >> cache/p3_run_k3_C2.log 2>&1
python3 - >> cache/p8_k3check.log 2>&1 <<'PY'
import os, sys
sys.path.insert(0, '.')
import p3_odd_action as P
P._fn = lambda kind, J: os.path.join('cache', 'p8_k3check_%s_J%d.pkl' % (kind, J))
for (J, D, x, kind) in ((5, 400.0, 0.3, 'C2'), (7, 400.0, 0.4, 'C2'), (9, 400.0, 0.4, 'C3H')):
    P.get(kind, J, D, [x*P.mphys(D, J)])
print('done k3check', flush=True)
PY
