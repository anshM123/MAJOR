#!/bin/bash
# PARITY_ODD_ROUTEA.md sect. 11: K3 (L1_DEG=8) Laplace values at the p8 validation points with the largest direct-K2 gaps (Delta = 400),
# separate cache cache/p8_k3check_<kind>_J<J>.pkl; runs after p_queueD.sh.
cd /home/user/MAJOR/ads4/code/spinfunc
export OMP_NUM_THREADS=1 L1_DEG=8
while pgrep -f "p_queueD.sh" > /dev/null; do sleep 30; done
python3 - >> cache/p8_k3check.log 2>&1 <<'PY'
import os, sys
sys.path.insert(0, '.')
import p3_odd_action as P
P._fn = lambda kind, J: os.path.join('cache', 'p8_k3check_%s_J%d.pkl' % (kind, J))
for (J, D, x, kind) in ((5, 400.0, 0.3, 'C2'), (7, 400.0, 0.4, 'C2'), (9, 400.0, 0.4, 'C3H')):
    P.get(kind, J, D, [x*P.mphys(D, J)])
print('done k3check', flush=True)
PY
