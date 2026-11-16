"""Referee: true threshold kappa_0 of the ++ island of F_R (fine beta grid), Dg = 100, 1000."""
import sys
sys.path.insert(0, '../adsfunc')
import numpy as np
from common import regge_parts, lam_of_gap
for Dg in [100, 1000]:
    Lam = lam_of_gap(Dg)
    last_neg = None
    for kap in np.arange(7.30, 7.80, 0.02):
        mn = np.inf
        for Lb in np.linspace(7.0, 11.0, 401):
            B = regge_parts('++', Lam, kap * Lam, Lb / Lam)[:4]
            mn = min(mn, sum(B) / sum(abs(x) for x in B))
        if mn < 0:
            last_neg = (kap, mn)
    print('Dg=%d: largest kappa on 0.02 grid with min F_R<0: %s' % (Dg, last_neg), flush=True)
