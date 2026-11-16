"""Referee: is the ++ island of the leading Regge action hit by PHYSICAL states (integer even J)?  And the margin of F_R just
above m = 7.36 Lam (vs. the size Lam^2/m^2 of the omitted O(nu^2/m^2) terms).  Uses adsfunc/common.py (validated in t1 and
against the referee's mpmath evaluation ref_check.py island)."""
import sys
sys.path.insert(0, '../adsfunc')
import numpy as np
from common import regge_parts, kin, tau_of_m, lam_of_gap

for Dg in [50, 100, 1000]:
    Lam = lam_of_gap(Dg)
    worst = (np.inf, None)
    for J in range(16, 40, 2):
        for mfl in np.arange(5.0, 8.0, 0.01):
            m = mfl * Lam
            tau = tau_of_m(m, J)
            if tau < Dg:
                continue
            _, beta = kin(tau, J)
            B2, B3, B4, B5, _ = regge_parts('++', Lam, m, beta)
            F = B2 + B3 + B4 + B5
            r = F / (abs(B2) + abs(B3) + abs(B4) + abs(B5))
            if r < worst[0]:
                worst = (r, (J, round(mfl, 2), F))
    print('Dg=%d ++ physical even J: min rel F_R = %.4e at (J, m/Lam, F_R) = %s' % (Dg, worst[0], worst[1]), flush=True)
    # margin just above 7.36 (continuous beta, as in t2)
    mins = []
    for kap in [7.36, 8.0, 9.0, 10.0]:
        mn = np.inf
        for Lb in np.linspace(0.3, 40, 800):
            beta = Lb / Lam
            m = kap * Lam
            if m * np.exp(-beta / 2) < Dg - 2:
                continue
            B2, B3, B4, B5, _ = regge_parts('++', Lam, m, beta)
            F = B2 + B3 + B4 + B5
            mn = min(mn, F / (abs(B2) + abs(B3) + abs(B4) + abs(B5)))
        mins.append((kap, mn, 1 / kap ** 2))
    print('   min rel F_R at m/Lam = kappa (++):', ', '.join('k=%.2f: %.2e (Lam^2/m^2=%.3f)' % x for x in mins), flush=True)
