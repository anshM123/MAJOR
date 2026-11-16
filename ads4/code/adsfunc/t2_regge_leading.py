"""t2: the leading Regge-limit AdS action of Omega_g3 (heavy formula C^pm_{k,nu} = (2P/m^{2k})(|c_s|^2 +- |c_o|^2)),
normal form  FR = B2 + B3 + B4 + B5  (see common.regge_parts), on a grid of (Lam, m/Lam, beta), both helicity channels.

Reported per (Lam, h):
  * kappa_0: smallest kappa on the grid such that FR > 0 for all m/Lam >= kappa and all beta on the grid;
  * the most negative point;
  * kappa(C0): the same with the envelope FR - C0 cosh(beta)/m^2 (|B2|+|B3|+|B4|+|B5|) (a MODEL of the unknown 1/m^2 terms).
The Regge formula is only valid for m >> Lam (nu <= Lam << m); points with m ~ Lam are shown for orientation only.
Only PHYSICAL points are used: twist tau = 2 + alpha >= Dg with alpha = m e^{-beta/2} (d = 3), i.e. beta <= 2 log(m/(Dg-2));
Lam = sqrt((Dg-2)(Dg-1)).  On physical points cosh(beta)/m^2 = (1/alpha^2 + 1/alphabar^2)/2 <= 1/(Dg-2)^2.
Floating point; not a certificate.
"""
import numpy as np
from common import curvT_all, SIGMA, lam_of_gap

import sys
DGS = [int(a) for a in sys.argv[1:]] or [6, 10, 20, 50, 100, 300, 1000]
KAP = np.logspace(0, 4, 121)
C0S = [0, 100, 300, 1000]

def betas(Lam):
    b1 = np.logspace(np.log10(0.3), np.log10(400.0), 140) / Lam      # b = Lam beta from 0.3 to 400
    b2 = np.linspace(0.02, 8.0, 160)
    bs = np.unique(np.concatenate([b1[b1 < 8], b2]))
    return bs

def run():
    res = {}
    for Dg in DGS:
        Lam = lam_of_gap(Dg)
        bs = betas(Lam)
        cs = [curvT_all(float(Lam), float(b)) for b in bs]
        for h in ['++', '+-']:
            s = SIGMA[h]
            out = {}
            worst = (np.inf, None, None)
            for C0 in C0S:
                kap_needed = 1.0
                for ib, b in enumerate(bs):
                    c = cs[ib]
                    m = KAP * Lam
                    u = 1 / KAP ** 2
                    B2 = 2 * c['psi2'] * np.ones_like(u)
                    B3 = s * u * c['psi3']
                    B4 = 10 * c['W'] * Lam ** 2 * u ** 2
                    B5 = 5 * s * u ** 3 * c['P0']
                    F = B2 + B3 + B4 + B5 - C0 * np.cosh(b) / m ** 2 * (np.abs(B2) + np.abs(B3) + np.abs(B4) + np.abs(B5))
                    phys = m * np.exp(-b / 2) >= Dg - 2
                    F = np.where(phys, F, np.inf)
                    neg = np.where(F <= 0)[0]
                    if neg.size:
                        kap_needed = max(kap_needed, KAP[min(neg.max() + 1, len(KAP) - 1)])
                    if C0 == 0:
                        k = np.argmin(F / (np.abs(B2) + np.abs(B3) + np.abs(B4) + np.abs(B5)))
                        r = F[k] / (abs(B2[k]) + abs(B3[k]) + abs(B4[k]) + abs(B5[k]))
                        if r < worst[0]:
                            worst = (r, KAP[k], b, F[k])
                out[C0] = kap_needed
            res[(Lam, h)] = (out, worst)
            print('Dg=%5d Lam=%8.2f %s  kappa needed (C0=%s): %s   min relative FR/|parts| = %.3e at m/Lam=%.2f beta=%.3f (FR=%.3e)'
                  % (Dg, Lam, h, '/'.join(map(str, C0S)), ' / '.join('%.3g' % out[C] for C in C0S), *worst), flush=True)
    return res

if __name__ == '__main__':
    run()
