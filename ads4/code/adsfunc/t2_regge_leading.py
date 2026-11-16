"""t2: the leading Regge-limit AdS action of Omega_g3 (heavy formula C^pm_{k,nu} = (2P/m^{2k})(|c_s|^2 +- |c_o|^2), ASSUMED at
k = 2, 3, 4, 5 and for odd k), normal form  FR = B2 + B3 + B4 + B5  (common.regge_parts), on (m/Lam, beta) grids, both channels.

Only PHYSICAL points: twist tau = 2 + alpha >= Dg, alpha = m e^{-beta/2} (d = 3), i.e. beta <= 2 log(m/(Dg-2));
Lam = sqrt((Dg-2)(Dg-1)); on physical points cosh(beta)/m^2 <= 1/(Dg-2)^2.

Envelope (a MODEL of the unknown corrections; revised after REFEREE_ADS4_FUNCTIONAL.md issue 1):
   FR - [C0 cosh(beta)/m^2 + C1 Lam^2/m^2] (|B2|+|B3|+|B4|+|B5|),
C0: the 1/m^2 corrections of the heavy formula; C1: the O(nu^2/m^2) corrections (nu <= Lam), C1 = 1 ~ the size of the flat O(u) terms.
Reported per (Dg, h, C0, C1): the LARGEST m/Lam on the grid at which the (enveloped) action is <= 0 ("last negative"; positivity holds on the
grid above it).  Grids (revised after issue 3): Lam*beta in [0.01, 400] (300 log points) plus step 0.01 on [6, 12] (++ island) and step 0.02 on
[12, 30] (+- near-zero region), plus 160 linear beta points in [0.02, 8]; m/Lam: 1200 log points in [1, 1e4] plus step 0.002 on [5, 12].
Floating point; not a certificate.
"""
import sys
import numpy as np
from common import curvT_all, SIGMA, lam_of_gap

DGS = [int(a) for a in sys.argv[1:]] or [6, 10, 20, 50, 100, 300, 1000]
KAP = np.unique(np.concatenate([np.logspace(0, 4, 1200), np.arange(5.0, 12.0, 0.002)]))
ENV = [(0, 0), (100, 0), (300, 0), (1000, 0), (0, 1), (100, 1), (300, 1), (1000, 1), (300, 3)]

def betas(Lam):
    b = np.concatenate([np.logspace(np.log10(0.01), np.log10(400.0), 300), np.arange(6.0, 12.0, 0.01), np.arange(12.0, 30.0, 0.02)])
    bs = np.concatenate([b / Lam, np.linspace(0.02, 8.0, 160)])
    return np.unique(bs[bs <= 8.0])

def run():
    for Dg in DGS:
        Lam = lam_of_gap(Dg)
        bs = betas(Lam)
        m = KAP * Lam
        u = 1 / KAP ** 2
        last = {(h, e): 0.0 for h in SIGMA for e in ENV}
        worst = {h: (np.inf, None) for h in SIGMA}
        for b in bs:
            c = curvT_all(float(Lam), float(b))
            phys = m * np.exp(-b / 2) >= Dg - 2
            for h, s in SIGMA.items():
                B2 = 2 * c['psi2'] * np.ones_like(u)
                B3 = s * u * c['psi3']
                B4 = 10 * c['W'] * Lam ** 2 * u ** 2
                B5 = 5 * s * u ** 3 * c['P0']
                F = B2 + B3 + B4 + B5
                N = np.abs(B2) + np.abs(B3) + np.abs(B4) + np.abs(B5)
                r = np.where(phys, F / N, np.inf)
                k = np.argmin(r)
                if r[k] < worst[h][0]:
                    worst[h] = (r[k], (round(KAP[k], 3), round(b, 5), round(Lam * b, 3)))
                for (C0, C1) in ENV:
                    G = np.where(phys, F - (C0 * np.cosh(b) / m ** 2 + C1 * Lam ** 2 / m ** 2) * N, np.inf)
                    neg = np.where(G <= 0)[0]
                    if neg.size:
                        last[(h, (C0, C1))] = max(last[(h, (C0, C1))], KAP[neg.max()])
        for h in SIGMA:
            print('Dg=%5d Lam=%8.2f %s  min FR/N = %+.3e at (m/Lam, beta, Lam*beta) = %s' % (Dg, Lam, h, *worst[h]))
            print('      last negative m/Lam for (C0,C1) = ' + '  '.join('%s:%.3f' % (e, last[(h, e)]) for e in ENV), flush=True)

if __name__ == '__main__':
    run()
