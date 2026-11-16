"""AdS4 uplift of a flat functional (ROBUST_FUNCTIONAL.md section 5): leading Regge action + envelopes, and the hybrid model.
Machinery of ADS4_FUNCTIONAL.md (t2, t3) with general weights.  Floating point; [NE]/[M], nothing certified.

Regge action (normal form, smearing scale Lam_s, u = (Lam_s/m)^2):
   FR = B2 + B3 + B4 + B5,  B2 = 2 int psi2 P_{-1/2+i Lam_s p}(cosh b),  B3 = sigma u int psi3 P,  B4 = 2a Lam_s^2 u^2 W(b),  B5 = a sigma u^3 P0(b).
Envelope:  FR - [C0 cosh(b)/m^2 + C1 u] (|B2|+|B3|+|B4|+|B5|) > 0.
Physical points: tau = 2 + m e^{-b/2} >= Dg (d = 3).  Lam = sqrt((Dg-2)(Dg-1)) = smallest heavy m; Lam_s = Lam/n.
Output per (Dg, n, h, envelope): the largest m/Lam (in units of the GAP scale Lam) at which the enveloped action is <= 0 on the grid
('last negative'); 1.000 or below means positive on every physical grid point.
Hybrid model [M]:  F_hyb = F_flat(J, u) + [Regge term - its flat limit at fixed bhat], term by term (ADS4_FUNCTIONAL.md section 2.4).
"""
import sys, json, os
import numpy as np
from rcommon import ReggeAction, SIGMA, C4, flat_parts, gl01, peval
from check import load

ENV = [(0, 0), (0, 1), (0, 3), (0, 10), (100, 1), (1000, 1), (1000, 3)]

def betas(Lam):
    b = np.concatenate([np.logspace(np.log10(0.01), np.log10(400.0), 300), np.arange(6.0, 12.0, 0.02), np.arange(12.0, 40.0, 0.05)])
    bs = np.concatenate([b / Lam, np.linspace(0.02, 8.0, 160)])
    return np.unique(bs[bs <= 8.0])

def regge_scan(F, Dg, n=1.0, verbose=True):
    Lam = np.sqrt((Dg - 2.0) * (Dg - 1.0))
    Ls = Lam / n
    RA = ReggeAction(F)
    KAP = np.unique(np.concatenate([np.logspace(0, 4, 1000), np.arange(1.0, 12.0, 0.004), np.arange(12, 40, 0.02)])) * n  # m/Lam_s
    m = KAP * Ls
    u = 1 / KAP ** 2
    last = {(h, e): 0.0 for h in SIGMA for e in ENV}
    worst = {h: (np.inf, None) for h in SIGMA}
    for b in betas(Ls):
        I2, I3, W, P0 = RA.transforms(Ls, b)
        phys = (m * np.exp(-b / 2) >= Dg - 2 - 1e-9) & (m >= Lam * (1 - 1e-12))
        for h, s in SIGMA.items():
            B2 = 2 * I2 * np.ones_like(u); B3 = s * u * I3; B4 = 2 * F.a * W * Ls ** 2 * u ** 2; B5 = F.a * s * u ** 3 * P0
            FR = B2 + B3 + B4 + B5
            N = np.abs(B2) + np.abs(B3) + np.abs(B4) + np.abs(B5)
            r = np.where(phys, FR / N, np.inf)
            k = np.argmin(r)
            if r[k] < worst[h][0]:
                worst[h] = (r[k], (round(KAP[k] / n, 3), round(b, 5)))
            for (C0, C1) in ENV:
                G = np.where(phys, FR - (C0 * np.cosh(b) / m ** 2 + C1 * u) * N, np.inf)
                neg = np.where(G <= 0)[0]
                if neg.size:
                    last[(h, (C0, C1))] = max(last[(h, (C0, C1))], KAP[neg.max()] / n)
    if verbose:
        for h in SIGMA:
            print('  Dg=%5d n=%.2f %s: min FR/N = %+.3e at (m/Lam, beta) = %s' % (Dg, n, h, *worst[h]))
            print('      last negative m/Lam: ' + '  '.join('%s:%.3f' % (e, last[(h, e)]) for e in ENV), flush=True)
    return last, worst

if __name__ == '__main__':
    name = sys.argv[1]
    F = load(name)
    for Dg in [int(x) for x in sys.argv[2:]] or [20, 100, 1000]:
        regge_scan(F, Dg)
