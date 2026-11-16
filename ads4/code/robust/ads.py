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


# ------------------------------------------------------------------ hybrid model [M]
def flatT(c, b):
    from rcommon import glang, nth
    th, wt = glang(nth(b))
    return (2 / np.pi) * np.sum(wt * c(b * np.sin(th)))

def jlist(h, mfl, jcap=1500):
    jmax = int(min(jcap, 150 * mfl))
    js = list(range(0, 41)); j = 41
    while j <= jmax:
        js.append(j); j = int(j * 1.04) + 1
    js = sorted(set(x for x in js if x <= jmax))
    return np.array(sorted(set(x + (x % 2) for x in js)) if h == '++' else [x for x in js if x >= 4])

def hybrid_scan(F, Dg, n=1.0, verbose=True):
    """criterion (A) of ADS4_FUNCTIONAL.md t3 with envelope (C0 cosh b/m^2 + C1 u) N_comp; returns last negative m/Lam per envelope."""
    Lam = np.sqrt((Dg - 2.0) * (Dg - 1.0)); Ls = Lam / n
    RA = ReggeAction(F)
    MFL = np.unique(np.concatenate([np.arange(1.0, 1.5, 0.02), np.logspace(np.log10(1.5), np.log10(60), 90)])) * n
    last = {(h, e): 0.0 for h in SIGMA for e in ENV}; worst = {h: (np.inf, None) for h in SIGMA}; nst = 0
    for mfl in MFL:
        u = 1 / mfl ** 2; m = mfl * Ls
        for h, s in SIGMA.items():
            js = jlist(h, mfl)
            Js, B2, B3, FW = flat_parts(F, h, u, int(js.max()))
            idx = {int(J): i for i, J in enumerate(Js)}
            for J in js:
                # physical: tau >= Dg  <=>  alpha >= Dg-2 with alpha*(alpha+2J+1) = m^2
                Bq = 2 * J + 1.0
                alpha = (-Bq + np.sqrt(Bq * Bq + 4 * m * m)) / 2
                if alpha < Dg - 2 - 1e-9:
                    continue
                abar = alpha + 2 * J + 1
                beta = np.arccosh((alpha ** 2 + abar ** 2) / (2 * alpha * abar))
                i = idx[int(J)]
                I2, I3, W, P0 = RA.transforms(Ls, beta)
                bh = (2 * J + 1) / mfl
                f2 = flatT(RA.c2, bh); f3 = flatT(RA.c3, bh)
                dB2 = 2 * (I2 - f2); dB3 = s * u * (I3 - f3)
                dFW = u ** 3 * (2 * F.a * m * m * W - 2 * F.a * J * (J + 1) - F.a / 2) + F.a * s * u ** 3 * (P0 - 1)
                Fh = B2[i] + B3[i] + FW[i] + dB2 + dB3 + dFW
                Nc = abs(B2[i] + dB2) + abs(B3[i] + dB3) + abs(FW[i] + dFW)
                nst += 1
                if Fh / Nc < worst[h][0]:
                    worst[h] = (Fh / Nc, (round(mfl / n, 3), int(J), round(beta, 4)))
                for (C0, C1) in ENV:
                    if Fh - (C0 * np.cosh(beta) / m ** 2 + C1 * u) * Nc <= 0:
                        last[(h, (C0, C1))] = max(last[(h, (C0, C1))], mfl / n)
    if verbose:
        for h in SIGMA:
            print('  hybrid Dg=%5d n=%.2f %s (%d states): min F_hyb/N = %+.3e at (m/Lam, J, beta) = %s' % (Dg, n, h, nst, *worst[h]))
            print('      last negative m/Lam: ' + '  '.join('%s:%.3f' % (e, last[(h, e)]) for e in ENV), flush=True)
    return last, worst

if __name__ == '__main__':
    name, mode = sys.argv[1], sys.argv[2]
    F = load(name)
    for Dg in [int(x) for x in sys.argv[3:]] or [20, 100, 1000]:
        print('== %s %s Dg=%d' % (name, mode, Dg), flush=True)
        (regge_scan if mode == 'regge' else hybrid_scan)(F, Dg)
