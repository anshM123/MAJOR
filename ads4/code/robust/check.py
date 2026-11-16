"""Dense floating-point check of a flat functional (ROBUST_FUNCTIONAL.md section 3).  [NE]

For each channel: min of FF/N over physical states (u = 1/m^2 in [1e-6, 1], J <= min(JCAP, (BMAX/sqrt(u) - 1)/2));
the envelope thresholds n(C) := sup{ m : FF - C u N < 0 }  (C = C1 of the AdS envelope C1 Lam^2/m^2, flat parts);
scaling limit 2 f^(bhat) / rho(bhat) on [0, 3000]; Lemma-M asymptotic data.
Usage: python3 check.py NAME   (NAME in cache/functionals.json, or 'Fg3')
"""
import sys, json, os
import numpy as np
from rcommon import Functional, paper_fg3, flat_parts, jmax_for, peval, gl01, SIGMA
from lp import rho, CACHE

CS = [0.1, 0.3, 1.0, 3.0, 10.0]

def load(name):
    if name == 'Fg3':
        return paper_fg3().normalised()
    d = json.load(open(os.path.join(CACHE, 'functionals.json')))[name]
    return Functional(d['psi2'], d['psi3'], d['a'], name)

def dense(F, BMAX=200.0, JCAP=3000, nu=140, verbose=True):
    us = np.unique(np.concatenate([np.logspace(-6, 0, nu), np.linspace(0.9, 1, 21)]))
    out = {}
    for h in ['++', '+-']:
        worst = (np.inf, None); thr = {C: 0.0 for C in CS}; worstpos = (np.inf, None)
        for u in us:
            Jm = jmax_for(u, BMAX, JCAP)
            Js, B2, B3, FW = flat_parts(F, h, u, Jm)
            if h == '++':
                s = Js % 2 == 0; Js, B2, B3, FW = Js[s], B2[s], B3[s], FW[s]
            FF = B2 + B3 + FW; N = np.abs(B2) + np.abs(B3) + np.abs(FW)
            r = FF / N
            i = np.argmin(r)
            if r[i] < worst[0]:
                worst = (r[i], dict(m=1 / np.sqrt(u), J=int(Js[i]), bhat=(2 * Js[i] + 1) * np.sqrt(u), FF=FF[i], N=N[i]))
            i = np.argmin(FF)
            if FF[i] < worstpos[0]:
                worstpos = (FF[i], dict(m=1 / np.sqrt(u), J=int(Js[i])))
            for C in CS:
                if np.any(FF - C * u * N < 0):
                    thr[C] = max(thr[C], 1 / np.sqrt(u))
        out[h] = dict(min_ratio=worst, min_FF=worstpos, thr=thr)
        if verbose:
            print('  %s: min FF/N = %.4e at %s' % (h, worst[0], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in worst[1].items()}))
            print('      min FF = %.4e at %s' % worstpos)
            print('      envelope thresholds n(C) = sup{m : FF < C u N}: ' + '  '.join('C=%g: %.3f' % (C, thr[C]) for C in CS), flush=True)
    return out

def scaling(F, verbose=True):
    from scipy.special import j0
    p, w = gl01(6000)
    bh = np.concatenate([np.linspace(0, 100, 5001), np.linspace(100, 3000, 6000)])
    f = np.array([2 * np.sum(w * peval(F.psi2, p) * j0(b * p)) for b in bh])
    r = f / rho(bh)
    i = np.argmin(r)
    c3 = -2 * F.psi2[2]          # coefficient of 1/bhat^3 in 2 f^
    if verbose:
        print('  scaling limit: min 2f^/rho = %.4f at bhat = %.3f ; 2f^ ~ %.5f/bhat^3 (rho tail %.5f/bhat^3); psi3(0) = %.5f, a = %.5f'
              % (r[i], bh[i], c3, 80 / 27.48509823509836, F.psi3[0], F.a))
    return dict(min_ratio=float(r[i]), at=float(bh[i]), c3=float(c3))

if __name__ == '__main__':
    for name in sys.argv[1:] or ['Fg3']:
        F = load(name)
        A, B, BA = F.AB()
        print('== %s: k=%.6f A=%.4f B=%.4f B_AdS=%.4f a=%.5f' % (name, F.k(), A, B, BA, F.a), flush=True)
        scaling(F)
        dense(F)
