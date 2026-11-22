"""
rd_summary (RD_CLOSE.md sections 3-5): tables from the rd_flow caches.
  * scalar far-region exponent on the crystal thimbles, c(r) = min s/Delta over thimble points with raw radius >= r (min over A+-/S);
  * EE far-region exponent cEE(r) = min [s - log(|Q|/Q_eff)]/Delta;
  * fit  Delta c(r) = c_inf(r) Delta - p log Delta + b  over Delta = 100, 300, 1000 (and 30 where the saddle is inside r);
  * the far-region EE log-weight at Delta_g = 300, 1000 against the failure thresholds log(2e-4), log(2e-5).
Usage: python3 rd_summary.py -> cache/rd_summary.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from rd_flow import RADII
np.seterr(all='ignore')

DS = (30, 100, 300, 1000); NUS = (2, 4, 6, 8, 12, 16)
tab = {}
for D in DS:
    for nu in NUS:
        f = os.path.join(HERE, 'cache', f'rd_flow_D{D}_nu{nu}.pkl')
        if not os.path.exists(f):
            continue
        R = pickle.load(open(f, 'rb'))
        cs = []; ce = []; qm = []
        for k, o in R.items():
            lQe = np.log(o['qnear'])
            sf = np.where(np.isfinite(o['sfirst']), o['sfirst'], np.inf)
            cs.append(np.min(sf, axis=0)/D)
            ce.append(np.min(np.where(np.isfinite(o['sfirst']), o['cee'] + lQe, np.inf), axis=0)/D)
            qm.append(np.max(np.where(np.isfinite(o['sfirst']), o['qmax'] - lQe, -np.inf), axis=0))
        tab[(D, nu)] = (np.min(cs, axis=0), np.min(ce, axis=0), np.max(qm, axis=0))

L = ['# rd_summary: far-region exponents on the crystal thimbles (A+- below the fold, S above; n3/n4 flow function, n0 = +1)',
     '#   c(r): scalar, min action drop/Delta at raw radius >= r;  cEE(r): EE (scalar x |Q|/Q_eff);  analytic leading value artanh(r)']
ri = [RADII.index(r) for r in (0.1, 0.2, 0.3, 0.5, 0.7)]
L.append('   Delta  nu  | c(0.1)  c(0.2)  c(0.3)  c(0.5)  c(0.7) | cEE(0.1) cEE(0.2) cEE(0.3) cEE(0.5) cEE(0.7)')
for D in DS:
    for nu in NUS:
        if (D, nu) not in tab:
            continue
        c, e, q = tab[(D, nu)]
        L.append(f'  {D:5d}  {nu:3d}  | ' + ' '.join(f'{c[i]:6.3f} ' for i in ri) + '| ' + ' '.join(f'{e[i]:7.3f} ' for i in ri))
L.append('   artanh(r):     ' + ' '.join(f'{np.arctanh(r):6.3f} ' for r in (0.1, 0.2, 0.3, 0.5, 0.7)))
L.append('')
L.append('# fit of the action drop: Delta c(r) = c_inf Delta - p log Delta + b (min over nu < 10, i.e. the worst nu), Delta = 100, 300, 1000')
for r in (0.1, 0.2, 0.3, 0.5, 0.7):
    k = RADII.index(r)
    for kind, idx in (('scal', 0), ('EE', 1)):
        y = []; Ds = []
        for D in (100, 300, 1000):
            vals = [tab[(D, nu)][idx][k]*D for nu in (2, 4, 6, 8) if (D, nu) in tab]
            if vals:
                y.append(min(vals)); Ds.append(D)
        Ds = np.array(Ds, float); y = np.array(y)
        A = np.stack([Ds, -np.log(Ds), np.ones_like(Ds)], axis=1)
        coef = np.linalg.solve(A, y) if len(Ds) == 3 else np.linalg.lstsq(A, y, rcond=None)[0]
        A2 = np.stack([Ds, np.ones_like(Ds)], axis=1)
        c2 = np.linalg.lstsq(A2, y, rcond=None)[0]
        L.append(f'   r = {r}: {kind:4s} drops {np.round(y, 1)} -> c_inf = {coef[0]:.4f} (p = {coef[1]:+.2f}, b = {coef[2]:+.2f}); '
                 f'linear fit c = {c2[0]:.4f}; artanh(r) = {np.arctanh(r):.4f}')
L.append('')
L.append('# far-region EE log-weight (max over the thimble points with raw radius >= r, relative to the near region; worst nu < 10)')
L.append('#   failure threshold: log(2e-4) = %.2f (Delta_g = 300), log(2e-5) = %.2f (1000)' % (np.log(2e-4), np.log(2e-5)))
for D in (300, 1000):
    row = []
    for r in (0.1, 0.2, 0.3, 0.5):
        k = RADII.index(r)
        w = max(-tab[(D, nu)][1][k]*D for nu in (2, 4, 6, 8))
        row.append(f'r >= {r}: {w:8.1f}')
    L.append(f'   Delta = {D}: ' + ';  '.join(row))
L.append('')
L.append('# max log(|Q|/Q_eff) met beyond r (polynomial growth of EE/scalar on the far thimble, worst nu < 10)')
for D in DS:
    row = []
    for r in (0.1, 0.2, 0.3, 0.5):
        k = RADII.index(r)
        row.append(f'{max(tab[(D, nu)][2][k] for nu in (2, 4, 6, 8)):+7.1f}')
    L.append(f'   Delta = {D:5d}: r >= 0.1 / 0.2 / 0.3 / 0.5: ' + ' / '.join(row))
open(os.path.join(HERE, 'cache', 'rd_summary.txt'), 'w').write('\n'.join(L) + '\n')
print('\n'.join(L))
