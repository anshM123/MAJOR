"""
tailaudit/ta_contour_rows.py [DV, double precision]: spot check of the CONTOUR_BOUNDS certified rows that TAIL_CERT uses as inputs:
per-cell int|fB| on the short detour, the vertical piece and the leg Y = 9/4 (cb_short, cb_vert, cb_horiz_2.25), recomputed by
Gauss-Legendre with the independent numpy f of ta_contour_bookkeeping.py, and the totals.  A row passes if computed <= certified
(within 1e-6 relative for double-precision noise; the RBF cancellation at Y = 9/4 is ~e^{26}, so ~1e-5 relative noise is expected).
Output: cache/ta_contour_rows.txt
"""
import os, sys, json, math, random
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import ta_contour_bookkeeping as TB      # only the f, B evaluators are used (module body re-runs its own check; cheap)
ADS4 = TB.ADS4
CB = os.path.join(ADS4, 'code', 'contour', 'cache')
out = []
def pr(s):
    print(s, flush=True); out.append(s)
GLX, GLW = np.polynomial.legendre.leggauss(30)
random.seed(7)
for fn, kind in (('cb_short.json', 'short'), ('cb_vert.json', 'vert'), ('cb_horiz_2.25.json', 'leg')):
    rows = json.load(open(os.path.join(CB, fn)))['rows']
    idx = list(range(len(rows))) if len(rows) < 200 else sorted(random.sample(range(len(rows)), 150) + [0, len(rows) - 1] +
                                                                  [i for i in range(len(rows)) if rows[i].get('t0', 0) < 45][::400])
    worst = 0.0; tot_c = 0.0; tot_r = 0.0
    for i in idx:
        r = rows[i]
        if kind == 'vert':
            y0, y1 = r['y0'], r['y1']; s = (y1 - y0)/2*(GLX + 1) + y0; nu = 0.25 + 1j*s; wts = (y1 - y0)/2*GLW
        else:
            Y = 0.75 if kind == 'short' else 2.25
            t0, t1 = r['t0'], r['t1']; s = (t1 - t0)/2*(GLX + 1) + t0; nu = s + 1j*Y; wts = (t1 - t0)/2*GLW
        v = float(np.sum(wts*np.abs(TB.f(nu)*TB.B(nu))))
        worst = max(worst, v/r['intfB'])
        tot_c += v; tot_r += r['intfB']
    pr(f'{kind}: {len(idx)} rows sampled: max computed/certified int|fB| = {worst:.6f};  sampled totals computed {tot_c:.6e} vs certified '
       f'{tot_r:.6e}')
open(os.path.join(HERE, 'cache', 'ta_contour_rows.txt'), 'w').write('\n'.join(out) + '\n')
