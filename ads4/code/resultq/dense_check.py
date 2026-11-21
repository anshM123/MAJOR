"""Dense off-grid check of LP functionals stored in cache/lp_g4.jsonl (line index argv[1], default: all with status 0).  [NE]
u grid disjoint from the LP grid (shifted log grid + fine near m = M), J up to max(3000, bhat 400).  Also the scaling limit
2 int psi2 J0(bhat p) for bhat in [0, 2000]."""
import sys, json, os
import numpy as np
from scipy.special import j0
from common_q import heavy_action, CACHE
from rcommon import jmax_for, gl01, peval
L = [json.loads(l) for l in open(os.path.join(CACHE, 'lp_g4.jsonl'))]
idx = [int(a) for a in sys.argv[1:]] or [i for i, r in enumerate(L) if r.get('status') == 0]
us = np.unique(np.concatenate([np.logspace(-5.3, 0, 97)[:-1] * 1.013, np.linspace(0.97, 1, 31), np.linspace(0.3, 0.97, 35)]))
us = us[us <= 1]
p, w = gl01(4000); bh = np.linspace(0, 2000, 40001)
for i in idx:
    r = L[i]
    psi2, psi3, a, c = np.array(r['psi2']), np.array(r['psi3']), r['a'], r['c']
    mn = (np.inf, None); neg = 0
    for u in us:
        for h in ['++', '+-']:
            Js, FF, N = heavy_action(psi2, psi3, a, c, h, u, jmax_for(u, 400.0, 3000))
            q = FF / N; j = np.argmin(q); neg += int((q < -1e-6).sum())
            if q[j] < mn[0]:
                mn = (q[j], (round(1 / np.sqrt(u), 4), int(Js[j]), h))
    S = 2 * (j0(np.outer(bh, p)) @ (w * peval(psi2, p)))
    print('line %d d=%d/%d r=%d delta=%.2f matter=%s b3=%s A=%.4f B=%.4f | dense min FF/N %.3e at %s, #neg(<-1e-6) %d | scaling min %.3e at bhat %.1f'
          % (i, r['d2'], r['d3'], r['r'], r['delta'], r['matter'], r['b3'], r['A'], r['B'], mn[0], mn[1], neg, S.min(), bh[np.argmin(S)]), flush=True)
