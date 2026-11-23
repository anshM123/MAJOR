"""
lsw/lsw_costs.py (LARGESPIN_W.md §4): exact per-unit costs (Lemma I) of every structure that the weight-shifted CHMRS expansion can
produce at O(1/m^2) and O(1/m^4), in de-shift (R1) and stripped (R2) form, relative to the per-helicity floor, on the H_LJ(1000)
sample states (same list as largespin/ls_struct.py) and on an extended near-gap list.

Output: cache/costs_1000.json, cache/costs_1000.txt
  per state: J, m/Lam, beta, bhat, floors per helicity, and for each structure S: c1 = |R1 cost|*m^2/floor, c2 = |R2 cost|*m^2/floor
  (i.e. the fraction of the floor consumed by a coefficient kappa = 1 in c2 units: dI = S/m^2).
Usage: python3 lsw_costs.py
"""
import json, os, sys, time
sys.dont_write_bytecode = True
import numpy as np
import lsw_common as W

t0 = time.time()
Js = [26, 30, 34, 38, 42, 46, 50, 64, 80, 100, 120, 150, 200, 250, 300, 344, 33, 51, 101, 201, 301]
rows = []
lines = []


def pr(s):
    print(s, flush=True); lines.append(s)


for J in Js:
    mg = W.C.mgap(W.Dg, J)
    for f in (1.0, 1.1, 1.2):
        mfl = mg*f
        if mfl > 1.3:
            continue
        m = mfl*W.L
        fl = W.floors(J, mfl)
        beta = list(fl.values())[0]['beta']
        S = W.structures(beta)
        cost = {k: W.act(v) for k, v in S.items()}
        r = dict(J=J, mfl=mfl, m=m, beta=beta, bhat=(2*J + 1)/mfl, fl={h: v for h, v in fl.items()},
                 cost={k: [c[0]/m**2, c[1]/m**2] for k, c in cost.items()})
        rows.append(r)
        fmin = min(v['floor'] for v in fl.values())
        pr('J=%3d m/L=%.3f beta=%.4f bhat=%6.1f floor(min h)=%.3e | ' % (J, mfl, beta, r['bhat'], fmin) +
           ' '.join('%s:%.1e/%.1e' % (k, abs(c[0])/m**2/fmin, abs(c[1])/m**2/fmin) for k, c in cost.items()) +
           '  (%.0fs)' % (time.time() - t0))

json.dump(rows, open(os.path.join(W.CACHE, 'costs_1000.json'), 'w'))
open(os.path.join(W.CACHE, 'costs_1000.txt'), 'w').write('\n'.join(lines) + '\n')
