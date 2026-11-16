"""Pareto scan: for each (d2, d3) and margin delta (criterion R1 with C = 0, u0 = 1, R2 with dsc = delta), the minimal A
(tie-break B) and the minimal A log(100) + B_AdS.  Output cache/pareto.json.  Floating point LP: [NE]."""
import sys, json, os
import numpy as np
from lp import Grid, solve, CACHE

res = []
degs = [(4, 4), (6, 6), (8, 8)] if len(sys.argv) < 2 else [tuple(int(t) for t in a.split(',')) for a in sys.argv[1:]]
for d2, d3 in degs:
    G = Grid(d2, d3, verbose=False)
    for delta in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6]:
        for LREF in [np.inf]:
            st, v, info = solve(G, delta, LREF=LREF)
            if st == 0:
                row = dict(d2=d2, d3=d3, delta=delta, LREF=str(LREF), A=info['A'], B=info['B'], B_AdS=info['B_AdS'], a=info['a'],
                           v=list(map(float, v)))
                print('d=(%d,%d) delta=%.2f LREF=%s: A=%.4f B=%.4f B_AdS=%.4f a=%.4f  [A log100 + B_AdS = %.3f]' % (
                    d2, d3, delta, LREF, info['A'], info['B'], info['B_AdS'], info['a'], info['A'] * np.log(100) + info['B_AdS']), flush=True)
                res.append(row)
            else:
                print('d=(%d,%d) delta=%.2f LREF=%s: status %d %s' % (d2, d3, delta, LREF, st, info.get('msg')), flush=True)
                if st == 2:
                    break
json.dump(res, open(os.path.join(CACHE, 'pareto_%s.json' % '_'.join('%d-%d' % d for d in degs)), 'w'))
