"""Pareto scan of the certifiable family (ROBUST_FUNCTIONAL.md section 2.2): criterion R1 (C = 0, u0 = 1) with margin delta, R2 (dsc = delta),
R3, matter, plus the Lemma-M proxy (BA = 120, u1 = 1/100, |psi3(0)| <= Y = 2, a >= 0.14 per unit k).  Minimise A (tie-break B).  [NE]"""
import sys, json, os
import numpy as np
from lp import Grid, solve, CACHE
d2, d3 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (4, 4)
G = Grid(d2, d3, verbose=False)
out = []
for delta in [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5]:
    st, v, info = solve(G, delta, mreg=dict(BA=120.0, u1=0.01, Y=2.0, amin=0.14, safety=1.4))
    if st != 0:
        print('d=(%d,%d) delta=%.2f: status %d %s' % (d2, d3, delta, st, info.get('msg')), flush=True); continue
    print('d=(%d,%d) delta=%.2f: A=%.4f B=%.4f B_AdS=%.4f a=%.4f  A log100 + B_AdS = %.3f' % (
        d2, d3, delta, info['A'], info['B'], info['B_AdS'], info['a'], info['A'] * np.log(100) + info['B_AdS']), flush=True)
    out.append(dict(d2=d2, d3=d3, delta=delta, A=info['A'], B=info['B'], B_AdS=info['B_AdS'], a=info['a'], v=list(map(float, v))))
    json.dump(out, open(os.path.join(CACHE, 'pareto2_d%d_%d.json' % (d2, d3)), 'w'))
