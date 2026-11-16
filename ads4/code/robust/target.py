"""Envelope-targeted LP (ROBUST_FUNCTIONAL.md section 2.3): for envelope strength C1 and rescaling n (smearing scale Lam_s = Lam/n,
so physical states have u = (Lam_s/m)^2 <= 1/n^2), impose  FF - (delta0 + C1 u) N >= 0 for u <= 1/n^2 (R1), R2 with dsc = delta0, R3, matter;
minimise the exact light-side value 16 L_AdS(Lam/n) = 16 int psi2/(p^2 + 9 n^2/(4 Lam^2)) (k = 1, LAM = 100).  The would-be AdS bound is |ghat3|^2 Lam^8 <= n^8 16 L_AdS(Lam/n) (+ O(Lam^-2) light-side terms).  [NE]"""
import sys, json, os
import numpy as np
from lp import Grid, solve, CACHE

d2, d3 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (4, 4)
delta0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.02
LAM = 100.0
G = Grid(d2, d3, verbose=False)
out = []
for C1 in [1.0, 3.0, 10.0]:
    for n in [1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0]:
        st, v, info = solve(G, delta0, C=C1, u0=1 / n ** 2, LAMS=LAM / n)
        if st != 0:
            print('C1=%g n=%.2f: status %d %s' % (C1, n, st, info.get('msg')), flush=True)
            continue
        from lp import LAdS_rows
        val = float(LAdS_rows(G, LAM / n) @ v[:G.n2])
        print('C1=%g n=%.2f: A=%.3f B=%.3f B_AdS=%.3f a=%.4f  16 L_AdS(Lam/n)=%.3f  n^8 x that = %.4g' % (
            C1, n, info['A'], info['B'], info['B_AdS'], info['a'], val, n ** 8 * val), flush=True)
        out.append(dict(C1=C1, n=n, delta0=delta0, d2=d2, d3=d3, A=info['A'], B=info['B'], B_AdS=info['B_AdS'], a=info['a'],
                        bound=n ** 8 * val, v=list(map(float, v))))
        break_ok = True
json.dump(out, open(os.path.join(CACHE, 'target_d%d_%d_%g.json' % (d2, d3, delta0)), 'w'))
