"""Lower bound on g4 in terms of |g3|^2 (tangent at the origin; CHLPS (4.5): 0.26; Maldacena et al. (6.13): 1/4).  [NE] for the LP,
[P given the forward sum rules] for the psi3 = 0 value 1/4.
Functional F = int psi3 B3 - a (dB4|_0 - B5|_0) + (1/2) B4|_0  (psi2 = 0 is forced: a positive functional with psi2 != 0 has
psi2'(0) > 0 and hence a log term, which cannot produce a tangent through the origin).  Heavy (m^4 x):
   sigma u int psi3 D_h + a (2JJ - c4 + sigma) u^3 + u^2 >= 0,   low:  g4/(8 pi G) >= s |g3|^2 M^2,  s = (2 mu3_4 + 4a)/8.
Usage: python3 slope.py d3 r [delta]"""
import sys, json
import numpy as np
from scipy.optimize import linprog
from common_q import Grid, SIGMA
d3, r = int(sys.argv[1]), int(sys.argv[2]); delta = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
G = Grid(2, d3, r)
n3 = G.n3; nv = n3 + 1
obj = np.zeros(nv); obj[:n3] = -2 * G.mu3_4 / 8; obj[-1] = -4 / 8
def rows(idx, s3):
    R = np.zeros((len(idx), nv)); R[:, :n3] = (1 - delta * s3)[:, None] * G.B3[idx]
    R[:, -1] = (1 - delta * np.sign(G.FWc[idx])) * G.FWc[idx]
    return -R, (1 - delta) * G.CB4[idx] / 2         # -(B3 + a FW) <= (1-delta) u^2
fixed, bf = [], []
for h in SIGMA:                                      # large bhat asymptotics
    A2, A3, Fc, Cb = G.asym[h]
    R = np.zeros((len(G.ab), nv)); R[:, :n3] = -A3; R[:, -1] = -Fc; fixed.append(R); bf.append(Cb / 2)
init = np.arange(0, len(G.u), max(1, len(G.u) // 3000))
cuts, cb = [], []
for s3 in (1, -1):
    R, b = rows(init, np.full(len(init), s3)); cuts.append(R); cb.append(b)
for it in range(60):
    A = np.concatenate(fixed + cuts); b = np.concatenate(bf + cb)
    sc = np.maximum(np.abs(A).max(axis=1), 1e-300)
    res = linprog(obj, A_ub=A / sc[:, None], b_ub=b / sc, bounds=[(-1e4, 1e4)] * n3 + [(0, 1e4)], method='highs')
    if res.status: print('status', res.status, res.message); sys.exit()
    v = res.x; B3 = G.B3 @ v[:n3]; FW = v[-1] * G.FWc; FF = B3 + FW + G.CB4 / 2; N = abs(B3) + abs(FW) + G.CB4 / 2
    g = FF - delta * N; viol = np.where(g < -1e-8 * N)[0]
    if len(viol) == 0: break
    o = viol[np.argsort(g[viol] / N[viol])][:3000]; R, bb = rows(o, np.sign(B3[o])); cuts.append(R); cb.append(bb)
print(json.dumps({'d3': d3, 'r': r, 'delta': delta, 'slope': -res.fun, 'a': v[-1], 'psi3_0': float(G.b3[0][0] * v[0]), 'it': it,
                  'minFF_N': float(np.min(FF / N))}))
