"""Q3 float minima of FF/N in the two regions the Q3* certificate does not cover (CERT_CLOSE.md §2.3) [NE]:
  (a) J >= 400, u <= 1/100, bhat >= 120 (both channels);  (b) +-, J >= 400, u >= 1/4  (Q3*: factored weights, q_backward split)."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'resultq')); sys.path.insert(0, HERE)
from common_q import heavy_action
import q_backward as QB
d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
p2 = np.array([float(x) for x in d['psi2']]); p3 = np.array([float(x) for x in d['psi3']]); a = float(d['a']); c = float(d['c'])
for h in ['++', '+-']:
    best = (9, None)
    for u in np.logspace(-6, -2, 121):
        Jm = int(min(20000, (600 / np.sqrt(u) - 1) / 2))
        Js, FF, N = heavy_action(p2, p3, a, c, h, u, Jm)
        b = (2 * Js + 1) * np.sqrt(u); k = (Js >= 400) & (b >= 120)
        if not np.any(k): continue
        q = FF[k] / N[k]; j = int(np.argmin(q))
        if q[j] < best[0]: best = (q[j], int(Js[k][j]), u, b[k][j])
    print('(a) %s, J >= 400, u <= 1/100, bhat >= 120 (bhat <= 600, J <= 2e4): min FF/N = %.5f at J = %d, u = %.3g, bhat = %.1f' % ((h,) + best), flush=True)
best = (9, None)
for u in np.concatenate([np.linspace(0.25, 0.99, 38), [0.995, 0.999, 0.9999, 1.0]]):
    Js, FF, N = QB.action('+-', u, 6000, with_rem=False)
    k = Js >= 400; q = FF[k] / N[k]; j = int(np.argmin(q))
    if q[j] < best[0]: best = (q[j], int(Js[k][j]), u)
print('(b) +-, J >= 400, u >= 1/4 (Q3*, J <= 6000): min FF/N = %.5f at J = %d, u = %.4g' % best, flush=True)
