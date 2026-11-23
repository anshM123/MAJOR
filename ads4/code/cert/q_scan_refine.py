"""Refine the Q3 ++ large-J minimum (float, [NE]): u on a fine grid around the q_scan_bhat minimum, all J by recurrence."""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
from common_q import heavy_action
d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
p2 = np.array([float(x) for x in d['psi2']]); p3 = np.array([float(x) for x in d['psi3']]); a = float(d['a']); c = float(d['c'])
best = (9, None)
for u in np.geomspace(4e-5, 3e-4, 81):
    Jm = int(min(30000, (200 / np.sqrt(u) - 1) / 2))
    Js, FF, N = heavy_action(p2, p3, a, c, '++', u, Jm)
    q = FF / N; j = int(np.argmin(q))
    if q[j] < best[0]: best = (q[j], int(Js[j]), u)
r, J, u = best
print('++ refined: min FF/N = %.6f at J = %d, u = %.5g, bhat = %.3f, m = %.2f M' % (r, J, u, (2 * J + 1) * np.sqrt(u), u ** -0.5), flush=True)
