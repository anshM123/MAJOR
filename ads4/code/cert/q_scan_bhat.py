"""Q3 flat ratio FF/N in the small-u, large-J region (float, [NE]): u on a log grid in [1e-6, 1e-2], all J up to
min(JCAP, (BMAX/sqrt(u) - 1)/2) by the recurrences of common_q.heavy_action (= rcommon.flat_integrals).  Usage: python3 q_scan_bhat.py [JCAP]"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
from common_q import heavy_action
d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
p2 = np.array([float(x) for x in d['psi2']]); p3 = np.array([float(x) for x in d['psi3']]); a = float(d['a']); c = float(d['c'])
JCAP = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
for h in ['++', '+-']:
    best = (9, None)
    for u in np.logspace(-6, -2, 161):
        Jm = int(min(JCAP, max(40, (200 / np.sqrt(u) - 1) / 2)))
        Js, FF, N = heavy_action(p2, p3, a, c, h, u, Jm)
        q = FF / N; j = int(np.argmin(q))
        if q[j] < best[0]: best = (q[j], int(Js[j]), u)
    r, J, u = best
    print('%s: min FF/N = %.6f at J = %d, u = %.4g, bhat = %.2f, m = %.2f M (JCAP %d)' % (h, r, J, u, (2 * J + 1) * np.sqrt(u), u ** -0.5, JCAP), flush=True)
