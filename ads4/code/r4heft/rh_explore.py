"""quick exploration: which s-channel residues d12 = -n survive per class (k = 2, 4), P1 angles."""
import sys, time
import numpy as np
from scipy.special import loggamma
import rh_common as R

k = int(sys.argv[1])
a, b, p2, p4 = R.P.config('P1')
t0 = time.time()
S = R.setup(k, a, b, p2, p4)
print('setup', len(S.t), f'{time.time()-t0:.0f}s', flush=True)
d13 = 2/3 + 0.7j
for name, f in [('D4', R.poly_D(4)), ('D5', R.poly_D(5)), ('D6', R.poly_D(6)), ('q0', R.poly_q(0))]:
    r = R.residues(S, f, d13, 14)
    for c in R.C.CLS:
        print(k, name, c, ' '.join('%.1e' % abs(x) for x in r[c]))
