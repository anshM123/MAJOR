import sys, math
sys.path.insert(0, '../analytic')
import crossing as X
import numpy as np
sys.path.insert(0, '.')
from fastS import S_vec, d_of
from common2 import n_of_L
tau = float(sys.argv[1]); sig = float(sys.argv[2])
G = lambda a: tau*math.acos(min(1, a/tau)) - a*math.acosh(max(1, tau/a))
print('G_998 =', G(998.0))
for c in (0.5, 0.7, 1.0, 1.5, 2.0):
    print('c', c, 'log bound', X.log_cross(sig, tau, 998.0, c))
for a in (998.0, 998.5, 1000.0, 1005., 1010.0, 1020.0, 1030.0, 1040.0):
    row = []
    for L in (3.0, 4.0, 4.5, 5.0, 6.0, 8.0):
        n = n_of_L(L, a); d = d_of(n, a, 1)
        S, lg = S_vec(n, a, 1, [complex(sig, tau)])
        row.append(math.log(abs(S[0])) + lg[0] + sig*math.log(d))
    print(f'a={a}: G_a={G(a):.3f}  log|S^| at L=3,4,4.5,5,6,8:', ' '.join(f'{x:.3f}' for x in row))
