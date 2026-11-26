# r8_flatJ0.py -- does the heavy angular formula continued to a J=0 state reproduce the paper's flat light-scalar
# polynomial F_D^(0) = 11.317027 + 20.771893 m^2 (supplement, sec:light, Omega^EE)?  Heavy formula at J=0 (b=0, n=0):
# B(1-r)^3 V = (2 pi G Lam^4/9)[F_0 (T+9/2) + F_1/2] + k=3 term (2 pi G Lam^4/9) Lam^2 Q_0, Q_0 = int psi3.
# Flat limit (Lam=1 units, T->m^2 Lam^2): F_0 m^2 + F_1/(2Lam^2) + int psi3,  F_0 = 2 int psi2,  F_1/(2 Lam^2) -> -int psi2 p^2.
import json
from fractions import Fraction as Fr
d = json.load(open('/home/user/MAJOR/ads4/prl/data/nf02_psi.json'))
k0 = Fr(d['k0']); p2 = [Fr(x) for x in d['psi2_exact']]; p3 = [Fr(x) for x in d['psi3_exact']]
mom = lambda c, j: sum(ci/(n+j+1) for n, ci in enumerate(c))/k0
print('2 int psi2               =', float(2*mom(p2, 0)), '  (paper m^2 coefficient 20.771893)')
print('-int psi2 p^2 + int psi3 =', float(-mom(p2, 2)+mom(p3, 0)), '  (paper constant 11.317027)')
