"""t5: the would-be light side of Omega_g3 IF the spinning graviton-exchange contribution to C^+_{2,nu} is the scalar one of
CHMRS (5.x) with 8 pi G -> 16 pi G:  -C^+_{2,nu}|_grav = 16 pi G/(nu^2 + 9/4)  [I: not computed for <TTTT>].
Then the B2 graviton term of the flat Theorem 2, 16 pi G L(eps), is replaced by 16 pi G L_AdS(Lam),
   L_AdS(Lam) = int_0^1 psi2(p)/(p^2 + a^2) dp,  a = 3/(2 Lam)  (exact quadrature)  = 65 log(1/a) - 1171/15 + O(a)  [D],
and the would-be bound (heavy positivity, k_AdS = k(0) and no other light terms assumed) reads
   |g3|^2 Lam^8 <= 16 L_AdS(Lam)/k(0),  k(0) = 330151/12012.
"""
import numpy as np
from scipy.integrate import quad
from common import PSI2, peval

k0 = 330151 / 12012
print('k(0) = %.6f ; 16*65/k(0) = %.4f ; 16*(1171/15)/k(0) = %.4f' % (k0, 16 * 65 / k0, 16 * 1171 / 15 / k0))
for Lam in [10, 30, 100, 300, 1000, 1e4]:
    a = 1.5 / Lam
    L = quad(lambda p: peval(PSI2, p) / (p * p + a * a), 0, 1, points=[a, 10 * a], limit=400, epsabs=1e-13)[0]
    asym = 65 * np.log(1 / a) - 1171 / 15
    print('Lam=%7g  L_AdS=%.6f  65log(2Lam/3)-1171/15=%.6f  diff=%.2e   16L/k0 = %.4f = 37.8387 log(Lam) %+.4f'
          % (Lam, L, asym, L - asym, 16 * L / k0, 16 * L / k0 - 16 * 65 / k0 * np.log(Lam)))
