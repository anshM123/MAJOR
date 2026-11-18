"""
e6: conditional bound of OMEGA_U item 7 at the exact-skeleton threshold:  |ghat3|^2 Lam^8 <= n^8 16 L_AdS(Lam/n)/k(0),
sqrt(t4^2 + alpha_T^2) <= 240 sqrt(.)/Lam^4,  L_AdS(L) = int_0^1 psi2(p)/(p^2 + 9/(4L^2)) dp  (M02 psi2 exact, k(0) = 10884217/396000)
-- identical formula to omegau/u7_summary.py `bound`.  Also the leading-log form 16 L_AdS/k(0) ~ A log(Lam/n) + B_AdS (A = 51.74, B = -86.30).
Usage: python3 e6_bound.py Dg n [n ...]  -> printed (and appended to cache/e6.txt)
"""
import os, sys, json
from fractions import Fraction as Fr
import numpy as np
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
Jd = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
c = [float(Fr(x)) for x in Jd['M02']['psi2_exact']]; k0 = 10884217/396000
Dg = int(sys.argv[1]); ns = [float(a) for a in sys.argv[2:]]
Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
lines = []
for n in ns:
    L = Lam/n
    LA = quad(lambda t: sum(ci*t**i for i, ci in enumerate(c))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]
    X = n**8*16*LA/k0
    s = (f"Dg={Dg} n={n:.3f}: 16 L_AdS/k(0) = {16*LA/k0:.4g} (A log(Lam/n)+B = {51.7366*np.log(L) - 86.2981:.4g}); "
         f"|ghat3|^2 Lam^8 <= {X:.4g};  sqrt(t4^2+aT^2) <= {240*np.sqrt(X)/Lam**4:.3g}")
    print(s); lines.append(s)
open(os.path.join(HERE, 'cache', 'e6.txt'), 'a').write("\n".join(lines) + "\n")
