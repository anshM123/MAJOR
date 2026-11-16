"""Conditional AdS4 bound constants (ROBUST_FUNCTIONAL.md section 5.4), given heavy positivity for m >= n Lam_s with Lam_s = Lam/n (i.e. for all
heavy states), the light side of ADS4_LIGHT.md (A = 16 x0/k, exact L_AdS, k_AdS = k(0), same inputs (F),(R),(P)):
    |ghat3|^2 Lam^8 <= n^8 * 16 L_AdS(Lam/n)/k  =  n^8 [A log(Lam/n) + B_AdS + O(n/Lam)]  =:  A' log Lam + B' + ...
A' = n^8 A, B' = n^8 (B_AdS - A log n).  Floating point (the LP-derived n are [NE]; the L_AdS integrals are elementary)."""
import numpy as np
from check import load
def LAdS(F, Lams):
    from numpy.polynomial.legendre import leggauss
    from rcommon import peval
    a = 1.5 / Lams
    edges = np.unique(np.concatenate([[0.0], a * np.logspace(-3, np.log10(1 / a), 80), [1.0]])); edges = edges[edges <= 1]
    x, w = leggauss(40); tot = 0.0
    for lo, hi in zip(edges[:-1], edges[1:]):
        p = lo + (x + 1) * (hi - lo) / 2; tot += np.sum(w * (hi - lo) / 2 * peval(F.psi2, p) / (p * p + a * a))
    return tot
rows = [('Fg3', 1.0, 'F_g3 if positive down to m = Lam (ADS4_LIGHT.md)'), ('Fg3', 8.03, 'F_g3, Regge threshold C1=1, Dg<=300'),
        ('Fg3', 22.1, 'F_g3, Regge threshold C1=1, Dg=1000'), ('Fg3', 10.05, 'F_g3, Regge threshold C1=10, Dg=100'),
        ('M01', 3.2, 'M01, C1=1 [NE threshold]'), ('M01', (100/9)**0.5, 'M01, C1=1, certified-flat threshold'),
        ('M02', (50/9)**0.5, 'M02, C1=1, certified-flat threshold (50/9)^1/2'), ('M02', (150/9)**0.5, 'M02, C1=3, certified-flat'),
        ('M02', (500/9)**0.5, 'M02, C1=10, certified-flat'), ('M02', 2.2, 'M02, C1=1 [NE threshold]'), ('M02', 4.0, 'M02, C1=3'), ('M02', 6.8, 'M02, C1=10'),
        ('M01', 5.2, 'M01, C1=3'), ('M01', 7.9, 'M01, C1=10'),
        ('T1', 1.75, 'T1 (targeted C1=1)'), ('T3', 2.25, 'T3 (targeted C1=3)'), ('T10', 4.0, 'T10 (targeted C1=10)')]
print('%-50s %6s %10s %10s %12s %12s %12s %12s' % ('case', 'n', "A'", "B'", 'Dg=20', 'Dg=100', 'Dg=1000', 'sqrt(t4^2+aT^2),Dg=100'))
for name, n, lab in rows:
    F = load(name); k = F.k(); A, B, BA = F.AB()
    Ap = n ** 8 * A; Bp = n ** 8 * (BA - A * np.log(n))
    vals = []
    for Dg in [20, 100, 1000]:
        Lam = np.sqrt((Dg - 2.0) * (Dg - 1.0))
        vals.append(n ** 8 * 16 * LAdS(F, Lam / n) / k)
    L100 = np.sqrt(98.0 * 99.0)
    print('%-50s %6.3f %10.4g %10.4g %12.4g %12.4g %12.4g %12.3g' % (lab, n, Ap, Bp, *vals, 240 * np.sqrt(vals[1]) / L100 ** 4))
