"""Referee: conditional AdS bound arithmetic n^8 * 16 L_AdS(Lam/n)/k with exact-ish quadrature (mpmath)."""
import mpmath as mp
from heavy_indep import weights
from fractions import Fraction as Fr
def val(name, n, Dg):
    psi2, psi3, a = weights(name)
    k = -2*sum(Fr(c)/(i+7) for i, c in enumerate(psi2)) + 2*sum(Fr(c)/(i+5) for i, c in enumerate(psi3)) + 4*a
    Lam = mp.sqrt((Dg-2)*(Dg-1)); Ls = Lam/n; aa = mp.mpf(9)/(4*Ls**2)
    f = lambda p: sum(mp.mpf(c.numerator)/c.denominator*p**i for i, c in enumerate(psi2))/(p*p+aa)
    L = mp.quad(f, [0, mp.sqrt(aa), 10*mp.sqrt(aa), 1])
    X = n**8*16*L/(mp.mpf(k.numerator)/k.denominator)
    return X, 240*mp.sqrt(X)/Lam**4
for nm, n in [('M02', 2.2), ('M02', 2.357), ('Fg3', 8.03), ('Fg3', 1.0), ('M01', 3.2)]:
    for Dg in (100, 1000):
        if nm == 'Fg3' and n == 8.03 and Dg == 1000: n2 = 22.1
        else: n2 = n
        X, t = val(nm, n2, Dg); print('%s n=%.3f Dg=%d: |g3|^2 Lam^8 <= %.4g ; sqrt(t4^2+aT^2) <= %.3g' % (nm, n2, Dg, X, t))
