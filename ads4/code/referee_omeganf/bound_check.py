"""Independent recomputation of the conditional bound sqrt(t4^2+aT^2) <= 240 sqrt(n^8 16 L_AdS(Lam/n)/k)/Lam^4 for NF02,
from the exact rational psi2 in functionals_noforward.json (expanded, mpmath high precision; not the factored form)."""
import json, os
import mpmath as mp
mp.mp.dps = 60
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_noforward.json')))['NF02']
from fractions import Fraction as Fr
s2 = [Fr(c) for c in d['psi2_exact']]; s3 = [Fr(c) for c in d['psi3_exact']]
mu = lambda c, j: sum(ci/(i + j + 1) for i, ci in enumerate(c))
k0 = -2*mu(s2, 6) + 2*mu(s3, 4)
print("k0 =", k0, float(k0))
c2 = [mp.mpf(c.numerator)/c.denominator/ (mp.mpf(k0.numerator)/k0.denominator) for c in s2]
psi2 = lambda t: mp.polyval(c2[::-1], t)
x0 = c2[1]; print("x0 =", x0, " A = 16 x0 =", 16*x0)
Lfin = mp.quad(lambda t: (psi2(t) - x0*t)/t**2, [0, 1]); B = 16*Lfin; print("B =", B, " B_AdS =", B - 16*x0*mp.log(mp.mpf(3)/2))
def LAdS(Ls):
    return mp.quad(lambda t: psi2(t)/(t*t + mp.mpf(9)/4/Ls**2), [0, 1.5/Ls, 10/Ls, 1])
for Dg in (100, 1000):
    Lam = mp.sqrt((Dg - 2)*(Dg - 1))
    for n in (1.0, 1.13, 1.4, 1.5, 2.1):
        X = n**8*16*LAdS(Lam/n)
        print(f"Dg={Dg} n={n}: n^8 16L/k = {float(X):.5g}  bound(Lam=sqrt((D-2)(D-1))) = {float(240*mp.sqrt(X)/Lam**4):.3g}  bound(Lam->Dg) = {float(240*mp.sqrt(X)/Dg**4):.3g}")
# A-only asymptotic estimate: 16 L ~ A log(Lam/(1.5 n)) + B
