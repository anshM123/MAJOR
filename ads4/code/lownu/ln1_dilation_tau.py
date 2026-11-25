"""LF1/LF4: tau for the dilated NF02 (exact rational moments; L_A by mpmath quadrature at 50 digits).
tau_a(Lam) = 240 sqrt(16 L(a Lam)/k) / (a Lam)^4,   L(Lam) = int_0^1 psi2/(p^2 + 9/(4 Lam^2)) dp,  k = -2 mu2_6 + 2 mu3_4."""
import mpmath as mp
from ln_common import nf02_exact, mom_exact, guard
mp.mp.dps = 50
p2, p3 = nf02_exact()
k = -2 * mom_exact(p2, 6) + 2 * mom_exact(p3, 4)
P2 = [mp.mpf(x.numerator) / x.denominator for x in p2]
def L(Lam):
    e2 = mp.mpf(9) / (4 * mp.mpf(Lam) ** 2)
    f = lambda p: mp.polyval(P2[::-1], p) / (p * p + e2)
    return mp.quad(f, [0, mp.mpf(3) / (2 * Lam), 0.05, 0.3, 1])
def tau(Lam, a=1):
    Lp = a * mp.mpf(Lam)
    return 240 * mp.sqrt(16 * L(Lp) / k) / Lp ** 4
out = []
print('k (unnormalised) =', float(k))
for Dg in (300, 1000, 2000, 4000, 10000):
    Lam = mp.sqrt((Dg - 2) * (Dg - 1))
    row = [Dg, float(Lam), float(L(Lam) / k)]
    s = 'Dg=%6d Lam=%9.2f L_A/k=%9.3f tau(a=1)=%.5e' % (Dg, Lam, L(Lam) / k, tau(Lam))
    for a in (0.025, 0.05, 0.075, 0.1, 0.107, 0.118, 0.163, 0.35):
        t = tau(Lam, a); row.append(float(t)); s += ' | a=%.3f %.3e' % (a, t)
    print(s); out.append(row)
import json
json.dump(out, open(guard('cache/ln1_tau.json'), 'w'))
