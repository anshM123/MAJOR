"""
poserr/curv_lemma.py (POS_ERRORS.md Lemma C): numerical check of the Mehler-Dirichlet comparison
   |P_{-1/2+i nu}(cosh b) - J0(2 nu sinh(b/2))| <= (1 - sech(b/2)) + min(2, nu (2 sinh(b/2) - b))     (all real nu >= 0, b > 0)
Proof (POS_ERRORS §2.1): P = (sqrt2/pi) int_0^b cos(nu t) dt / sqrt(cosh b - cosh t); substitute s = 2 sinh(t/2):
   P = (2/pi) int_0^B cos(nu t(s)) ds / (cosh(t/2) sqrt(B^2 - s^2)),  B = 2 sinh(b/2);  J0(nu B) = (2/pi) int_0^B cos(nu s) ds / sqrt(B^2 - s^2);
   |cos(nu t)/cosh(t/2) - cos(nu s)| <= (1 - sech(t/2)) + min(2, nu (s - t)), both increasing in t; (2/pi) int ds/sqrt(B^2-s^2) = 1.
Here: 50-digit mpmath values of P (Legendre function of the first kind, type 3) and J0 on a grid; prints the worst ratio |diff|/bound.
"""
import os, sys
import mpmath as mp
mp.mp.dps = 30
worst = (0, None)
rows = []
for b in [1e-3, 3e-3, 0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 4.0]:
    for nu in [0.0, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0]:
        if nu*b > 400:                                    # mpmath cost; the bound is 2+ there and |P|,|J0| <= 1 trivially
            continue
        P = mp.legenp(-mp.mpf(1)/2 + 1j*nu, 0, mp.cosh(b), type=3).real
        J = mp.besselj(0, 2*nu*mp.sinh(b/2))
        bd = (1 - 1/mp.cosh(b/2)) + min(2, nu*(2*mp.sinh(b/2) - b))
        r = abs(P - J)/bd if bd > 0 else 0
        rows.append((b, nu, float(P), float(J), float(abs(P - J)), float(bd), float(r)))
        if r > worst[0]:
            worst = (float(r), (b, nu))
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'curv_lemma.txt')
with open(out, 'w') as f:
    f.write('%8s %8s %14s %14s %12s %12s %8s\n' % ('beta', 'nu', 'P', 'J0', '|P-J0|', 'bound', 'ratio'))
    for r in rows:
        f.write('%8.3g %8.3g %14.6e %14.6e %12.4e %12.4e %8.4f\n' % r)
    f.write('worst ratio |P - J0|/bound = %.4f at (beta, nu) = %s  (must be <= 1)\n' % worst)
print(open(out).read().splitlines()[-1])
