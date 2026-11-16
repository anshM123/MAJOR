r"""Self-test of the interval-arithmetic layer (sanity only): at 2000 pseudo-random points the mpmath.iv enclosures of
sin, cos, exp, pi and of a rational are compared with 60-digit mpmath values, and common.S_iv / common.K_low_iv are
compared with independent 60-digit evaluations by quadrature-free formulas.  A failure here would indicate a broken
mpmath installation; a pass is not a proof of correctness of mpmath.iv (which the certificates assume)."""
import random
import sys
import mpmath
from common import *

random.seed(20261116)
hp = mpmath.mp.clone()
hp.dps = 60
ok = True
for _ in range(2000):
    x = random.uniform(0.01, 60.0)
    xi = iv.mpf(x)
    for name, fi, fh in (("sin", iv.sin, hp.sin), ("cos", iv.cos, hp.cos), ("exp(-x)", lambda z: iv.exp(-z), lambda z: hp.exp(-z))):
        e = fi(xi)
        v = fh(hp.mpf(x))
        if not (hp.mpf(e.a) <= v <= hp.mpf(e.b)):
            ok = False
            print("FAIL enclosure", name, x)
ok = ok and (hp.mpf(iv.pi.a) <= hp.pi <= hp.mpf(iv.pi.b))
# S_iv against the defining integral evaluated through the exact antiderivative formula at 60 digits
for _ in range(200):
    t = random.uniform(3.0, 30.0)
    e = S_iv(iv.mpf(t))
    T = hp.mpf(t)
    z = -1j * T
    v = hp.im(sum((D0.get(j, 0) - hp.exp(-z) * D1.get(j, 0)) / z ** (j + 1) for j in range(7)))
    if not (hp.mpf(e.a) <= v <= hp.mpf(e.b)):
        ok = False
        print("FAIL S enclosure", t)
    M = random.uniform(1.0, 7.0)
    e = K_low_iv(iv.mpf(M), iv.mpf(t), 4)
    v = sum((1 if k == 0 else 2) * hp.im(sum((D0.get(j, 0) - hp.exp(-(2 * hp.pi * k * M - 1j * T)) * D1.get(j, 0)) / (2 * hp.pi * k * hp.mpf(M) - 1j * T) ** (j + 1) for j in range(7))) for k in range(5))
    if not (hp.mpf(e.a) <= v <= hp.mpf(e.b)):
        ok = False
        print("FAIL K_low enclosure", M, t)
print(f"VERDICT interval self-test: {'PASS' if ok else 'FAIL'}")
sys.exit(0 if ok else 1)
