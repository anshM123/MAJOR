"""Q3 flat action at m = M (u = 1) and near it, +- channel, large J (float, [NE]).  The stored psi2, psi3 are IEEE doubles; their
(1-p)^3 factor holds only to rounding: psi = (1-p)^3 q(p) + rem(p), rem of degree <= 2 computed EXACTLY (rationals).  Each part is evaluated
accurately (factored form; rem exactly), so the effect of rem (~1e-10) at the backward peak of d_44 (~J^8) is resolved.
Prints FF/N for the stored (binary) functional and for the 'exact-factor' functional (rem dropped)."""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'robust'))
from rcommon import nodes, flat_integrals
d = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
a = float(d['a']); c = float(d['c'])
def split(coeffs):
    P = [Fr(float(x)) for x in coeffs]
    q = list(P); quo = []
    # divide by (1-p)^3 = -(p-1)^3: synthetic division by (p - 1) three times, keep remainders exactly
    rem_poly = P
    import sympy as sp
    p = sp.Symbol('p'); Pe = sum(sp.Rational(x.numerator, x.denominator) * p ** i for i, x in enumerate(P))
    Q_, R_ = sp.div(sp.Poly(Pe, p), sp.Poly((1 - p) ** 3, p))
    return [float(x) for x in reversed(Q_.all_coeffs())], [sp.Rational(x) for x in reversed(R_.all_coeffs())]
q2, r2 = split(d['psi2']); q3, r3 = split(d['psi3'])
print('rem psi2 =', [float(x) for x in r2], ' rem psi3 =', [float(x) for x in r3])
def action(h, u, Jmax, with_rem=True):
    n = Jmax + 16; p, w = nodes(n)
    s = 1 if h == '++' else -1; c4 = 1 if h == '++' else 41
    f2 = (1 - p) ** 3 * np.polynomial.polynomial.polyval(p, q2); f3 = (1 - p) ** 3 * np.polynomial.polynomial.polyval(p, q3)
    if with_rem:
        f2 = f2 + sum(float(r) * p ** i for i, r in enumerate(r2)); f3 = f3 + sum(float(r) * p ** i for i, r in enumerate(r3))
    W = np.vstack([w * f2 * (2 - u * p * p), w * s * u * f3]).T
    Js, R = flat_integrals(h, u, Jmax, W, p)
    if h == '++': k = Js % 2 == 0; Js, R = Js[k], R[k]
    FW = a * (2 * Js * (Js + 1.0) - c4 + s) * u ** 3
    FF = R[:, 0] + R[:, 1] + FW - 2 * c * u * u; N = np.abs(R[:, 0]) + np.abs(R[:, 1]) + np.abs(FW) + 2 * c * u * u
    return Js, FF, N
Jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
for h in ['+-', '++']:
    for u in [1.0, 0.9999, 0.999, 0.99, 0.9, 0.6, 0.3]:
        out = []
        for wr in (True, False):
            Js, FF, N = action(h, u, Jmax, wr); q = FF / N; j = int(np.argmin(q))
            out.append('min FF/N = %+.5f at J = %d (first J with FF<0: %s)' % (q[j], Js[j], (int(Js[np.argmax(FF < 0)]) if np.any(FF < 0) else '-')))
        print('%s u=%-7g binary: %s | exact-factor: %s' % (h, u, out[0], out[1]), flush=True)
