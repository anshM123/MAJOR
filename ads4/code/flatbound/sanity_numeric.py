"""Floating-point sanity checks (NOT part of any proof).

 1. heavy action from the exact u-polynomial vs Gauss-Legendre quadrature with scipy Legendre/Jacobi values;
 2. Clebsch-Gordan identity (1-y)^8 dtilde44 = sum c_j P_{J+j}, at integer, non-integer and complex J;
 3. Mehler expansion H_w(L) = w(0)/b - w''(0)/(2b^3) + O(b^-4) and the size of b^4 * remainder;
 4. scan of the heavy action: minima of m^4 F (and of b^3 m^4 F) per channel, location of the near-zero.
"""
import numpy as np, mpmath as mp
from numpy.polynomial.legendre import leggauss
from scipy.special import eval_legendre, eval_jacobi
from fractions import Fraction as Fr
from common import poly_u, cg_coeffs, PHI_CONST

xg, wg = leggauss(800); pg = (xg + 1) / 2; wg = wg / 2
psi2 = (1 - pg) ** 3 * pg * (65 + 155 * pg + 47 * pg ** 2)
psi3 = 2 * (1 - pg) ** 3 * (81 - 1927 * pg + 4601 * pg ** 2)

def FF_quad(J, m, h, pess=True):
    u = 1 / m ** 2; x = 1 - 2 * u * pg ** 2
    D = eval_legendre(J, x) if h == '++' else eval_jacobi(J - 4, 0, 8, x)
    s = 1 if h == '++' else -1
    Phi = 10 * J * (J + 1) - (PHI_CONST[h] if pess else (5 if h == '++' else 205))
    return np.sum(wg * (psi2 * (2 - u * pg ** 2) + s * u * psi3) * D) + Phi * u ** 3

print('1. exact polynomial vs quadrature')
mx = 0
for h, J in [('++', 0), ('++', 10), ('++', 60), ('+-', 4), ('+-', 9), ('+-', 45)]:
    c = poly_u(J, h)
    for m in [Fr(1), Fr(13, 10), Fr(27, 10), Fr(9)]:
        u = 1 / m ** 2
        v1 = float(sum(ci * u ** i for i, ci in enumerate(c))); v2 = FF_quad(J, float(m), h)
        mx = max(mx, abs(v1 - v2) / max(1e-3, abs(v2)))
print('   max relative deviation %.2e' % mx)

print('2. CG identity')
mp.mp.dps = 40
for J in [mp.mpf(4), mp.mpf(7), mp.mpf('10.37'), mp.mpc('25.3', '4.1')]:
    c = cg_coeffs(J)
    for y in [mp.mpf('0.1'), mp.mpf('0.6'), mp.mpc('0.2', '0.15')]:
        lhs = (1 - y) ** 8 * mp.hyp2f1(4 - J, J + 5, 1, y)
        rhs = sum(c[j + 4] * mp.hyp2f1(-(J + j), J + j + 1, 1, y) for j in range(-4, 5))
        print('   J=%s y=%s  |lhs-rhs| = %.1e' % (mp.nstr(J, 5), mp.nstr(y, 3), float(abs(lhs - rhs))))

print('3. Mehler expansion, ++ channel weight at m=10 (b^4 * remainder)')
for J in [150, 300, 600, 1000]:
    m = 10.0; u = 1 / m ** 2; b = (2 * J + 1) / m
    I = np.sum(wg * (psi2 * (2 - u * pg ** 2) + u * psi3) * eval_legendre(J, 1 - 2 * u * pg ** 2))
    main = 162 * u / b + (80 - 21250 * u) / b ** 3
    print('   J=%d b=%.1f  b^4*(I-main) = %.1f   (certified bound C_R <= 6860)' % (J, b, (I - main) * b ** 4))

print('4. scan (pessimistic Phi): minima of m^4 F')
def jac_all(nmax, a, b, x):
    P = np.zeros((nmax + 1, len(x))); P[0] = 1
    if nmax >= 1: P[1] = (a + 1) + (a + b + 2) * (x - 1) / 2
    for n in range(2, nmax + 1):
        c = 2 * n + a + b
        P[n] = (((c - 1) * (a * a - b * b) + (c - 2) * (c - 1) * c * x) * P[n - 1] - 2 * (n + a - 1) * (n + b - 1) * c * P[n - 2]) / (2 * n * (n + a + b) * (c - 2))
    return P
for h in ['++', '+-']:
    best = (1e9, 0, 0); bestb3 = (1e9, 0, 0)
    for m in np.concatenate([np.linspace(1, 2, 101), np.linspace(2, 40, 761), np.linspace(40, 120, 81)]):
        Jmax = int(max(60, 3 * m ** 2)); u = 1 / m ** 2; x = 1 - 2 * u * pg ** 2
        if h == '++':
            D = jac_all(Jmax, 0, 0, x); Js = np.arange(0, Jmax + 1); s = 1
        else:
            D = jac_all(Jmax - 4, 0, 8, x); Js = np.arange(4, Jmax + 1); s = -1
        v = D @ (wg * (psi2 * (2 - u * pg ** 2) + s * u * psi3)) + (10 * Js * (Js + 1) - PHI_CONST[h]) * u ** 3
        if h == '++': v = v[::2]; Js = Js[::2]
        k = np.argmin(v)
        if v[k] < best[0]: best = (v[k], m, Js[k])
        bb = (2 * Js + 1) / m; w = v * np.maximum(bb, 1) ** 3; k = np.argmin(w)
        if w[k] < bestb3[0]: bestb3 = (w[k], m, Js[k])
    print('   %s: min m^4F = %.4e at m=%.3f J=%d ; min max(1,b)^3 m^4F = %.4f at m=%.3f J=%d' % (h, *best, *bestb3))

print('5. exterior Mehler formula P_L(cos t) = (sqrt2/pi) int_t^pi sin((L+1/2)f)/sqrt(cos t - cos f) df')
mp.mp.dps = 30
for L, t in [(0, 0.7), (3, 1.1), (17, 0.4), (40, 2.5)]:
    t = mp.mpf(t)
    val = mp.sqrt(2) / mp.pi * mp.quad(lambda f: mp.sin((L + mp.mpf(1) / 2) * f) / mp.sqrt(mp.cos(t) - mp.cos(f)), [t, (t + mp.pi) / 2, mp.pi])
    print('   L=%d theta=%s  |difference| = %.1e' % (L, mp.nstr(t, 3), float(abs(val - mp.legendre(L, mp.cos(t))))))
