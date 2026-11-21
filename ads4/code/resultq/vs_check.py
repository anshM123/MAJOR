"""Virasoro-Shapiro (type II superstring) MHV graviton amplitude in CHLPS conventions (arXiv:2201.06602 (4.2), p = 0):
f = 8 pi G/(stu) * G(1-s)G(1-t)G(1-u)/(G(1+s)G(1+t)G(1+u)), M = 1 = first massive level (M^2 = 4/alpha').
Checks: (i) g4/(8 pi G) = 2 zeta(3), g5 = 0, |g3| = 0; (ii) residue of f at s = 1 is t-independent (pure J = 4, +- channel, positive)
and the t = 1 residue is s-independent (pure J = 0, ++ channel); (iii) placement against the Stage-0 bounds."""
import json, os
import mpmath as mp
mp.mp.dps = 30
R = lambda s, t: mp.gamma(1 - s) * mp.gamma(1 - t) * mp.gamma(1 + s + t) / (mp.gamma(1 + s) * mp.gamma(1 + t) * mp.gamma(1 - s - t))
for eps in [mp.mpf('1e-3'), mp.mpf('1e-4')]:
    s, t = eps, 2 * eps; u = -s - t
    print('(R-1)/(stu) at eps=%s: %s   2 zeta3 = %s' % (eps, mp.nstr((R(s, t) - 1) / (s * t * u), 12), mp.nstr(2 * mp.zeta(3), 12)))
# g5: f = 8piG/(stu) + g4 + g5 t + ... ; fully symmetric -> linear term vanishes; check numerically
s, t = mp.mpf('1e-3'), mp.mpf('-3e-3'); u = -s - t
c1 = ((R(s, t) - 1) / (s * t * u) - 2 * mp.zeta(3))
print('next term size (should be O(eps^2), no linear t term):', mp.nstr(c1, 5))
# residues at s = 1 (u = -1 - t): Res f ~ 8 pi G * (-1)/(t u) * G(1-t)G(1-u)/(G(1+t)G(1+u))
for tt in [mp.mpf('-0.3'), mp.mpf('-0.7')]:
    uu = -1 - tt
    res = (1 / (tt * uu)) * mp.gamma(1 - tt) * mp.gamma(1 - uu) / (mp.gamma(1 + tt) * mp.gamma(1 + uu))
    print('s=1 residue/(8 pi G) at t=%s: %s' % (tt, mp.nstr(res, 12)))
z3 = float(2 * mp.zeta(3))
from math import log
rows = {}
for name, A, B in [('robust full (d12, delta .05)', 11.448, -12.979), ('robust minimal no-B3 (d12, delta .05)', 15.005, -15.973),
                   ('optimum full (d12, delta 0)', 10.112, -11.476), ('robust full + light matter', 14.268, -15.576),
                   ('CHLPS optimal (4.4) w/ matter', 12.3, -13.5), ('CHLPS example F_g4', 14.143, -15.507)]:
    vals = [A * log(2 * L / 3) + B for L in (100, 300, 1000)]
    Lc = 1.5 * mp.e ** ((z3 - B) / A)
    print('%-40s bound at Lam=100/300/1000: %6.2f %6.2f %6.2f | 2zeta3/bound = %.3f %.3f %.3f | VS on boundary at Lam=%.2f'
          % (name, *vals, *[z3 / v for v in vals], float(Lc)))

# (iv) next contacts: f(8piG)^-1 = 1/(stu) + g4 + g5 t + g6 t^2 - g6' s u + ...  (VS: expect g6 = g6' = 2 zeta5, g5 = 0)
def contact(s, t):
    u = -s - t
    return (R(s, t) - 1) / (s * t * u)
e = mp.mpf('1e-4')
c0 = contact(e * 0, e) if False else None
# fit g6, g6' from two directions at small scale
def quad(s, t):
    return (contact(s, t) - 2 * mp.zeta(3))
x = mp.mpf('1e-3')
q1 = quad(x, x)        # t = x, s = x, u = -2x: g6 x^2 - g6' (x)(-2x) = (g6 + 2 g6') x^2
q2 = quad(-2 * x, x)   # s = -2x, t = x, u = x:  g6 x^2 - g6' (-2x^2) = (g6 + 2 g6') x^2  (same); use t = -2x, s = x, u = x:
q3 = quad(x, -2 * x)   # t = -2x: 4 g6 x^2 - g6' x^2
g6p = (4 * q1 - q3) / (9 * x * x); g6 = (q1 / (x * x) - 2 * g6p)
print('VS: g6/(8piG) = %s  g6p/(8piG) = %s  (2 zeta5 = %s);  g6 M^4/g4 = %.4f (forward bound <= 1);  g6p/g6 = %.4f (CHLPS: 0 <= . <= 3)'
      % (mp.nstr(g6, 8), mp.nstr(g6p, 8), mp.nstr(2 * mp.zeta(5), 8), float(g6 / (2 * mp.zeta(3))), float(g6p / g6)))
