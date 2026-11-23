"""Float prototype (NOT a certificate) for the NF02 +- corner (J >= 400, 3/8 <= u <= 1): split w = w chi + w (1 - chi),
chi = 1 on [0,a], C^4 smoothstep 1 -> 0 on [a,c], 0 on [c,1].  Estimates
  CR  = Lemma-M constant of w_A (1-up^2)^-8 (formula of FLAT_GRAVITON_BOUND §3 (M)),
  T   = int |L^2 w_B| min(Bn, (1-up^2)^-4) dp  (crude bound on (4 u lam)^2 |I_B|),
  gam = -w''(0)/2,  and the margin  gam b/(1+ds)^3 - CR/(1-ds)^4 - (b^2/(b^2-81u))^2 T  at J = 401."""
import json, sys
import numpy as np, sympy as sp
from math import comb
from numpy.polynomial import Polynomial as Pl
from fractions import Fraction as Fr
D = json.load(open('/home/user/MAJOR/ads4/code/robust/functionals_noforward.json'))['NF02']
Xp = Pl([float(Fr(c)) for c in D['X']]); Yp = Pl([float(Fr(c)) for c in D['Y']])
P_ = Pl([0, 1])
ps, us = sp.symbols('p u'); F = sp.Function('f')
Ls = lambda g: sp.diff((sp.diff(g, ps) - g / ps) * (1 - us * ps**2) + 16 * us * ps * g, ps)
KL = int(sys.argv[2]) if len(sys.argv) > 2 else 3
LK = F(ps)
for _ in range(KL): LK = Ls(LK)
LK = sp.expand(LK); LKm1 = F(ps)
for _ in range(KL - 1): LKm1 = Ls(LKm1)
dLKm1 = sp.expand(sp.diff(LKm1, ps))
dsyms = sp.symbols('f0:%d' % (2 * KL + 1))
sub = lambda e: e.subs({sp.Derivative(F(ps), (ps, k)): dsyms[k] for k in range(2 * KL, 0, -1)}).subs(F(ps), dsyms[0])
L2f = sp.lambdify((ps, us) + dsyms, sub(LK), 'numpy')
dLm1f = sp.lambdify((ps, us) + dsyms, sub(dLKm1), 'numpy')
NDER = 2 * KL + 1
rs = [sp.lambdify((ps, us), sp.diff((1 - us * ps**2) ** -8, ps, k), 'numpy') for k in range(5)]
def smooth(a, c, k=5):
    t = Pl([-a / (c - a), 1 / (c - a)]); S = Pl([0])
    for j in range(k, 2 * k): S = S + comb(2 * k - 1, j) * t ** j * (1 - t) ** (2 * k - 1 - j)
    return S
def omega(u, s2, s3, dl):
    return 2 * (1 - dl * s2) * Xp + u * (-(1 - dl * s2) * P_ ** 2 * Xp - (1 - dl * s3) * Yp)
def est(J, u, a, c, s2, s3, dl, NP=20001):
    om = omega(u, s2, s3, dl); w = P_ * om * Pl([1, -1]) ** 5
    S = smooth(a, c)
    gam = -w.deriv(2)(0.0) / 2
    # C_R of w_A (1-up^2)^-8
    n4 = npw3 = n2 = s2m = s1m = 0.0
    for lo, hi, f in [(0, a, w), (a, c, w * (1 - S))]:
        pp = np.linspace(lo, hi, NP)
        fd = [f.deriv(k)(pp) if k else f(pp) for k in range(5)]
        rd = [np.asarray(rs[k](pp, u), float) * np.ones_like(pp) for k in range(5)]
        g = [sum(comb(k, i) * fd[i] * rd[k - i] for i in range(k + 1)) for k in range(5)]
        n4 += np.trapezoid(np.abs(g[4]), pp); npw3 += np.trapezoid(np.abs(pp * g[3]), pp); n2 += np.trapezoid(np.abs(g[2]), pp)
        s2m = max(s2m, np.abs(g[2]).max()); s1m = max(s1m, np.abs(g[1]).max())
    CR = 2 / np.pi * (2 / 3 * n4 + 6 * u * npw3 + 4 * u * n2) + u ** 1.5 * (3 * np.pi / 4 * s2m + s1m)
    # T: int |L^2 w_B| E ; near p = 1 use a geometric grid; factored evaluation of w = p (1-p)^5 om
    n = J - 4; Bn = comb(n + 8, 8)
    T = 0.0
    for lo, hi, piece in [(a, c, 'mid'), (c, 1.0, 'hi')]:
        if piece == 'mid':
            pp = np.linspace(lo, hi, NP); f = w * S
            fd = [f.deriv(k)(pp) if k else f(pp) for k in range(NDER)]
        else:
            dd = np.unique(np.concatenate([np.geomspace(1e-10, 1 - lo, NP), [1 - lo]])); pp = 1 - dd[::-1]; dd = dd[::-1]
            # derivatives of p (1-p)^5 om(p) in factored (delta) form via product rule
            A = [P_ * om] ; [A.append(A[-1].deriv()) for _ in range(NDER)]
            Av = [Ak(pp) for Ak in A]
            from math import perm
            Bv = [(-1) ** k * perm(5, k) * dd ** (5 - k) if k <= 5 else 0 * dd for k in range(NDER)]   # d^k/dp^k (1-p)^5
            fd = [sum(comb(k, i) * Av[i] * Bv[k - i] for i in range(k + 1)) for k in range(NDER)]
            # boundary term at p = 1 (k-th integration by parts): w_{K-1}'(1) (1-u) y(1)
            f1 = [A[0](1.0) * 0] * NDER
            Av1 = [Ak(1.0) for Ak in A]
            fd1 = [sum(comb(k, i) * Av1[i] * ((-1) ** (k - i) * perm(5, k - i) if k - i == 5 else 0) for i in range(k + 1)) for k in range(NDER)]
            bnd = abs(dLm1f(1.0, u, *fd1)) * (1 - u) * min(Bn, (1 - u) ** -4.0 if u < 1 else Bn) if KL % 2 == 1 or True else 0
        v = np.abs(L2f(pp, u, *fd)) * np.minimum(Bn, (1 - u * pp * pp) ** -4.0)
        T += np.trapezoid(v, pp)
    b = (2 * J + 1) * np.sqrt(u); ds = 8 * np.sqrt(u) / b
    fac = (b * b / (b * b - 81 * u)) ** 2 / (b * b - 81 * u) ** (KL - 2)
    T = fac * (T + (bnd if KL >= 3 else 0))
    marg = gam * b / (1 + ds) ** 3 - CR / (1 - ds) ** 4 - T
    return gam, CR, T, b, marg
if __name__ == '__main__':
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 401
    for a, c in [(0.3, 0.6), (0.4, 0.7), (0.3, 0.7), (0.5, 0.75), (0.2, 0.6)]:
        worst = None
        for s2 in (1, -1):
            for s3 in (1, -1):
                for u in [0.375, 0.6, 0.8, 0.95, 0.99, 0.999, 1.0]:
                    g, CR, T, b, m = est(J, u, a, c, s2, s3, 0.2)
                    if worst is None or m / (g * b) < worst[0]:
                        worst = (m / (g * b), s2, s3, u, g * b, CR, T)
        print('a=%.2f c=%.2f  worst margin/(gam b) = %.3f at s2=%+d s3=%+d u=%.3f  (gam b=%.3g CR=%.3g T=%.3g)' % ((a, c) + worst), flush=True)
