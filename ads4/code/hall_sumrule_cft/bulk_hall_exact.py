"""
Bulk check (AdS4 black brane, parity-odd cubic coupling L > lam Q, Q = C C C~), conventions of
../scope_sumrule/hall_uv_sign.py (imported read-only; it prints its own derivation of the source on import).

(1) EXACT large-Omega constant.  At large Euclidean Omega the Hall response is localised at z = 1/r ~ 1/Omega, so the
    UV constant is the term linear in the black-brane mass m (the T-exchange: linear in <T>_beta).  Scaling z = zeta/Omega
    isolates it: the O(m) source is S1(z) = -144 Omega^3 z^4 (all other terms are O(m^2) = O(Omega^-3)).  With the
    pure-AdS regular mode Phi0 = 1 + Omega z, the sourced equation  z X'' - (2 + 2 Omega z) X' + 2 Omega X = -z S1
    has the unique polynomially bounded solution with X(0) = 0, X = -30 z^3 - 30 Omega z^4 - 18 Omega^2 z^5, so
        Pi(oo) = g''(0)/2 = -90,   S_H = Pi/2 = -45 per unit lam,   independent of Omega.
    (hall_uv_sign.py found -44.99 +- 0.03 by fitting.)  With lam = -alpha3~/6 and |alpha_T| = 240 |alpha3~|
    (DICTIONARY_T4.md sec. 4, graded [D] there), |S_H| = 7.5 |alpha3~| = |alpha_T|/32.

(2) Convergent (spectral) form, real frequency.  delta h^R(w)/eps = Pi(Omega = -i w)/2.  Check
        S_H = (1/(pi eps)) Int_{-oo}^{oo} eta_H(w) dw,   eta_H(w) = -Im delta h^R(w)/w,
    numerically, including eta_H(0)/eps = 54 (= 4 pi eta_H/s = 108 with s = 2 pi eps on the brane).
"""
import os
import sys
sys.dont_write_bytecode = True  # do not leave __pycache__ in the imported (read-only) directories

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'scope_sumrule'))
print('--- importing ../scope_sumrule/hall_uv_sign.py (its own output follows) ---')
import hall_uv_sign as H  # noqa: E402
print('--- end of import ---\n')

z, Om, zeta = H.z, H.Om, sp.Symbol('zeta', positive=True)

# ---------------- (1) exact UV constant
Phi0 = 1 + Om*z
ode0 = lambda y: z*sp.diff(y, z, 2) - (2 + 2*Om*z)*sp.diff(y, z) + 2*Om*y   # pure-AdS part of hall_uv_sign's operator
assert sp.simplify(ode0(Phi0)) == 0
S = sp.expand(sum(c*sp.diff(Phi0, z, k) for k, c in enumerate(H.coef)))
Ssc = sp.expand(S.subs(z, zeta/Om))
lead = sp.limit(Ssc*Om, Om, sp.oo)
S1 = sp.expand(lead.subs(zeta, Om*z)/Om)
rest = sp.expand(Ssc - lead/Om)
print('source with Phi0 = 1 + Omega z:  S =', sp.factor(S))
print('O(m) part (leading at fixed zeta = Omega z):  S1 =', S1, ';  remainder is O(Omega^%s)'
      % sp.degree(sp.expand(rest.subs(Om, 1/sp.Symbol('t'))), sp.Symbol('t')).__neg__())
a = sp.symbols('a0:8')
X = sum(a[k]*z**k for k in range(8))
eqs = sp.Poly(sp.expand(ode0(X) + z*S1), z).all_coeffs() + [X.subs(z, 0)]
sol = sp.solve(eqs, a, dict=True)
assert len(sol) == 1
Xs = sp.expand(X.subs(sol[0]))
print('polynomial (regular) solution X =', Xs)
g = sp.diff(Xs, z) - Om*Xs                                     # g = (1 - z^3) X' - Omega X at O(m)
Pi_inf = sp.simplify(sp.diff(g, z, 2).subs(z, 0)/2)
print('g(0) =', g.subs(z, 0), ", g'(0) =", sp.diff(g, z).subs(z, 0), ',  Pi(oo) = g\'\'(0)/2 =', Pi_inf,
      '=>  S_H = Pi/2 =', Pi_inf/2, 'per unit lam (exact, Omega-independent)')
SHa = (Pi_inf/2)*sp.Rational(-1, 6)
print('with lam = -alpha3~/6:  S_H = %s alpha3~ ;  |alpha_T| = 240|alpha3~|  =>  |S_H| = |alpha_T|/%s'
      % (SHa, sp.Integer(240)/abs(SHa)))
print('approach to the limit (full Chebyshev numerics of hall_uv_sign.solve, N = 120; roundoff ~1e-2 at Omega = 40):')
for o in (5.0, 10.0, 20.0, 40.0):
    print('   Omega = %-5g Pi(oo) = %.5f' % (o, H.solve(o, 120)[2]))

# ---------------- (2) spectral form on the real axis
print('\nreal frequency: delta h^R(w)/eps = Pi(-i w)/2 ;  eta_H(w)/eps = -Im Pi(-i w)/(2 w)')


def etaH(w, N=80):
    return -np.imag(H.solve(-1j*w, N)[2])/(2*w)


# Chebyshev roundoff grows with N (derivatives at z = 0), so N = 60-90 is optimal; above w ~ 30 the signal
# (eta_H ~ 1e-3) is below the roundoff floor, so the integral is cut at W = 24 and the tail fitted to c/w^4
# (next thermal operator [TT], Delta = 6).
for w in (1e-3, 0.5, 1, 2, 4, 8, 12, 16, 24):
    print('   w = %-6g eta_H/eps = %+.8f   (N=60)   %+.8f (N=90)' % (w, etaH(w, 60), etaH(w, 90)))
W = 24.0
edges = [0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 12, 16, W]
xg, wg = np.polynomial.legendre.leggauss(24)
tot = 0.0
for lo, hi in zip(edges[:-1], edges[1:]):
    ws = 0.5*(hi - lo)*xg + 0.5*(hi + lo)
    tot += 0.5*(hi - lo)*np.dot(wg, np.array([etaH(w) for w in ws]))
tw = np.array([12.0, 14.0, 16.0, 20.0, 24.0])
tv = np.array([etaH(w) for w in tw])
print('   tail: eta_H/eps * w^4 at w = 12..24:', np.round(tv*tw**4, 1))
c4 = np.mean(tv*tw**4)
tail = c4/(3*W**3)
SH = (2/np.pi)*(tot + tail)
print('   (1/(pi eps)) Int_R eta_H dw = (2/pi)[%.6f + tail %.4f] = %.4f   (exact UV constant: -45)' % (tot, tail, SH))
print('   eta_H(0)/eps -> 54 (4 pi eta_H/s = 108), yet the integral is negative: eta_H(w) changes sign; no positivity.')
