"""t1: validation of the numerical kernels (floating point).
 (a) curvT against direct quadrature of the hypergeometric definition of P_{-1/2+i nu}(cosh beta);
 (b) W(beta), P0(beta) against mpmath;
 (c) flatT against direct J_0 quadrature, and against UNIFORM_POSITIVITY.md numbers (fhat(9.755) = 0.023158);
 (d) flat heavy action (FlatGrid) against the exact polynomial of FLAT_GRAVITON_BOUND.md Lemma P
     (19/315 at +-, J=4, m=M, completed; 3904/315 at ++, J=0, m=M), and the near-zero m^4F ~ 2.67e-5 at m~21.79, J=245.
"""
import sys
import numpy as np, mpmath as mp
from scipy.special import j0
from common import *

mp.mp.dps = 25
print('(a) curvT vs hypergeometric definition')
for Lam, beta in [(6.0, 1.0), (20.0, 0.3), (3.0, 2.5), (50.0, 0.05)]:
    c = curvT_all(Lam, beta)
    z = mp.sinh(beta / 2) ** 2
    for name, poly in [('psi2', PSI2), ('psi3', PSI3)]:
        f = lambda p: mp.mpf(float(peval(poly, float(p)))) * mp.re(mp.hyp2f1(0.5 - 1j * Lam * p, 0.5 + 1j * Lam * p, 1, -z))
        ref = mp.quad(f, mp.linspace(0, 1, 9))
        print('   Lam=%5.1f beta=%4.2f %s: %.12f  ref %.12f  diff %.1e' % (Lam, beta, name, c[name], ref, abs(c[name] - ref)))
print('(b) W, P0')
for beta in [0.1, 1.0, 3.0]:
    c = curvT_all(1.0, beta)
    z = mp.sinh(beta / 2) ** 2
    P = lambda nu2: mp.re(mp.hyp2f1(0.5 - 1j * mp.sqrt(nu2), 0.5 + 1j * mp.sqrt(nu2), 1, -z))
    Wref = -mp.diff(P, 0)
    print('   beta=%.1f  W=%.12f ref %.12f   P0=%.12f ref %.12f   sinh^2(b/2)=%.6f' % (beta, c['W'], Wref, c['P0'], P(0), float(z)))
print('(c) flatT')
xg, wg = np.polynomial.legendre.leggauss(400); pg = (xg + 1) / 2; wg = wg / 2
for b in [0.5, 9.755, 40.0]:
    f = flatT_all(b)
    ref = np.sum(wg * peval(PSI2, pg) * j0(b * pg))
    print('   b=%6.3f  fhat=%.10f  direct %.10f' % (b, f['psi2'], ref))
print('(d) flat heavy action')
G = FlatGrid(1600)
Js, B2, B3, FW, ab = G.parts('+-', 10, 1.0)
v = (B2 + B3 + FW)[0]
print('   +- J=4 m=1: %.10f   19/315 = %.10f   (abs norm %.3f)' % (v, 19 / 315, ab[0]))
Js, B2, B3, FW, ab = G.parts('++', 4, 1.0)
print('   ++ J=0 m=1: %.10f   3904/315 = %.10f' % ((B2 + B3 + FW)[0], 3904 / 315))
mm = 21.79
Js, B2, B3, FW, ab = G.parts('+-', 300, 1 / mm ** 2)
k = np.where(Js == 245)[0][0]
print('   +- J=245 m=21.79: m^4F = %.4e  (FLAT_GRAVITON_BOUND: 2.67e-5); parts B2=%.5f B3=%.5f FW=%.5f; abs norm %.4f'
      % ((B2 + B3 + FW)[k], B2[k], B3[k], FW[k], ab[k]))
