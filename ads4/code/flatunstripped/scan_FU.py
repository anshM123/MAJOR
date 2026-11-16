"""[NE] floating-point scan of the heavy action of F_U (naive flat limit of Omega_g3 with C_k -> B'_k), normal form m^4 F:
  ++ : int [psi2 (2-uq) + u psi3] /(1-uq)^2 P_J dp + (10JJ - 20) u^3
  +- : int [psi2 (2-uq) - u psi3] /(1-uq)^2 d44 dp + (10JJ - 190) u^3,   d44 = (1-uq)^4 P^(0,8)_{J-4},  q = p^2, u = 1/m^2.
Compared with F_g3 (stripped, completed) on the same grid.  Not part of any proof.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(600); pg = (xg + 1) / 2; wg = wg / 2
psi2 = (1 - pg)**3 * pg * (65 + 155*pg + 47*pg**2); psi3 = 2*(1 - pg)**3*(81 - 1927*pg + 4601*pg**2)
def jac_all(nmax, a, b, x):
    P = np.zeros((nmax + 1, len(x))); P[0] = 1
    if nmax >= 1: P[1] = (a + 1) + (a + b + 2)*(x - 1)/2
    for n in range(2, nmax + 1):
        c = 2*n + a + b
        P[n] = (((c-1)*(a*a-b*b) + (c-2)*(c-1)*c*x)*P[n-1] - 2*(n+a-1)*(n+b-1)*c*P[n-2]) / (2*n*(n+a+b)*(c-2))
    return P
res = {}
for h in ['++', '+-']:
    for name in ['F_U', 'F_g3']:
        best = (1e9, 0, 0); bestb = (1e9, 0, 0)
        for m in np.concatenate([np.linspace(1, 2, 201), np.linspace(2, 40, 761), np.linspace(40, 100, 61)]):
            u = 1/m**2; q = pg**2; x = 1 - 2*u*q; Jmax = int(max(60, 3*m**2))
            if h == '++':
                D = jac_all(Jmax, 0, 0, x); Js = np.arange(0, Jmax + 1); s = 1
            else:
                D = jac_all(Jmax - 4, 0, 8, x); Js = np.arange(4, Jmax + 1); s = -1
            JJ = Js*(Js + 1)
            if name == 'F_U':
                wgt = (psi2*(2 - u*q) + s*u*psi3) / (1 - u*q)**2
                if h == '+-': wgt = wgt * (1 - u*q)**4          # d44 = (1-uq)^4 dtilde44
                fw = (10*JJ - (20 if h == '++' else 190))*u**3
            else:
                wgt = psi2*(2 - u*q) + s*u*psi3
                fw = (10*JJ - (0 if h == '++' else 210))*u**3
            v = D @ (wg*wgt) + fw
            if h == '++': v = v[::2]; Js = Js[::2]
            k = np.argmin(v)
            if v[k] < best[0]: best = (v[k], m, Js[k])
            bb = (2*Js + 1)/m; ww = v*np.maximum(bb, 1)**3; k = np.argmin(ww)
            if ww[k] < bestb[0]: bestb = (ww[k], m, Js[k])
        print('%s %s: min m^4F = %.4e at m=%.3f J=%d ; min max(1,b)^3 m^4F = %.4f at m=%.3f J=%d' % (name, h, *best, *bestb), flush=True)
