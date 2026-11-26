"""
fullpos4/fp4_s55.py (FULLPOS4.md §1) [DV]: the author's §55 contour, literal recipe (coordinator's message, 26 Nov).
  x = e^{i al}, y = e^{i ga}, xi = (y-x)/(xy-1) (rational; no fresh sqrt), z = xi^2, beta = Log((1+xi)/(1-xi)).
  W(al, ga) = xi^{s-1} P(beta) / (sin al (1-xi^2)^p) e^{i S al - i a ga}
  I_U: al = A - i eps, ga = C - i eps,  0 <= A <= C <= pi/2      (dC dA)
  I_T: al = A - i eps, ga = pi/2 - i t, 0 <= A <= pi/2, t >= eps  (dt dA)
  I_L: al = A - i eps, ga = A - i(eps + v), -pi/2 <= A <= 0, v >= 0 (dv dA)
  nu > 1 (outgoing P+):       K_p = (C/pi) Re(I_U - i I_T - i I_L)
  nu <= 1 (complete conical): K_p = (C/2pi)(I_U - i I_T - i I_L)
Diagonal: C = A + (pi/2 - A) t^4, v = V t^4 (Jacobians), never P+(0).
Quadrature: Gauss-Legendre panels in A (graded at A = 0) and in t.
Usage: fp4_s55.py n N b p nu [nA] [nt]
"""
import sys, math, cmath
sys.dont_write_bytecode = True
import numpy as np
import mpmath as mp
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from fp4_pf import Pplus, Pfull, gl_panels, graded, Cconst
import fp4_kern as K


def W(al, ga, p, nu, s, S, a, outgoing):
    x = cmath.exp(1j*al); y = cmath.exp(1j*ga)
    xi = (y - x)/(x*y - 1)
    if outgoing and abs(xi) < 1e-13:      # t^4-map: these nodes carry weight O(t^3) ~ 1e-10 x log; dropped
        return 0j
    beta = cmath.log((1 + xi)/(1 - xi))
    P = Pplus(nu, beta) if outgoing else Pfull(nu, beta)
    return xi**(s - 1)*P/(cmath.sin(al)*(1 - xi*xi)**p)*cmath.exp(1j*S*al - 1j*a*ga)


def pieces(n, N, b, p, nu, nA=24, nt=24, V=None):
    s = 2*b + 1; a = N - n; J = n + b; S = 2*J + 1 + a; eps = 1.0/S
    outgoing = nu > 1
    V = V or 60.0/a
    # A grid: graded at 0 (1/sin al with |al| >= eps)
    Ap, wAp = gl_panels(graded(0.0, math.pi/2, eps*1e-2, 1.5), nA)
    tq, wtq = gl_panels(graded(0.0, 1.0, 1e-3, 1.5), nt)          # t in [0,1] for the t^4 maps
    IU = IT = IL = 0j; aU = aT = aL = 0.0
    for A, wa in zip(Ap, wAp):
        al = A - 1j*eps
        # I_U: C = A + (pi/2 - A) t^4
        L = math.pi/2 - A
        sU = 0j; sa = 0.0
        for t, wt in zip(tq, wtq):
            C = A + L*t**4; jac = 4*L*t**3
            if t == 0 or jac == 0:
                continue
            v = W(al, C - 1j*eps, p, nu, s, S, a, outgoing)*jac*wt
            sU += v; sa += abs(v)
        IU += sU*wa; aU += sa*wa
        # I_T: t >= eps, map t = eps + u, u in [0, oo): graded panels to 60/a + 2
        uq, wuq = gl_panels(graded(0.0, 80.0/max(a - nu, 1.0) + 2.0, 1e-3, 1.5), nt)
        sT = 0j; sa = 0.0
        for u, wu in zip(uq, wuq):
            v = W(al, math.pi/2 - 1j*(eps + u), p, nu, s, S, a, outgoing)*wu
            sT += v; sa += abs(v)
        IT += sT*wa; aT += sa*wa
        # I_L at A' = -A: v = V t^4
        alL = -A - 1j*eps
        sL = 0j; sa = 0.0
        for t, wt in zip(tq, wtq):
            vv = V*t**4; jac = 4*V*t**3
            if t == 0:
                continue
            v = W(alL, alL - 1j*vv, p, nu, s, S, a, outgoing)*jac*wt
            sL += v; sa += abs(v)
        IL += sL*wa; aL += sa*wa
    Cc = Cconst(n, N, s)
    comb = IU - 1j*IT - 1j*IL
    Kc = Cc/math.pi*comb.real if outgoing else Cc/(2*math.pi)*comb
    absI = (Cc/math.pi if outgoing else Cc/(2*math.pi))*(aU + aT + aL)
    return dict(IU=IU, IT=IT, IL=IL, K=Kc, absI=absI, C=Cc)


if __name__ == '__main__':
    n, N, b, p, nu = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
    nA = int(sys.argv[6]) if len(sys.argv) > 6 else 24; nt = int(sys.argv[7]) if len(sys.argv) > 7 else 24
    mp.mp.dps = 30
    ex = K.exact(p, nu, n, mp.mpf(N), 2*b + 1)
    r = pieces(n, N, b, p, nu, nA, nt)
    Kc = r['K']
    print(f'n={n} N={N} b={b} p={p} nu={nu}: I_U={r["IU"]:.8f} I_T={r["IT"]:.8f} I_L={r["IL"]:.8f}  contour={Kc}  '
          f'exact={mp.nstr(ex, 16)}  rel.err={abs(complex(Kc) - float(ex))/abs(float(ex)):.2e}  int|W|/|K|={r["absI"]/abs(float(ex)):.2e}',
          flush=True)
