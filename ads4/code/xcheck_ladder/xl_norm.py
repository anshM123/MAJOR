"""
xcheck_ladder/xl_norm.py: independent (mpmath) implementation of the per-structure flat normalisation N_s and the hybrid action, as
defined in code/spinfunc/ee1000_common.hyb / code/largespin/ls_core.hyb_fast (the convention of every direct row, UT_v2 Sec. 2.5):

   y = u p^2, u = 1/mfl^2 (mfl = m/Lam),   Dt_++(y) = P_J(1-2y),   Dt_+-(y) = P^(0,8)_{J-4}(1-2y)
   B2 = int_0^1 psi2 (2-y) Dt dp,          B3 = sigma u int_0^1 psi3 Dt dp                       (flat NF02 parts)
   c2 = 2 [OmLR_2(beta) - OmFlat_2(bhat)],  c3 = sigma u [OmLR_3(beta) - OmFlat_3(bhat)]          (leading-Regge curvature)
   OmLR_k(beta) = int psi_k(p) P_{-1/2+i Lam p}(cosh beta) dp,  OmFlat_k(bhat) = int psi_k J0(bhat p) dp,
   sinh(beta/2) = (J+1/2)/m,  bhat = (2J+1)/mfl
   N_s = |B2 + c2| + |B3 + c3|,   H_hyb = B2 + B3 + c2 + c3.
psi_k = NF02 (prl/data/nf02_psi.json, exact rationals / k0).  P_{-1/2+i nu}(cosh beta) = 2F1(1/2-i nu, 1/2+i nu; 1; -sinh^2(beta/2)).
"""
import sys, os
sys.dont_write_bytecode = True
import mpmath as mp
import xl_common as X

_P = None
_GL = {}


def psis():
    global _P
    if _P is None:
        p2, p3 = X.nf02_exact()
        _P = ([mp.mpf(c.numerator)/c.denominator for c in p2], [mp.mpf(c.numerator)/c.denominator for c in p3])
    return _P


def poly(c, p):
    return mp.polyval(c[::-1], p)


def gl01(n):
    if n not in _GL:
        xs, ws = mp.gauss_quadrature(n, 'legendre') if hasattr(mp, 'gauss_quadrature') else _gl(n)
        _GL[n] = [((x + 1)/2, w/2) for x, w in zip(xs, ws)]
    return _GL[n]


def _gl(n):
    # Golub-Welsch-free: Newton on Legendre roots
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi*(i - mp.mpf(1)/4)/(n + mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2*k - 1)*x*p1 - (k - 1)*p0)/k
            dp = n*(x*p1 - p0)/(x*x - 1)
            dx = p1/dp
            x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps - 3):
                break
        xs.append(x); ws.append(2/((1 - x*x)*dp*dp))
    return xs, ws


def jacobi08(n, x):
    """P^(0,8)_n(x) by the three-term recurrence (a = 0, b = 8)."""
    a, b = 0, 8
    if n == 0:
        return mp.mpf(1)
    p0 = mp.mpf(1); p1 = (a - b + (a + b + 2)*x)/2
    for k in range(2, n + 1):
        c1 = 2*k*(k + a + b)*(2*k + a + b - 2)
        c2 = (2*k + a + b - 1)*((2*k + a + b)*(2*k + a + b - 2)*x + a*a - b*b)
        c3 = 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
        p0, p1 = p1, (c2*p1 - c3*p0)/c1
    return p1


def Pcon(nu, sh2):
    """conical function P_{-1/2+i nu}(cosh beta), sh2 = sinh^2(beta/2)."""
    return mp.re(mp.hyp2f1(mp.mpf(1)/2 - 1j*nu, mp.mpf(1)/2 + 1j*nu, 1, -sh2))


def hybrid(Dg, J, mfl, h, npn=240, dps=30):
    """returns dict(B2, B3, c2, c3, Ns, Hhyb) for structure h in {'++', '+-'}; mfl may be a string."""
    with mp.workdps(dps):
        Lam = mp.sqrt(X.lam_int(Dg)); mfl = mp.mpf(mfl)
        u = 1/mfl**2; s = 1 if h == '++' else -1
        p2, p3 = psis()
        m = mfl*Lam
        sh2 = ((J + mp.mpf(1)/2)/m)**2
        bh = (2*J + 1)/mfl
        B2 = B3 = L2 = L3 = F2 = F3 = mp.mpf(0)
        for p, w in gl01(npn):
            y = u*p*p
            Dt = mp.legendre(J, 1 - 2*y) if h == '++' else jacobi08(J - 4, 1 - 2*y)
            a2 = poly(p2, p); a3 = poly(p3, p)
            B2 += w*a2*(2 - y)*Dt; B3 += w*a3*Dt
            P = Pcon(Lam*p, sh2); j0 = mp.besselj(0, bh*p)
            L2 += w*a2*P; L3 += w*a3*P; F2 += w*a2*j0; F3 += w*a3*j0
        B3 *= s*u
        c2 = 2*(L2 - F2); c3 = s*u*(L3 - F3)
        Ns = abs(B2 + c2) + abs(B3 + c3)
        return dict(B2=B2, B3=B3, c2=c2, c3=c3, Ns=Ns, Hhyb=B2 + B3 + c2 + c3, OmLR2=L2, OmLR3=L3)
