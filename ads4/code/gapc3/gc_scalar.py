"""
gc_scalar (GAP_C3.md §3): the l1_crystal integrand in LOG-FACTORED form, for Taylor-model certification on large polydiscs.

Same function as spinfunc/l1_crystal.build (scalar k = 2 functional, Q_J prefactor), rewritten so that every logarithm / fractional
power acts on an ELEMENTARY factor f_k (a polynomial or simple rational function of the 7 variables), relative to its value at the
crystal:
    Phi(z) - Phi(z*)  = (i xi/2) sum_k a_k Log(f_k/f_k*) + (1/2) sum_k b_k Log(f_k/f_k*)          (A- and B-factors)
    g(z) / g(z*)      = prod_k (f_k/f_k*)^{p_k}  x  rational factors x  Q_J(z)/Q_J(z*)                (amplitude e^G Q_J)
Each (f/f*)^p, Log(f/f*) is the analytic continuation from the crystal (the branch convention of l1_tps: principal log at z*, then
series); on a polydisc on which |f/f* - 1| < 1 it is the principal one of 1 + t.  This avoids the large intermediate products of the
original code (e.g. the denominator of A, ~ -775 at x = 0.5), whose Taylor-model ranges are badly overestimated.

Backends: 'tm' (gc_ball Taylor models, certified) and 'pt' (numpy complex points, for validation against n1_pw / l1_crystal).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import legendre as npleg


class PT:
    """point backend: arrays of complex values, index 0 = crystal."""
    @staticmethod
    def const(c, like):
        return np.full(like.shape, c, dtype=complex)

    @staticmethod
    def mul(a, b):
        return a*b

    @staticmethod
    def recip(a):
        return 1.0/a

    @staticmethod
    def rel(f):                          # t = f/f* - 1
        return f/f[..., :1] - 1.0

    @staticmethod
    def log1(t):
        return np.log1p(t)

    @staticmethod
    def pow1(t, p):
        return np.exp(p*np.log1p(t))

    @staticmethod
    def val0(f):
        return f[..., 0]

    @staticmethod
    def norm(f):
        return f/f[..., :1]


class TMB:
    """Taylor-model backend (gc_ball)."""
    def __init__(self, B):
        self.B = B

    def const(self, c, like=None):
        return self.B.const(c)

    def mul(self, a, b):
        return self.B.mul(a, b)

    def recip(self, a):
        return self.B.recip(a)

    def rel(self, f):
        B = self.B
        f0 = f[0]
        t = f/f0
        t[0] = t[0] - B.A(1)
        return t

    def _onepluslog(self, t, p=None):
        B = self.B
        one = B.const(1.0)
        if p is None:
            return B.log(one + t)        # continuation from the crystal: principal log(1) = 0 at the centre
        return B.powr(one + t, p)

    def log1(self, t):
        return self._onepluslog(t)

    def pow1(self, t, p):
        return self._onepluslog(t, p)

    def val0(self, f):
        return f[0]

    def norm(self, f):
        return f/f[0]


def legendre_poly_coeffs(J):
    c = npleg.leg2poly(np.eye(J + 1)[J])
    return [(p, c[p]) for p in range(J + 1) if abs(c[p]) > 0]


def build(zv, nu, Delta, J, K, y00, n0=1.0, xi=None):
    """returns dict: dPhi = Phi - Phi*, grel = g/g*, plus diagnostics.  y00: value of y0 at the crystal (from l1_crystal, principal)."""
    xt0, ux, vx, uy, vy, xp, yp = zv
    M = K.mul
    one = K.const(1.0, xt0)
    def dot(p, q):
        return -0.5*(M(p[0], q[1]) + M(p[1], q[0])) + M(p[2], q[2])
    def sub(p, q):
        return [p[i] - q[i] for i in range(3)]
    def msq(p, q):
        r = sub(p, q); return -dot(r, r)
    def scal(c, p):
        return [M(c, p[i]) for i in range(3)]
    uv = M(uy, vy); yp2 = M(yp, yp)
    F2 = M(uv, one - M(xt0, xt0)) - yp2
    F3 = uv - yp2
    tuy, tvy, tF2, tF3 = K.rel(uy), K.rel(vy), K.rel(F2), K.rel(F3)
    y0 = y00*M(M(K.pow1(tuy, 0.5), K.pow1(tvy, 0.5)), M(K.pow1(tF2, 0.5), K.pow1(tF3, -0.5)))
    u0 = M(uy, xt0) + M(M(y0, yp), K.recip(vy))
    v0 = -M(vy, xt0) + M(M(y0, yp), K.recip(uy))
    x0 = [u0, v0, y0]
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = dot(x, x)
    rx = K.recip(-xsq)
    x4 = scal(rx, x); x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    e = [K.const(1.0, xt0), K.const(1.0, xt0), K.const(0.0, xt0)]
    z0 = [K.const(1.0, xt0), K.const(0.0, xt0), K.const(0.0, xt0)]
    x10, x30, x20, x40 = sub(x1, x0), sub(x3, x0), sub(x2, x0), sub(x4, x0)
    m24, m20, m40 = msq(x2, x4), msq(x2, x0), msq(x4, x0)
    Nn = 2*(M(dot(z0, x40), -m20) - M(dot(z0, x20), -m40))
    c13 = 2*dot(x10, x30)
    c10 = 2*dot(z0, x10); c30 = -2*dot(z0, x30)
    m35, m45, m12, m15, m25, m34 = msq(x3, e), msq(x4, e), msq(x1, x2), msq(x1, e), msq(x2, e), msq(x3, x4)
    m13, m14, m23 = msq(x1, x3), msq(x1, x4), msq(x2, x3)
    xy = dot(x, y); x2y2 = M(xsq, dot(y, y))
    fac = dict(uy=uy, vy=vy, F2=F2, F3=F3, mxsq=-xsq, Nn=Nn, c13=c13, c10=c10, c30=c30, m24=m24, m20=m20, m40=m40,
               m35=m35, m45=m45, m12=m12, m15=m15, m25=m25, m34=m34, m13=m13, m14=m14, m23=m23,
               fp1=one + 2*xy + x2y2, fp2=one - 2*xy + x2y2, fp3=one - x2y2)
    t = {k: K.rel(v) for k, v in fac.items()}
    L = {}
    def lg(k):
        if k not in L:
            L[k] = K.log1(t[k])
        return L[k]
    # phase: log A = 2 log Nn + log c13 - log m24 - log m40 - log m20 - log c10 - log c30;  log B = m35 m45 m12 / (m15 m25 m34)
    dlA = 2*lg('Nn') + lg('c13') - lg('m24') - lg('m40') - lg('m20') - lg('c10') - lg('c30')
    dlB = lg('m35') + lg('m45') + lg('m12') - lg('m15') - lg('m25') - lg('m34')
    if xi is None:
        xi = nu/Delta
        dPhi = (0.5j*xi)*dlA + 0.5*dlB
    else:                                 # xi: acb ball (an interval of xi values, Taylor models then cover all of them)
        from flint import acb
        dPhi = (acb(0, 0.5)*xi)*dlA + 0.5*dlB
    # amplitude e^G / e^G*: powers of elementary factors (rational pieces exactly), fsub as a rational function
    P = {'uy': 0.5, 'vy': 0.5, 'F3': -0.5, 'F2': -0.5, 'mxsq': -6.0, 'fp1': 1.0, 'fp2': 1.0, 'fp3': 1.0,
         'Nn': -0.5, 'm24': 0.75, 'm20': -0.25, 'm40': -0.25, 'c13': 1.25, 'c10': -0.25, 'c30': -0.25,
         'm12': -(6 + J)/2, 'm15': -J/2, 'm25': -J/2, 'm34': -(3 + J)/2, 'm35': -(3 + J)/2, 'm45': -(3 + J)/2}
    grel = None
    for k, p in P.items():
        if p == 0:
            continue
        f = K.pow1(t[k], p)
        grel = f if grel is None else M(grel, f)
    # f_sub = (vp - up)/(up vp) = m13 m24 (m14 m23 - m12 m34)/(m12 m34 m14 m23)
    num = M(m14, m23) - M(m12, m34)
    fsub = M(M(M(m13, m24), num), K.recip(M(M(m12, m34), M(m14, m23))))
    grel = M(grel, K.norm(fsub))
    Q = None; ZZWW = None
    if J > 0:
        x15, x25, x53, x54 = sub(x1, e), sub(x2, e), sub(e, x3), sub(e, x4)
        Z = sub(scal(m25, x15), scal(m15, x25))
        W = sub(scal(m45, x53), scal(m35, x54))
        ZW = dot(Z, W); ZZWW = M(dot(Z, Z), dot(W, W))
        for (p, c) in legendre_poly_coeffs(J):
            term = K.const(c, xt0)
            for _ in range(p):
                term = M(term, ZW)
            for _ in range((J - p)//2):
                term = M(term, ZZWW)
            Q = term if Q is None else Q + term
    return dict(dPhi=dPhi, grel=grel, Q=Q, ZZWW=ZZWW, fac=fac, t=t)
