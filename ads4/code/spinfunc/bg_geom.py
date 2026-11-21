"""
bg_geom (BELOW_FOLD2.md): the CLS-D crystal integrand (l1_crystal.build, scalar; w0_ee_common EE) re-built on the a-series algebra
of bg_ser, on the real section D (bf_core) in the variables p1, p2, q1, q2, s, phi, alpha (n0 = +1, psi = alpha - phi in (-pi, pi]).

Exact geometry (BELOW_FOLD §1.2):  x = a (p1^2, p2^2, s p1 p2),  y = a (q1^2, q2^2, q1 q2 sf),
  x~0 = cf cos(psi), y0 = a q1 q2 sin(psi), u0 = a q1^2 (cf cs + sf ss), v0 = a q2^2 (sf ss - cf cs),
  measure x meas-factor = a^5 p1 p2 dX1 dX2 ds dY1 dY2 dphi dpsi (the celestial Jacobian, exact).
assemble(): integrand = [const] * atoms^(E0 + i nu E1) * e^{-(X1+X2+Y1+Y2)} * sum_j (i nu/2)^j sum_k a^k Corr[j][k].
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from bg_ser import Poly, ASer, series_fn, lmul, NV, NAMES

M = Poly.mono


def vec(c, d, K):
    return [ASer.mono(d, p, K) if isinstance(p, Poly) else ASer.const(p, K) for p in c]


def dot(p, q):
    return (p[0]*q[1] + p[1]*q[0]).scale(-0.5) + p[2]*q[2]


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def neg(p):
    return [-c for c in p]


def scal(c, p):
    return [c*pi for pi in p]


def msq(p, q):
    r = sub(p, q)
    return -dot(r, r)


def geometry(K):
    """positions x0..x4, e, z0 as ASer 3-vectors (lightcone components (u, v, perp))."""
    cs = M(1, ca=1, cf=1) + M(1, sa=1, sf=1)                      # cos(psi), psi = alpha - phi
    ss = M(1, sa=1, cf=1) + M(-1, ca=1, sf=1)                     # sin(psi)
    x = [ASer.mono(1, M(1, p1=2), K), ASer.mono(1, M(1, p2=2), K), ASer.mono(1, M(1, s=1, p1=1, p2=1), K)]
    y = [ASer.mono(1, M(1, q1=2), K), ASer.mono(1, M(1, q2=2), K), ASer.mono(1, M(1, q1=1, q2=1, sf=1), K)]
    u0 = ASer.mono(1, M(1, q1=2)*(M(1, cf=1)*cs + M(1, sf=1)*ss), K)
    v0 = ASer.mono(1, M(1, q2=2)*(M(1, sf=1)*ss - M(1, cf=1)*cs), K)
    y0 = ASer.mono(1, M(1, q1=1, q2=1)*ss, K)
    x0 = [u0, v0, y0]
    xsq = dot(x, x)
    rx = (-xsq).recip('-x^2')
    x4 = scal(rx, x); x2 = neg(x4); x1 = neg(y); x3 = y
    e = [ASer.const(1.0, K), ASer.const(1.0, K), ASer.const(0.0, K)]
    z0 = [ASer.const(1.0, K), ASer.const(0.0, K), ASer.const(0.0, K)]
    return dict(x=x, y=y, x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, e=e, z0=z0, xsq=xsq)


def scalar_pieces(K):
    """the J = 0 CLS-D scalar integrand (l1_crystal.build, parts g) as lists:
    logs  [(name, f, p0, p1)]  meaning f^(p0 + (i nu/2) p1);   polys [(name, f)] multiplied (integer power 1);
    plus the geometry dict (for spin / B)."""
    G = geometry(K)
    x, y, x0, x1, x2, x3, x4, e, z0, xsq = (G[k] for k in ('x', 'y', 'x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0', 'xsq'))
    x10, x30, x20, x40 = sub(x1, x0), sub(x3, x0), sub(x2, x0), sub(x4, x0)
    m24, m20, m40 = msq(x2, x4), msq(x2, x0), msq(x4, x0)
    Nn = (dot(z0, x40)*(-m20) - dot(z0, x20)*(-m40)).scale(2.0)
    c13 = dot(x10, x30).scale(2.0)
    c10 = dot(z0, x10).scale(2.0); c30 = dot(z0, x30).scale(-2.0)
    m35, m45, m12, m15, m25, m34 = msq(x3, e), msq(x4, e), msq(x1, x2), msq(x1, e), msq(x2, e), msq(x3, x4)
    m13, m14, m23 = msq(x1, x3), msq(x1, x4), msq(x2, x3)
    xy = dot(x, y); x2y2 = xsq*dot(y, y)
    FP1 = 1.0 + xy.scale(2.0) + x2y2; FP2 = 1.0 - xy.scale(2.0) + x2y2; FP3 = 1.0 - x2y2
    up = m12*m34*(m13*m24).recip('m13 m24'); vp = m14*m23*(m13*m24).recip('m13 m24')
    upvp = up*vp
    fsub = (vp - up)*upvp.recip('u\'v\'')
    A_num = Nn*Nn*c13
    A_den = m24*m40*m20*c10*c30
    logs = [('A', A_num, 0.0, 1.0), ('A', A_den, 0.0, -1.0),
            ('xsq', -xsq, -6.0, 0.0), ('FP1', FP1, 1.0, 0.0), ('FP2', FP2, 1.0, 0.0), ('FP3', FP3, 1.0, 0.0),
            ('dual', Nn, -0.5, 0.0), ('dual', m24, 0.75, 0.0), ('dual', m20, -0.25, 0.0), ('dual', m40, -0.25, 0.0),
            ('cel', c13, 1.25, 0.0), ('cel', c10, -0.25, 0.0), ('cel', c30, -0.25, 0.0),
            ('3pt', m12, -3.0, 0.0), ('3pt', m34, -1.5, 0.0), ('3pt', m35, -1.5, 0.0), ('3pt', m45, -1.5, 0.0),
            ('meas', ASer.mono(5, M(1, p1=1, p2=1), K), 1.0, 0.0)]
    polys = [('sub', fsub)]
    named = dict(m12=m12, m34=m34, m13=m13, m24=m24, m14=m14, m23=m23, m15=m15, m25=m25, m35=m35, m45=m45, m20=m20, m40=m40,
                 Nn=Nn, c13=c13, c10=c10, c30=c30, up=up, vp=vp, FP1=FP1, FP2=FP2, FP3=FP3, A_num=A_num, A_den=A_den)
    return logs, polys, G, named


def bfac(K):
    """exp((Delta/2) log B + X1 + X2 + Y1 + Y2) at Delta = 1/a, B = B_x B_y, B_x = (1 - a sx + a^2 Qx)/(1 + a sx + a^2 Qx)
    (sx = X1 + X2, Qx = X1 X2 (1 - s^2); y: sy = Y1 + Y2, Qy = Y1 Y2 cf^2), as an ASer (d = 0)."""
    tot = ASer(0, [Poly() for _ in range(K + 1)])
    for (sv, Qv) in ((M(1, p1=2) + M(1, p2=2), M(1, p1=2, p2=2, w=1)), (M(1, q1=2) + M(1, q2=2), M(1, q1=2, q2=2, cf=2))):
        Kp = K + 1
        em = ASer(0, [Poly(), -sv, Qv] + [Poly() for _ in range(Kp - 2)])
        ep = ASer(0, [Poly(), sv, Qv] + [Poly() for _ in range(Kp - 2)])
        f = series_fn(em, 'log').c
        g = series_fn(ep, 'log').c
        h = [(fi - gi).scale(0.5) for fi, gi in zip(f, g)]       # (1/2a) f(a): coefficient of a^{k} is h[k+1]
        assert len(h[0]) == 0
        corr = [h[k + 1] for k in range(K + 1)]
        corr[0] = corr[0] + sv                                    # + (X1 + X2): remove e^{-S}
        tot = ASer(0, [t + c for t, c in zip(tot.c, corr)])
    assert np.max(np.abs(tot.c[0].ev(_pts()))) < 1e-12
    tot.c[0] = Poly()
    return series_fn(tot, 'exp')


def _pts():
    from bg_ser import _PTS
    return _PTS


def spin_t(G, K):
    """t = Z.W / sqrt(Z^2 W^2) (l1_crystal spin angle) -> (t0, delta ASer with d = 0, c[0] = 0)."""
    x1, x2, x3, x4, e = (G[k] for k in ('x1', 'x2', 'x3', 'x4', 'e'))
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    x15, x25, x53, x54 = sub(x1, e), sub(x2, e), sub(e, x3), sub(e, x4)
    Z = sub(scal(m25, x15), scal(m15, x25))
    W = sub(scal(m45, x53), scal(m35, x54))
    ZW = dot(Z, W); ZZWW = dot(Z, Z)*dot(W, W)
    c0, e0, eps = ZZWW.rel('Z^2W^2')
    assert np.all(e0 % 2 == 0), e0
    r = series_fn(eps, 'pow', -0.5)
    inv = ASer(-ZZWW.d//2, [p.mul_mono(-e0//2, c0**-0.5) for p in r.c])
    assert ZZWW.d % 2 == 0
    t = ZW*inv
    assert t.d == 0, t.d
    t0c, t0e = t.lead('t0')
    assert np.all(t0e == 0)
    dl = ASer(0, [Poly()] + [t.c[k] for k in range(1, t.K + 1)])
    return complex(t0c), dl


def assemble(logs, polys, K, nu_pow=2):
    """-> dict(E0, E1, lconst=(l0, l1), corr={j: [Poly]*K+1}) with integrand = exp(l0 + i nu l1) atoms^(E0 + i nu E1)
    e^{-S} sum_j (i nu/2)^j sum_k a^(dc + k) corr[j][k] (dc: overall a-power of the corrections, dropped)."""
    E0 = np.zeros(NV); E1 = np.zeros(NV); l0 = 0j; l1 = 0j
    L0 = [Poly() for _ in range(K + 1)]; L1 = [Poly() for _ in range(K + 1)]
    for (nm, f, p0, p1) in logs:
        c0, e0, eps = f.trunc(K).rel(nm) if f.K >= K else (None, None, None)
        if c0 is None:
            raise ValueError(f'{nm}: depth {f.K} < {K}')
        E0 += p0*e0; E1 += 0.5*p1*e0
        l0 += p0*np.log(c0); l1 += 0.5*p1*np.log(c0)
        lg = series_fn(eps, 'log').c
        L0 = [x + y.scale(p0) for x, y in zip(L0, lg)]
        L1 = [x + y.scale(p1) for x, y in zip(L1, lg)]
    base = series_fn(ASer(0, L0), 'exp').c
    for (nm, f) in polys:
        base = lmul(base, f.c[:K + 1], K)
    corr = {0: base}
    pw = [Poly.const(1.0)] + [Poly() for _ in range(K)]
    fac = 1.0
    for j in range(1, nu_pow + 1):
        pw = lmul(pw, L1, K); fac *= j
        if all(len(q) == 0 for q in pw):
            break
        corr[j] = [q.scale(1.0/fac) for q in lmul(pw, base, K)]
    return dict(E0=E0, E1=E1, lconst=(l0, l1), corr=corr)
