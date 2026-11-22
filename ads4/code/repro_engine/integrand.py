"""repro_engine/integrand.py -- point-wise EXACT (all orders in 1/Delta) EE k = 4 CLS-D integrand on the real section D, J = 0,
own numpy implementation (no project module imported).

Definition (read from spinfunc/bg_geom.py + bg_ee.py, which re-build l1_crystal.build + w0_ee_common.make_build_even_EE):
 variables X1, X2, Y1, Y2 > 0, s in (-1,1), phi in (-pi/2, pi/2), alpha in (-pi, pi]; a = 1/Delta; psi = alpha - phi; lightcone
 3-vectors (u, v, perp), x.y = -(u v' + v u')/2 + x_perp y_perp;
   x = a(X1, X2, s sqrt(X1 X2)),  y = a(Y1, Y2, sqrt(Y1 Y2) sin phi),
   x0 = a(Y1 cos(psi - phi), -Y2 cos(alpha), sqrt(Y1 Y2) sin psi),  x4 = x/(-x^2), x2 = -x4, x1 = -y, x3 = y, e = (1,1,0), z0 = (1,0,0);
 m_ij = -(x_i - x_j)^2;  Nn = 2(z0.x40 (-m20) - z0.x20 (-m40)), c13 = 2 x10.x30, c10 = 2 z0.x10, c30 = -2 z0.x30;
 FP1,2 = 1 +- 2 x.y + x^2 y^2, FP3 = 1 - x^2 y^2;  u' = m12 m34/(m13 m24), v' = m14 m23/(m13 m24);
 A = Nn^2 c13/(m24 m40 m20 c10 c30);
 integrand = A^{i nu/2} (-x^2)^-6 FP1 FP2 FP3 Nn^-1/2 m24^3/4 m20^-1/4 m40^-1/4 c13^5/4 c10^-1/4 c30^-1/4 sqrt(X1 X2)
             m12^-5 m34^-7/2 m35^-3/2 m45^-3/2 V042^-4 (v'-u') pref R B^{Delta/2},
   R = ((v'-u')^2 - 1)/(4 u'v')  (Lemma B, k = 4),  B = Bx By, Bx = (1 - a(X1+X2) + a^2 X1X2(1-s^2))/(1 + a(X1+X2) + a^2 X1X2(1-s^2)),
   By likewise with Y1Y2 cos^2 phi;  measure dX1 dX2 dY1 dY2 ds dphi dalpha  (overall a-powers and constant phases dropped: they cancel
   in the state/reference ratio through the nu = 0.02 normalisation).
 pref = the EE polarisation prefactor:  <T(x1, z=x10) T(x2, .) O(e)>_{Delta} (x)_{x2} E(.,.)^2 (x)_{x4} <T(x3, z=x30) T(x4, .) O(e)>_{3-Delta},
   E(z2, z4) = V042^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)  (X_ij = (x_i - x_j)^2), V, H as in CPPR (see tto_conserved.py),
   <TTO> numerators = sum_c coef_c(Delta) {V1^2 V2^2, V1 V2 H12, H12^2} with coef from tto_conserved.py (own derivation).
 Contraction at x2, x4 (spin 2, d = 3): OWN method -- the quadratic forms are made traceless (M -> M - tr(g^-1 M) g/3) and contracted
 with the inverse metric, T^{mn} S_{mn}  (normalisation <(a.z)^2, (b.z)^2> = (a.b)^2; the project uses a null-circle Fourier
 contraction instead).

Non-integer powers are taken of |f|: every such factor has a fixed sign on the sampled region (checked by sign_report), so the
difference is a constant phase common to state and reference.
"""
import numpy as np

GL = np.array([[0.0, -0.5, 0.0], [-0.5, 0.0, 0.0], [0.0, 0.0, 1.0]])      # metric (lightcone comps)
GI = np.linalg.inv(GL)


def coef(Delta):
    """conserved <TTO_Delta> coefficients (V1^2V2^2, V1V2H12, H12^2), J = 0 (tto_conserved.py)."""
    d = float(Delta)
    den = d*d - 6*d + 6
    return np.array([d*(d + 2)/(2*den), 2*d*(d - 4)/den, 1.0])


def dot(p, q):
    return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def sq(p, q):
    r = sub(p, q); return dot(r, r)


def Vf(xi, z, xj, xk):
    """V_{i,jk}(z) = (z.(xj-xi) (xi-xk)^2 - z.(xk-xi) (xi-xj)^2)/(xj-xk)^2."""
    return (dot(z, sub(xj, xi))*sq(xi, xk) - dot(z, sub(xk, xi))*sq(xi, xj))/sq(xj, xk)


def Hf(xi, zi, xj, zj):
    """H_ij = (zi.zj)(xi-xj)^2 + 2 zi.(xj-xi) zj.(xi-xj)."""
    return dot(zi, zj)*sq(xi, xj) + 2*dot(zi, sub(xj, xi))*dot(zj, sub(xi, xj))


BAS = [np.eye(3)[m] for m in range(3)]


def cov_lin(fun):
    """covariant components f_m = f(B_m) of a function linear in z (z = z^m B_m)."""
    return np.stack([fun(BAS[m]) for m in range(3)], axis=0)            # (3, n)


def traceless_raise(Mc):
    """Mc (3,3,n) covariant symmetric -> contravariant traceless T^{mn}."""
    tr = np.einsum('ab,abn->n', GI, Mc)
    Mh = Mc - tr[None, None, :]*GL[:, :, None]/3.0
    return np.einsum('am,bn,mnk->abk', GI, GI, Mh)


def tto_form(xa, za, xb, x5, cf):
    """<T(xa, za) T(xb, z) O(x5)> numerator as a covariant quadratic form in z (3,3,n)."""
    V1 = Vf(xa, za, xb, x5)
    v = cov_lin(lambda B: Vf(xb, [B[0], B[1], B[2]], x5, xa))
    h = cov_lin(lambda B: Hf(xa, za, xb, [B[0], B[1], B[2]]))
    vv = np.einsum('mn,kn->mkn', v, v); vh = 0.5*(np.einsum('mn,kn->mkn', v, h) + np.einsum('mn,kn->mkn', h, v))
    hh = np.einsum('mn,kn->mkn', h, h)
    return cf[0]*V1**2*vv + cf[1]*V1*vh + cf[2]*hh


def point(X1, X2, Y1, Y2, s, phi, al, Delta, want=('G', 'lA'), Dcoef=None):
    """exact integrand at points (arrays, same shape). Returns dict:
       lG  = log|integrand without A^{i nu/2}| (real),  sG = its sign,  lA = log A (A > 0 checked),  plus diagnostics."""
    a = 1.0/Delta
    p1, p2, q1, q2 = np.sqrt(X1), np.sqrt(X2), np.sqrt(Y1), np.sqrt(Y2)
    cf, sf, ca, sa = np.cos(phi), np.sin(phi), np.cos(al), np.sin(al)
    cs = ca*cf + sa*sf; ss = sa*cf - ca*sf
    x = [a*X1, a*X2, a*s*p1*p2]
    y = [a*Y1, a*Y2, a*q1*q2*sf]
    x0 = [a*Y1*(cf*cs + sf*ss), a*Y2*(sf*ss - cf*cs), a*q1*q2*ss]
    xsq = -a*a*X1*X2*((1.0 - s)*(1.0 + s))           # = dot(x, x), without the cancellation at s -> +-1
    x4 = [c/(-xsq) for c in x]; x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    one = np.ones_like(X1); zer = np.zeros_like(X1)
    e = [one, one, zer]; z0 = [one, zer, zer]
    msq = lambda p, q: -sq(p, q)
    x10, x30, x20, x40 = sub(x1, x0), sub(x3, x0), sub(x2, x0), sub(x4, x0)
    m24, m20, m40 = msq(x2, x4), msq(x2, x0), msq(x4, x0)
    Nn = 2.0*(dot(z0, x40)*(-m20) - dot(z0, x20)*(-m40))
    c13 = 2.0*dot(x10, x30); c10 = 2.0*dot(z0, x10); c30 = -2.0*dot(z0, x30)
    m35, m45, m12, m34 = msq(x3, e), msq(x4, e), msq(x1, x2), msq(x3, x4)
    m13, m14, m23 = msq(x1, x3), msq(x1, x4), msq(x2, x3)
    xy = dot(x, y); x2y2 = xsq*(-a*a*Y1*Y2*cf*cf)       # y^2 = -a^2 Y1 Y2 cos^2 phi
    FP1 = 1.0 + 2*xy + x2y2; FP2 = 1.0 - 2*xy + x2y2; FP3 = 1.0 - x2y2
    up = m12*m34/(m13*m24); vp = m14*m23/(m13*m24)
    # v' - u' without the 1/a^2 cancellation: x3 = -x1, x2 = -x4 => m34 = m12, m23 = m14, m14 - m12 = -4 x4.y, m14 + m12 = -2(x4^2 + y^2)
    vmu = 8.0*dot(x4, y)*(dot(x4, x4) + dot(y, y))/(m13*m24)
    R = (vmu*vmu - 1.0)/(4*up*vp)
    v042 = Vf(x0, z0, x4, x2)
    # EE prefactor
    Dc = Delta if Dcoef is None else Dcoef            # Dcoef: <TTO> coefficients at another Delta (order-0 diagnostics only)
    F = tto_form(x1, x10, x2, e, coef(Dc))
    Gt = tto_form(x3, x30, x4, e, coef(3.0 - Dc))
    Fu = traceless_raise(F); Gu = traceless_raise(Gt)
    inv = 1.0/(sq(x2, x0)*sq(x4, x0))
    H24 = np.stack([np.stack([Hf(x2, list(BAS[m]), x4, list(BAS[k])) for k in range(3)]) for m in range(3)])   # (3,3,n)
    V204 = cov_lin(lambda B: Vf(x2, [B[0], B[1], B[2]], x0, x4))
    V402 = cov_lin(lambda B: Vf(x4, [B[0], B[1], B[2]], x0, x2))
    E = v042**2*inv*(H24 - 2*np.einsum('mn,kn->mkn', V204, V402))                                            # E_{mk}
    pref = np.einsum('mnz,mkz,nlz,klz->z', Fu, E, E, Gu)
    # B^{Delta/2}
    sx, Qx = X1 + X2, X1*X2*(1 - s*s); sy, Qy = Y1 + Y2, Y1*Y2*cf*cf
    lB = (np.log1p(-a*sx + a*a*Qx) - np.log1p(a*sx + a*a*Qx) + np.log1p(-a*sy + a*a*Qy) - np.log1p(a*sy + a*a*Qy))/(2*a)
    # assemble log|.| and sign
    facs = [(-xsq, -6.0), (FP1, 1.0), (FP2, 1.0), (FP3, 1.0), (Nn, -0.5), (m24, 0.75), (m20, -0.25), (m40, -0.25),
            (c13, 1.25), (c10, -0.25), (c30, -0.25), (p1*p2, 1.0), (m12, -5.0), (m34, -3.5), (m35, -1.5), (m45, -1.5),
            (v042, -4.0), (vmu, 1.0), (pref, 1.0), (R, 1.0)]
    lG = lB.copy(); sG = np.ones_like(lB); signs = {}
    names = ['-xsq', 'FP1', 'FP2', 'FP3', 'Nn', 'm24', 'm20', 'm40', 'c13', 'c10', 'c30', 'meas', 'm12', 'm34', 'm35', 'm45',
             'v042', 'vmu', 'pref', 'R']
    for nm, (f, pw) in zip(names, facs):
        lG = lG + pw*np.log(np.abs(f))
        sgn = np.sign(f)
        signs[nm] = sgn
        if float(pw).is_integer():
            sG = sG*sgn**int(pw)
    Anum = Nn*Nn*c13; Aden = m24*m40*m20*c10*c30
    lA = np.log(np.abs(Anum)) - np.log(np.abs(Aden))
    signs['A'] = np.sign(Anum*Aden)
    return dict(lG=lG, sG=sG, lA=lA, signs=signs, pref=pref, R=R, vmu=vmu, lB=lB)


if __name__ == '__main__':
    rng = np.random.default_rng(0)
    n = 6
    X1, X2, Y1, Y2 = (rng.gamma(k, size=n) for k in (4, 3, 4, 3.5))
    s = rng.uniform(-0.9, 0.9, n); phi = rng.uniform(-1.4, 1.4, n); al = rng.uniform(-3, 3, n)
    for Dl in (300.0, 2e4):
        P = point(X1, X2, Y1, Y2, s, phi, al, Dl)
        print(Dl, P['lG'], P['sG'], P['lA'])
        print({k: v for k, v in P['signs'].items()})
