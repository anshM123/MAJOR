"""
l1 (M6 step 1): CLS-D scalar functional Psi_{2,nu} (d = 3, Delta_phi = 3, subtraction f = (v'-u')/(u'v')) on a block G_{Delta,J}:
the exact gauge-fixed 7-dim integrand (CLS-D eq. scalarfunctional_plus_gaugefixed_nozintegral, delta functions solved, d = 3 so
S^{d-3} = {n0 = +-1}), as TPS around the scattering crystal, and its Laplace expansion to relative O(1/Delta^2).

Variables z = (xt0, ux, vx, uy, vy, xp, yp)  [CLS-D (t0x0tilde_definition), t0~ = 0 and |y0| from the delta functions].
Gauge: x1 = -y, x3 = y, x4 = -x2 = x/(-x^2), x5 = e = (1,1,0), z0 = (1,0,0); lightcone x^2 = -u v + xp^2.
Fast exponent (CLS-D bulkpoint_integrand_factor with m -> Delta, exact):  lam F = (i nu/2) log A + (Delta/2) log B.
Slow G = log of everything else that depends on z (measure, (-x^2)^{-2d}, FP determinant, dual / celestial / 3-pt remainders,
subtraction factor).  z-independent factors: A_{k,nu}, sin^2, vol's (nu-only or state constants), 1/beta_{Delta,J} (exact, added).
The spin-J contraction of the two 3-pt tensor structures is Q_J = (Z^2 W^2)^{J/2} P_J(Z.W/sqrt(Z^2 W^2)) (d = 3), up to a J-only
constant.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import legendre as npleg
from scipy.special import loggamma
import l1_tps as T

D = 3


def dot(p, q):
    return -0.5*(T.mul(p[0], q[1]) + T.mul(p[1], q[0])) + T.mul(p[2], q[2])


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def neg(p):
    return [-p[i] for i in range(3)]


def scal(c, p):
    return [T.mul(c, p[i]) if isinstance(c, np.ndarray) else c*p[i] for i in range(3)]


def cpt(u, v, t):
    return [T.const(u), T.const(v), T.const(t)]


def legendre_poly_coeffs(J):
    """P_J(t) = sum_i c_i t^{J-2i}: return list of (power, coef)."""
    c = npleg.leg2poly(np.eye(J + 1)[J])
    return [(p, c[p]) for p in range(J + 1) if abs(c[p]) > 0]


def build(zv, nu, Delta, J, n0=1.0, parts=False, pref=False):
    """zv: list of 7 TPS. Returns (lamF, G) TPS (and dict of pieces if parts)."""
    xt0, ux, vx, uy, vy, xp, yp = zv
    x = [ux, vx, xp]; y = [uy, vy, yp]
    uv = T.mul(uy, vy)
    yp2 = T.mul(yp, yp)
    one = T.const(1.0)
    y0sq = T.div(T.mul(uv, T.mul(uv, one - T.mul(xt0, xt0)) - yp2), uv - yp2)
    y0 = n0*T.powr(y0sq, 0.5)
    u0 = T.mul(uy, xt0) + T.div(T.mul(y0, yp), vy)
    v0 = -T.mul(vy, xt0) + T.div(T.mul(y0, yp), uy)
    x0 = [u0, v0, y0]
    xsq = dot(x, x)
    rx = T.recip(-xsq)
    x4 = scal(rx, x); x2 = neg(x4); x1 = neg(y); x3 = y
    e = cpt(1.0, 1.0, 0.0)
    z0 = cpt(1.0, 0.0, 0.0)
    def msq(p, q):                       # -(p - q)^2
        r = sub(p, q); return -dot(r, r)
    x10 = sub(x1, x0); x30 = sub(x3, x0); x20 = sub(x2, x0); x40 = sub(x4, x0)
    m24, m20, m40 = msq(x2, x4), msq(x2, x0), msq(x4, x0)
    Nn = 2*(T.mul(dot(z0, x40), -m20) - T.mul(dot(z0, x20), -m40))     # 2 z0.x40 x20^2 - 2 z0.x20 x40^2
    N2 = T.mul(Nn, Nn)
    c13 = 2*dot(x10, x30)
    c10 = 2*dot(z0, x10); c30 = -2*dot(z0, x30)
    A = T.div(T.mul(N2, c13), T.mul(T.mul(m24, m40), T.mul(m20, T.mul(c10, c30))))
    m35, m45, m12, m15, m25, m34 = msq(x3, e), msq(x4, e), msq(x1, x2), msq(x1, e), msq(x2, e), msq(x3, x4)
    B = T.div(T.mul(T.mul(m35, m45), m12), T.mul(T.mul(m15, m25), m34))
    lA, lB = T.log(A), T.log(B)
    lamF = (0.5j*nu)*lA + (0.5*Delta)*lB
    # ---- slow pieces (logs)
    g = {}
    g['meas'] = 0.5*T.log(uv) - 0.5*T.log(uv - yp2) - 0.5*T.log(T.mul(uv, one - T.mul(xt0, xt0)) - yp2)
    g['xsq'] = -2*D*T.log(-xsq)
    xy = dot(x, y); x2y2 = T.mul(xsq, dot(y, y))
    g['FP'] = T.log(one + 2*xy + x2y2) + T.log(one - 2*xy + x2y2) + (D - 2)*T.log(one - x2y2)
    g['dual'] = -0.25*T.log(N2) + 0.75*T.log(m24) - 0.25*T.log(m20) - 0.25*T.log(m40)
    g['cel'] = 0.25*T.log(T.div(c13, T.mul(c10, c30))) + T.log(c13)
    m13, m14, m23 = msq(x1, x3), msq(x1, x4), msq(x2, x3)
    up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
    fsub = T.div(vp - up, T.mul(up, vp))
    g['sub'] = T.log(fsub)
    g['3pt'] = (-(6 + J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25)) - ((3 + J)/2)*(T.log(m34) + T.log(m35) + T.log(m45))
    if J > 0:
        x15, x25, x53, x54 = sub(x1, e), sub(x2, e), sub(e, x3), sub(e, x4)
        Z = sub(scal(m25, x15), scal(m15, x25))          # x15 |x25^2| - x25 |x15^2|   (overall sign (-1)^J irrelevant)
        W = sub(scal(m45, x53), scal(m35, x54))
        ZW = dot(Z, W); ZZWW = T.mul(dot(Z, Z), dot(W, W))
        Q = T.const(0.0)
        for (p, c) in legendre_poly_coeffs(J):
            term = T.const(c)
            for _ in range(p):
                term = T.mul(term, ZW)
            for _ in range((J - p)//2):
                term = T.mul(term, ZZWW)
            Q = Q + term
        if not pref:
            g['QJ'] = T.log(Q)
    G = sum(g.values())
    if pref:                                  # Q_J as polynomial prefactor (uniform near zeros of P_J(1-2x^2))
        return lamF, G, (Q if J > 0 else T.const(1.0))
    if parts:
        return lamF, G, dict(g=g, lA=lA, lB=lB, fsub=fsub, y0=y0, Q=(Q if J > 0 else None))
    return lamF, G


def saddle(nu, Delta):
    a = (Delta - np.sqrt(Delta**2 - nu**2))/nu
    return np.array([0.0, -1j*a, 1j*a, 1j*a, -1j*a, 0.0, 0.0], complex), a


def log_beta(Delta, J, d=D):
    """log |beta_{Delta,J}| (CLS-D shadow_coeff_scalar)."""
    return ((J - 1)*np.log(2) + (d - 2)/2*np.log(np.pi) + loggamma(Delta + 2 - d) + 2*loggamma((J + Delta)/2)
            + 2*loggamma((Delta + 2 - d - J)/2) - loggamma(Delta + (2 - d)/2) - loggamma(J + Delta) - loggamma(Delta + 2 - d - J)).real


def laplace(nu, Delta, J, n0=1.0, check=False):
    """Laplace expansion around the crystal.  Returns dict: logY0 (complex, incl. -log beta), K1, K2 (relative corrections)."""
    zs, a = saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF, G = build(zv, nu, Delta, J, n0)
    gF, HF = T.grad_hess(lamF)
    R = T.csym_cholesky(-HF)                  # -H = R R^T  ->  L = R^{-T}:  L^T (-H) L = 1
    L = np.linalg.inv(R).T
    for _ in range(2):                        # refine (round-off) on the xi-pass Hessian
        xv = T.linsub(zs, L)
        lamF2, G2 = build(xv, nu, Delta, J, n0)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    h = [T.homog(lamF2, k) for k in range(T.DEG + 1)]
    gg = [T.homog(G2, k) for k in range(T.DEG + 1)]
    A1 = h[3] + gg[1]; A2 = h[4] + gg[2]; A3 = h[5] + gg[3]; A4 = h[6] + gg[4]
    A1s = T.mul(A1, A1)
    K1 = T.expect(A2) + 0.5*T.expect(A1s)
    K2 = (T.expect(A4) + 0.5*T.expect2(A2, A2) + T.expect2(A1, A3) + 0.5*T.expect2(A1s, A2) + T.expect2(A1s, A1s)/24.0)
    logdetL = np.log(np.linalg.det(L))
    logY0 = lamF[0] + G[0] + 3.5*np.log(2*np.pi) + logdetL - log_beta(Delta, J)
    out = dict(logY0=logY0, K1=K1, K2=K2, a=a, gradF=gF, HF=HF, lamF0=lamF[0], G0=G[0],
               quad_check=np.max(np.abs(T.homog(lamF2, 2)[T.DEGS == 2] - np.where(
                   np.array([max(e) == 2 for e in T.EXP[T.DEGS == 2]]), -0.5, 0.0))))
    return out


def laplace_pref(nu, Delta, J, n0=1.0):
    """As laplace(), with the spin contraction Q_J kept as a prefactor (l1_tps.expand_pref).  Returns logY0 (without Q),
    and S = (S0, S2, S4): Y = exp(logY0) (S0 + S2 + S4)."""
    zs, a = saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF, G, Q = build(zv, nu, Delta, J, n0, pref=True)
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2, G2, Q2 = build(T.linsub(zs, L), nu, Delta, J, n0, pref=True)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    lamF2, G2, Q2 = build(T.linsub(zs, L), nu, Delta, J, n0, pref=True)
    S = T.expand_pref(lamF2, G2, Q2)
    logY0 = lamF[0] + G[0] + 3.5*np.log(2*np.pi) + np.log(np.linalg.det(L)) - log_beta(Delta, J)
    return dict(logY0=logY0, S=S)
