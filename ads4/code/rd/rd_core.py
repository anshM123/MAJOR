"""
rd_core (RD_CLOSE.md): vectorised, branch-free evaluation of the CLS-D crystal integrand (J = 0, sheet n0 = +1) at COMPLEX raw points,
with exact holomorphic first derivatives (forward-mode dual numbers), for

  * the scalar flow function of the n3/n4 crystal thimbles (n1_direct.phi_w, re-implemented here):
        Phi_w = lamF + G_scal + logJac(w),   lamF = (i nu/2) log A + (Delta/2) log B,
    written as  Phi = sum_k e_k log f_k + logJac  with every f_k a RATIONAL function of the raw point (and of y0 = sqrt(y0sq),
    continued by proximity), so that d Phi = sum_k e_k df_k/f_k is single valued: no log replay is needed for gradient flows;
  * the EE (k = 4) / scalar ratio  Q = Phi_EE / Phi_scal  (a single-valued rational function, RD_CLOSE section 1):
        Q = m12^-2 m34^-2 V042^-4 pref ((v'-u')^2 - 1)/4,
    pref the EE polarisation contraction (repro_engine/integrand.py functions, imported read-only).

Raw variables z = (xt0, ux, vx, uy, vy, xp, yp) as in spinfunc/l1_crystal.py; w-coordinates as in spinfunc/n1_direct.zmap
(u_x, v_x, u_y, v_y = c_i e^{w_i}, c = crystal; xt0 = tanh w0; others linear).
No project file is modified; spinfunc/l1_crystal (saddle) and repro_engine/integrand (EE contraction) are imported read-only.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
ADS = os.path.dirname(HERE)
import numpy as np

NV = 7


# ----------------------------------------------------------------------------------------------- dual numbers (value, d/dz_0..6)
class Dual:
    __slots__ = ('v', 'd')
    __array_priority__ = 100

    def __init__(self, v, d):
        self.v = v; self.d = d

    @staticmethod
    def _c(o):
        return o if isinstance(o, Dual) else None

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v + o.v, self.d + o.d)
        return Dual(self.v + o, self.d)
    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v - o.v, self.d - o.d)
        return Dual(self.v - o, self.d)

    def __rsub__(self, o):
        return Dual(o - self.v, -self.d)

    def __neg__(self):
        return Dual(-self.v, -self.d)

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v*o.v, self.d*o.v + o.d*self.v)
        return Dual(self.v*o, self.d*o)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, Dual):
            iv = 1.0/o.v
            return Dual(self.v*iv, (self.d - o.d*(self.v*iv))*iv)
        return Dual(self.v/o, self.d/o)

    def __rtruediv__(self, o):
        iv = 1.0/self.v
        return Dual(o*iv, -self.d*(o*iv*iv))


def var(z):
    """z: (7, N) complex -> list of 7 Duals."""
    N = z.shape[1]
    out = []
    for i in range(NV):
        d = np.zeros((NV, N), complex); d[i] = 1.0
        out.append(Dual(z[i].astype(complex), d))
    return out


def val(x):
    return x.v if isinstance(x, Dual) else x


def csqrt(x, ref):
    """sqrt with the sign chosen closest to ref (continuation by proximity)."""
    s = np.sqrt(val(x) + 0j)
    if ref is not None:
        flip = np.abs(s - ref) > np.abs(-s - ref)
        s = np.where(flip, -s, s)
    if isinstance(x, Dual):
        return Dual(s, x.d/(2*s))
    return s


# ----------------------------------------------------------------------------------------------- lightcone geometry (l1_crystal)
def dot(p, q):
    return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]


def sub(p, q):
    return [p[i] - q[i] for i in range(3)]


def msq(p, q):
    r = sub(p, q); return -dot(r, r)


def geometry(z, y0ref=None, n0=1.0, y0val=None, xsqval=None):
    """z: list of 7 (Dual or arrays). Returns dict of the rational building blocks (same objects as l1_crystal.build).
    y0val: y0 given explicitly (e.g. multiprecision evaluation); xsqval: x^2 given explicitly (avoids the cancellation at s -> +-1)."""
    xt0, ux, vx, uy, vy, xp, yp = z
    x = [ux, vx, xp]; y = [uy, vy, yp]
    uv = uy*vy; yp2 = yp*yp
    y0sq = uv*(uv*(1.0 - xt0*xt0) - yp2)/(uv - yp2)
    if y0val is not None:
        y0 = y0val
    else:
        y0 = csqrt(y0sq, None if y0ref is None else n0*y0ref)
        if y0ref is None:
            y0 = y0*n0
    u0 = uy*xt0 + y0*yp/vy
    v0 = -(vy*xt0) + y0*yp/uy
    x0 = [u0, v0, y0]
    xsq = dot(x, x) if xsqval is None else xsqval
    rx = 1.0/(-xsq)
    x4 = [rx*c for c in x]; x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    one = np.ones_like(val(ux)); zer = np.zeros_like(val(ux))
    e = [one, one, zer]; z0 = [one, zer, zer]
    x10, x30, x20, x40 = sub(x1, x0), sub(x3, x0), sub(x2, x0), sub(x4, x0)
    m24, m20, m40 = msq(x2, x4), msq(x2, x0), msq(x4, x0)
    Nn = 2*(dot(z0, x40)*(-m20) - dot(z0, x20)*(-m40))
    c13 = 2*dot(x10, x30); c10 = 2*dot(z0, x10); c30 = -2*dot(z0, x30)
    m35, m45, m12, m15, m25, m34 = msq(x3, e), msq(x4, e), msq(x1, x2), msq(x1, e), msq(x2, e), msq(x3, x4)
    m13, m14, m23 = msq(x1, x3), msq(x1, x4), msq(x2, x3)
    xy = dot(x, y); x2y2 = xsq*dot(y, y)
    FP1 = 1.0 + 2*xy + x2y2; FP2 = 1.0 - 2*xy + x2y2; FP3 = 1.0 - x2y2
    vmuN = 8.0*dot(x4, y)*(dot(x4, x4) + dot(y, y))          # (v'-u') m13 m24  (x3 = -x1, x2 = -x4)
    return dict(x=x, y=y, x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, e=e, z0=z0, y0=y0, y0sq=y0sq, uv=uv, yp2=yp2, xt0=xt0,
                negxsq=-xsq, m24=m24, m20=m20, m40=m40, Nn=Nn, c13=c13, c10=c10, c30=c30, m35=m35, m45=m45, m12=m12, m15=m15,
                m25=m25, m34=m34, m13=m13, m14=m14, m23=m23, FP1=FP1, FP2=FP2, FP3=FP3, vmuN=vmuN,
                mA=uv*(1.0 - xt0*xt0) - yp2, mB=uv - yp2)


def scalar_factors(nu, Delta):
    """(name, exponent) list: Phi_z(scalar, J = 0) = sum e_k log f_k  (l1_crystal.build: lamF + G; constants dropped).
    A = Nn^2 c13/(m24 m40 m20 c10 c30);  B = m35 m45 m12/(m15 m25 m34);  fsub = (v'-u')/(u'v') = vmuN m13 m24/(m12 m34 m14 m23)."""
    iv = 0.5j*nu; hD = 0.5*Delta
    E = {}
    def add(k, c):
        E[k] = E.get(k, 0.0) + c
    # lamF
    add('Nn', 2*iv); add('c13', iv); add('m24', -iv); add('m40', -iv); add('m20', -iv); add('c10', -iv); add('c30', -iv)
    add('m35', hD); add('m45', hD); add('m12', hD); add('m15', -hD); add('m25', -hD); add('m34', -hD)
    # G: meas, xsq, FP, dual, cel, sub, 3pt
    add('uv', 0.5); add('mB', -0.5); add('mA', -0.5)
    add('negxsq', -6.0)
    add('FP1', 1.0); add('FP2', 1.0); add('FP3', 1.0)
    add('Nn', -0.5); add('m24', 0.75); add('m20', -0.25); add('m40', -0.25)
    add('c13', 1.25); add('c10', -0.25); add('c30', -0.25)
    add('vmuN', 1.0); add('m13', 1.0); add('m24', 1.0); add('m12', -1.0); add('m34', -1.0); add('m14', -1.0); add('m23', -1.0)
    add('m12', -3.0); add('m34', -1.5); add('m35', -1.5); add('m45', -1.5)
    return [(k, v) for k, v in E.items() if v != 0]


LOGV = (1, 2, 3, 4)


def zmap(w, c):
    """w (7, N) complex -> raw z (7, N), logJac (N), dz/dw (7, N)  (n1_direct.zmap with TANH0)."""
    z = np.empty_like(w); dz = np.empty_like(w)
    lj = np.zeros(w.shape[1], complex)
    for i in range(7):
        if i in LOGV:
            z[i] = c[i]*np.exp(w[i]); dz[i] = z[i]; lj = lj + w[i] + np.log(c[i])
        elif i == 0:
            e2 = np.exp(2*w[0]); z[0] = (e2 - 1)/(e2 + 1); dz[0] = 4*e2/(e2 + 1)**2; lj = lj + np.log(4*e2/(e2 + 1)**2)
        else:
            z[i] = w[i]; dz[i] = 1.0
    return z, lj, dz


def dlogjac(w):
    g = np.zeros_like(w)
    for i in LOGV:
        g[i] = 1.0
    e2 = np.exp(2*w[0]); g[0] = 2.0 - 4*e2/(e2 + 1)          # d/dw0 log(4 e2/(e2+1)^2)
    return g


class Flow:
    """scalar flow function at (nu, Delta), J = 0, n0 = +1, in the w-coordinates of crystal c."""

    def __init__(self, nu, Delta, c):
        self.nu, self.Delta, self.c = float(nu), float(Delta), np.asarray(c, complex)
        self.fac = scalar_factors(nu, Delta)

    def grad(self, w, y0ref):
        """dPhi/dw (7, N) and the geometry values; y0ref: previous y0 (N) for the sqrt continuation."""
        z, lj, dz = zmap(w, self.c)
        g = geometry(var(z), y0ref)
        gz = np.zeros_like(w)
        for k, e in self.fac:
            f = g[k]
            gz = gz + e*f.d/f.v
        gw = gz*dz + dlogjac(w)
        vals = {k: val(v) for k, v in g.items() if not isinstance(v, list)}
        return gw, vals

    def logabs_parts(self, vals):
        """sum_k Re(e_k) log|f_k| (the part of Re Phi without the arg terms; arg terms come only from i nu/2 log A)."""
        out = 0.0
        for k, e in self.fac:
            out = out + np.real(e)*np.log(np.abs(vals[k]))
        return out

    def logA(self, vals):
        return 2*np.log(vals['Nn']) + np.log(vals['c13']) - np.log(vals['m24']) - np.log(vals['m40']) - np.log(vals['m20']) \
            - np.log(vals['c10']) - np.log(vals['c30'])

    def logB(self, vals):
        return np.log(vals['m35']) + np.log(vals['m45']) + np.log(vals['m12']) - np.log(vals['m15']) - np.log(vals['m25']) \
            - np.log(vals['m34'])

    def hess(self, w, y0ref, h=1e-4):
        """holomorphic Hessian d^2 Phi/dw^2 (7, 7) at ONE point w (7,), 4-point stencil on the exact gradient."""
        om = np.array([1, 1j, -1, -1j])
        W = np.repeat(w[:, None], 28, axis=1)
        for i in range(7):
            for j in range(4):
                W[i, 4*i + j] += h*om[j]
        G, _ = self.grad(W, np.full(28, y0ref))
        H = np.zeros((7, 7), complex)
        for i in range(7):
            H[:, i] = sum(G[:, 4*i + j]/om[j] for j in range(4))/(4*h)
        return 0.5*(H + H.T)


def crystal(nu, Delta):
    a = (Delta - np.sqrt(Delta**2 - nu**2))/nu
    return np.array([0.0, -1j*a, 1j*a, 1j*a, -1j*a, 0.0, 0.0], complex), a


# ----------------------------------------------------------------------------------------------- EE / scalar ratio Q
_RE = None


def _repro():
    global _RE
    if _RE is None:
        sys.path.insert(0, os.path.join(ADS, 'repro_engine'))
        import integrand as RE                       # read-only import (own numpy code of the reproducer)
        _RE = RE
    return _RE


def Qratio(gv, Delta, parts=False):
    """Q = Phi_EE/Phi_scal = m12^-2 m34^-2 V042^-4 pref ((v'-u')^2 - 1)/4 at the points of geometry values gv (plain arrays).
    pref/V042^4 is formed as F.Ehat.Ehat.G with Ehat = E/v042^2 (exact; avoids the v042 cancellation)."""
    RE = _repro()
    x1, x2, x3, x4, x0, e, z0 = gv['x1'], gv['x2'], gv['x3'], gv['x4'], gv['x0'], gv['e'], gv['z0']
    x10, x30 = sub(x1, x0), sub(x3, x0)
    F = RE.tto_form(x1, x10, x2, e, RE.coef(Delta))
    Gt = RE.tto_form(x3, x30, x4, e, RE.coef(3.0 - Delta))
    Fu = RE.traceless_raise(F); Gu = RE.traceless_raise(Gt)
    inv = 1.0/(RE.sq(x2, x0)*RE.sq(x4, x0))
    BAS = RE.BAS
    H24 = np.stack([np.stack([RE.Hf(x2, list(BAS[m]), x4, list(BAS[k])) for k in range(3)]) for m in range(3)])
    V204 = RE.cov_lin(lambda Bv: RE.Vf(x2, [Bv[0], Bv[1], Bv[2]], x0, x4))
    V402 = RE.cov_lin(lambda Bv: RE.Vf(x4, [Bv[0], Bv[1], Bv[2]], x0, x2))
    Eh = inv*(H24 - 2*np.einsum('mn,kn->mkn', V204, V402))
    prefh = np.einsum('mnz,mkz,nlz,klz->z', Fu, Eh, Eh, Gu)         # = pref / v042^4
    vmu = gv['vmuN']/(gv['m13']*gv['m24'])
    Q = prefh*(vmu*vmu - 1.0)/(4.0*gv['m12']**2*gv['m34']**2)
    if parts:
        return Q, dict(prefh=prefh, vmu=vmu)
    return Q


def geom_plain(z, y0ref=None):
    """geometry with plain complex arrays; z (7, N)."""
    g = geometry([z[i].astype(complex) for i in range(7)], y0ref)
    return g


def raw_from_D(X1, X2, s, Y1, Y2, phi, al, Delta, n0=1.0):
    """raw z (7, N) of the real-section D point (bf_core conventions, psi = alpha - n0 phi)."""
    a = 1.0/Delta
    psi = al - n0*phi
    rho = np.sqrt(Y1*Y2 + 0j)
    z = np.empty((7, np.size(X1)), complex)
    z[0] = np.cos(phi)*np.cos(psi); z[1] = a*X1; z[2] = a*X2; z[3] = a*Y1; z[4] = a*Y2
    z[5] = a*s*np.sqrt(X1*X2 + 0j); z[6] = a*rho*np.sin(phi)
    y0 = n0*a*rho*np.sin(psi)
    return z, y0


def Qratio_mp(z, y0, Delta, dps=60):
    """Q at ONE raw point (z (7,) complex, y0 complex) in multiprecision (mpmath, dps digits): the double-precision contraction
    cancels catastrophically near s -> +-1 (x^2 -> 0; repro_engine E3) and near u_y v_y - y_p^2 -> 0."""
    import mpmath as mp
    with mp.workdps(dps):
        zo = [np.array([mp.mpc(complex(v))], dtype=object) for v in z]
        g = geometry(zo, y0val=np.array([mp.mpc(complex(y0))], dtype=object))
        Q = Qratio(g, mp.mpf(Delta))
        return complex(Q[0])


def Q_condition(z):
    """crude double-precision loss estimate of Qratio at raw points z (7, N): ~ 1e-16 * max(1/|w_x|, 1/|w_y|)^4,
    w_x = 1 - xp^2/(ux vx) (x^2 / (-ux vx)), w_y = 1 - yp^2/(uy vy)."""
    wx = np.abs(1.0 - z[5]**2/(z[1]*z[2])); wy = np.abs(1.0 - z[6]**2/(z[3]*z[4]))
    return 1e-16*np.maximum(1.0/wx, 1.0/wy)**4
