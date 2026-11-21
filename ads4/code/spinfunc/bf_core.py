"""
bf_core (BELOW_FOLD.md): the CLS-D crystal integral evaluated WITHOUT a saddle, on the real "future-cone" section D, by
randomised QMC with an importance density that is the exact separable leading-Regge (Delta -> oo) density.

Section D (raw l1_crystal variables z = (xt0, ux, vx, uy, vy, xp, yp), all REAL), with a = 1/Delta (Regge scaling):
    ux = a X1, vx = a X2, xp = a s sqrt(X1 X2)          X1, X2 > 0, s in (-1, 1)      (x future-timelike)
    uy = a Y1, vy = a Y2, yp = a rho sin(phi)           Y1, Y2 > 0, rho = sqrt(Y1 Y2), phi in (-pi/2, pi/2)
    xt0 = cos(phi) cos(psi), y0 = n0 a rho sin(psi)     psi in (0, pi), n0 = +-1 (summed)      (celestial sphere)
    d^7 z = a^6 sqrt(X1 X2) rho cos^2(phi) sin(psi) dX1 dX2 ds dY1 dY2 dphi dpsi.
On D every factor of the integrand is real and A > 0 (so A^{i nu/2} is a pure phase), B^{Delta/2} -> exp(-(X1+X2+Y1+Y2)):
the integral is absolutely convergent, non-oscillatory at nu -> 0, and has no saddle (hence no fold).
Leading-Regge density (bf_lead.py, checked against the code to 1e-7): q_x^{7/4} X2^{-1/2} (X1Y2 + X2Y1 - 2 Xp Yp) Y1^{9/4} Y2^{7/4}
cos^{9/2}(phi) |sin(psi + n0 phi)|^{-1/2} e^{-S} (times phases); the sampler below draws from a product density close to it.

The integrand itself is the project's exact integrand, evaluated point-wise by n1_pw (l1_crystal.build for the scalar,
l9_points.build_spin for the graviton kinds), read-only.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.stats import qmc, gamma as GD, beta as BD
from scipy.special import gammaln, betaln

# sampling exponents (shape parameters of the product density)
KX1, KX2, KY1, KY2 = 4.0, 3.0, 4.0, 3.5       # Gamma shapes
BS = 1.75                                      # s ~ (1-s^2)^BS
CPHI = 4.5                                     # phi: Beta(CPHI+1, CPHI+1) on (phi + pi/2)/pi
APSI = 0.5                                     # alpha = psi + n0 phi (mod pi): pi * Beta(APSI, APSI)


class Sample:
    """QMC sample of D (scaled variables) with log importance density; m = log2(points), seed = scrambling seed."""

    def __init__(self, m, seed, kx=(KX1, KX2), ky=(KY1, KY2), theta=0.0, scale=1.0, cphi=CPHI):
        u = qmc.Sobol(d=7, scramble=True, seed=seed).random_base2(m)
        u = np.clip(u, 1e-14, 1 - 1e-14)
        self.theta = float(theta); sc = float(scale)/np.cos(self.theta)     # radial rays rotated by -+theta: |e^{-X}| = e^{-r cos theta}
        self.sc = sc
        self.X1 = GD.ppf(u[:, 0], kx[0])*sc; self.X2 = GD.ppf(u[:, 1], kx[1])*sc
        b = BD.ppf(u[:, 2], BS + 1, BS + 1); self.s = 2*b - 1
        self.Y1 = GD.ppf(u[:, 3], ky[0])*sc; self.Y2 = GD.ppf(u[:, 4], ky[1])*sc
        c = BD.ppf(u[:, 5], cphi + 1, cphi + 1); self.phi = np.pi*(c - 0.5)
        self.al = np.pi*BD.ppf(u[:, 6], APSI, APSI)
        lp = (GD.logpdf(self.X1/sc, kx[0]) + GD.logpdf(self.X2/sc, kx[1]) + BD.logpdf(b, BS + 1, BS + 1) - np.log(2.0)
              + GD.logpdf(self.Y1/sc, ky[0]) + GD.logpdf(self.Y2/sc, ky[1]) - 4*np.log(sc) + BD.logpdf(c, cphi + 1, cphi + 1) - np.log(np.pi)
              + BD.logpdf(self.al/np.pi, APSI, APSI) - np.log(np.pi))
        self.logp = lp
        self.n = u.shape[0]

    def raw(self, a, n0):
        """raw z (n x 7, complex) at scale a for sheet n0; and log of the Jacobian d^7z / (scaled measure) (a^6 dropped)."""
        psi = np.mod(self.al - n0*self.phi, np.pi)
        rho = np.sqrt(self.Y1*self.Y2)
        sq = np.sqrt(self.X1*self.X2)
        Z = np.empty((self.n, 7), complex)
        Z[:, 0] = np.cos(self.phi)*np.cos(psi)
        e = np.exp(1j*self.theta)
        Z[:, 1] = a*self.X1/e; Z[:, 2] = a*self.X2*e; Z[:, 3] = a*self.Y1*e; Z[:, 4] = a*self.Y2/e
        Z[:, 5] = a*self.s*sq; Z[:, 6] = a*rho*np.sin(self.phi)
        lj = np.log(sq) + np.log(rho) + 2*np.log(np.cos(self.phi)) + np.log(np.sin(psi))
        return Z, lj, psi


class SampleR(Sample):
    """as Sample, but the radial pairs in light-cone 'rapidity' form: X1 = r e^{eta}, X2 = r e^{-eta} (same for Y), so that the boost
    direction (where P_J(t) ~ e^{J|eta_x - eta_y|} and graviton polarisation factors grow) is sampled with a heavy (logistic) tail:
    r ~ Gamma(kr), eta ~ logistic(mu, s).  dX1 dX2 = 2 r dr deta."""

    def __init__(self, m, seed, kr=(6.0, 6.0), mu=(0.25, 0.25), sc=(0.6, 0.6), theta=0.0):
        u = qmc.Sobol(d=7, scramble=True, seed=seed).random_base2(m)
        u = np.clip(u, 1e-14, 1 - 1e-14)
        self.theta = float(theta); c = 1.0/np.cos(self.theta)
        lp = 0.0
        out = []
        for k, (iu, ie) in enumerate(((0, 1), (3, 4))):
            r = GD.ppf(u[:, iu], kr[k])*c/2
            eta = mu[k] + sc[k]*np.log(u[:, ie]/(1 - u[:, ie]))
            z = (eta - mu[k])/sc[k]
            lp = lp + GD.logpdf(2*r/c, kr[k]) + np.log(2/c) - z - 2*np.log1p(np.exp(-z)) - np.log(sc[k]) - np.log(2*r)
            out.append((r*np.exp(eta), r*np.exp(-eta)))
        (self.X1, self.X2), (self.Y1, self.Y2) = out
        b = BD.ppf(u[:, 2], BS + 1, BS + 1); self.s = 2*b - 1
        cc = BD.ppf(u[:, 5], CPHI + 1, CPHI + 1); self.phi = np.pi*(cc - 0.5)
        self.al = np.pi*BD.ppf(u[:, 6], APSI, APSI)
        self.logp = (lp + BD.logpdf(b, BS + 1, BS + 1) - np.log(2.0) + BD.logpdf(cc, CPHI + 1, CPHI + 1) - np.log(np.pi)
                     + BD.logpdf(self.al/np.pi, APSI, APSI) - np.log(np.pi))
        self.n = u.shape[0]


_PW = None


def pw():
    global _PW
    if _PW is None:
        import n1_pw
        _PW = n1_pw
    return _PW


def _ref_point(a, theta=0.0):
    """a fixed interior point of D_theta (index 0 of every batch: all logs are continued from it)."""
    e = np.exp(1j*theta)
    z = np.array([0.3, a*2.0/e, a*1.5*e, a*2.5*e, a*1.7/e, a*0.2, a*0.3], complex)
    return z


def scalar_parts(S, Dl, J, n0, chunk=1 << 14):
    """per point: lA, lB, G (+ log Jacobian), real section; J enters G through Q_J (3-pt powers cancel)."""
    PW = pw(); a = 1.0/Dl
    Z, lj, psi = S.raw(a, n0)
    out = {k: np.empty(S.n, complex) for k in ('lA', 'lB', 'G')}
    z0 = _ref_point(a, S.theta)
    for c0 in range(0, S.n, chunk):
        zz = np.vstack([z0[None, :], Z[c0:c0 + chunk]])
        PW.set_batch(zz.shape[0])
        lF, G, P = PW.Cr.build([zz[:, i].copy() for i in range(7)], 1.0, Dl, J, n0=float(n0), parts=True)
        out['lA'][c0:c0 + chunk] = P['lA'][1:]; out['lB'][c0:c0 + chunk] = P['lB'][1:]
        out['G'][c0:c0 + chunk] = G[1:] + lj[c0:c0 + chunk]
    return out


_EE = [None]


def ee_build():
    """point-wise EE builder (w0_ee_common.make_build_even_EE on the n1_pw backend; n0 = +1 hard-wired there)."""
    if _EE[0] is None:
        import w0_ee_common as W0
        PW = pw()
        _EE[0] = W0.make_build_even_EE(PW, PW.S, PW.Cr, PW.A6)
    return _EE[0]


def grav_parts(S, Dl, J, kind='EE', chunk=1 << 12):
    """per point (sheet n0 = +1; the graviton builders are n0-even): lA, lB, Gs (+ log Jacobian), prefactor pref[a, c, pt] in the
    conserved-structure basis of l5_spin.conserved (1 x 1 for J = 0, 2), and the Lemma-B factor R (k = 4 = k = 2 x R)."""
    PW = pw(); a = 1.0/Dl
    import l5_spin as S5
    mons, sol, _ = S5.conserved(Dl, J)
    _, solt, _ = S5.conserved(3.0 - Dl, J)
    Z, lj, psi = S.raw(a, 1)
    z0 = _ref_point(a, S.theta)
    build = ee_build() if kind == 'EE' else (lambda zl, nu, D, J_, m_, s_, st_: PW.P9.build_spin(zl, nu, D, J_, m_, s_, st_, kind))
    lA = np.empty(S.n, complex); lB = np.empty(S.n, complex); G = np.empty(S.n, complex); pref = None
    for c0 in range(0, S.n, chunk):
        zz = np.vstack([z0[None, :], Z[c0:c0 + chunk]])
        PW.set_batch(zz.shape[0])
        zl = [zz[:, i].copy() for i in range(7)]
        lF, Gs, pr_ = build(zl, 1.0, Dl, J, mons, sol, solt)
        _, _, P = PW.Cr.build(zl, 1.0, Dl, 0, n0=1.0, parts=True)
        lA[c0:c0 + chunk] = P['lA'][1:]; lB[c0:c0 + chunk] = P['lB'][1:]
        G[c0:c0 + chunk] = Gs[1:] + lj[c0:c0 + chunk]
        if pref is None:
            pref = np.empty(pr_.shape[:2] + (S.n,), complex)
        pref[..., c0:c0 + chunk] = pr_[..., 1:]
    up, vp, R = cross_ratios(Z)
    return dict(lA=lA, lB=lB, G=G, pref=pref, R=R)


def t_ZW(Z):
    """t = Z.W/sqrt(Z^2 W^2) of l1_crystal (spin-J angle), point-wise from raw z (n x 7); J-dependence of the scalar = P_J(t)."""
    def dot(p, q):
        return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    ux, vx, uy, vy, xp, yp = (Z[:, i] for i in range(1, 7))
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = dot(x, x); x4 = [c/(-xsq) for c in x]; x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    e = [1.0, 1.0, 0.0]
    sub = lambda p, q: [p[i] - q[i] for i in range(3)]
    msq = lambda p, q: -dot(sub(p, q), sub(p, q))
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    x15, x25, x53, x54 = sub(x1, e), sub(x2, e), sub(e, x3), sub(e, x4)
    Zv = [m25*x15[i] - m15*x25[i] for i in range(3)]
    Wv = [m45*x53[i] - m35*x54[i] for i in range(3)]
    return dot(Zv, Wv)/np.sqrt(dot(Zv, Zv)*dot(Wv, Wv))


def cross_ratios(Z):
    """u', v' of x1..x4 (l1_crystal gauge) point-wise from raw z; Lemma B factor R = ((v'-u')^2 - 1)/(4u'v') (k = 4 = k = 2 x R)."""
    def dot(p, q):
        return -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    ux, vx, uy, vy, xp, yp = (Z[:, i] for i in range(1, 7))
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = dot(x, x); x4 = [c/(-xsq) for c in x]; x2 = [-c for c in x4]; x1 = [-c for c in y]; x3 = y
    sub = lambda p, q: [p[i] - q[i] for i in range(3)]
    msq = lambda p, q: -dot(sub(p, q), sub(p, q))
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    up = m12*m34/(m13*m24); vp = m14*m23/(m13*m24)
    return up, vp, ((vp - up)**2 - 1)/(4*up*vp)


def combine(parts, nus, Dl, logp, shift=None):
    """sum over points of exp((i nu/2) lA + (Dl/2) lB + G - logp - shift) for each nu; returns (vals[nu], shift)."""
    base = 0.5*Dl*parts['lB'] + parts['G'] - logp
    if shift is None:
        shift = np.max(base.real[np.isfinite(base.real)])
    w = np.exp(base - shift)
    bad = ~np.isfinite(w) | ~np.isfinite(parts['lA'])
    w = np.where(bad, 0.0, w)              # far region (raw |x| ~ 1, weight ~ e^{-Delta}): dropped
    lA = np.where(bad, 0.0, parts['lA'])
    vals = np.array([np.sum(w*np.exp(0.5j*nu*lA)) for nu in nus])
    return vals, shift
