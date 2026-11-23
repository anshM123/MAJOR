"""Continuum k = 2 de-shift weights (CERT_CLOSE.md §3): Result D's Omega^EE (NF02 psi2/k0) and Result Q's smeared k = 2 piece (Q3 psi2),
built in the SHIFTED form of DESHIFT_CLOSE §2.1, k = 2 version.

Identity (chi-normalisation of prl/data/README §4; X = chi (R_- + R_+)): for real even chi,
    W[chi](t) := int_0^oo chi [R_- cos((nu + i)t) + R_+ cos((nu - i)t)] dnu  =  Phi(t) := int_0^oo psi(nu) cos(nu t) dnu.
R_+(nu) = R_-(-nu) (exact), so the integrand is even and W = 1/2 int_R chi R_- cos((nu+i)t) + (R_+ term).  Shift nu = mu - i (R_- term) and
nu = mu + i (R_+ term).  If chi R_- is analytic in -1 < Im nu < 0 (continuous up to the boundary, integrable, decaying), then for EVERY real t
    W(t) = 1/2 int_R E(mu) cos(mu t) dmu,      E(mu) = f(mu) + f(-mu) = 2 Re f(mu),   f(mu) := chi(mu - i) R_-(mu - i),
and no residue term arises (R_-'s only pole in the closed lower strip is nu = 0, cancelled by chi(0) = 0; R_- has a DOUBLE ZERO at
-i/2, so chi may carry poles of order <= 2 at +-i/2: these are the "common zeros" of k = 2).  Hence, with psi_e(mu) := psi(|mu|),
    |W(t) - Phi(t)| <= eps_W := int_0^oo |E(mu) - psi(mu)| dmu      for all t >= 0  (and |image - target| <= eps_W P_{-1/2}(cosh beta)).
The target psi(|mu|) has a |mu| kink at 0 (psi'(0) = kappa != 0), which no chi analytic in the closed strip can produce: the exact solution
is singular at nu = +-i.  We include the singular functions  s_k(nu) = nu^(2k) l(nu) exp(-nu^2/(2 s^2)),  l = (nu^2+1) log(nu^2+1)
(principal branch; analytic in |Im nu| < 1, branch points at +-i), whose boundary values reproduce |mu| and mu^2 log|mu| exactly.
Weight:  chi(nu) = Q(nu) B(nu) + sum_k y_k s_k(nu),   Q = nu^2/(nu^2 + 1/4)^m,   B = Gaussian RBFs (even) on two scales.
Coefficients by weighted least squares (E-rows on a mu-grid, weak gauge rows towards the stored grid weight's real part); then eps_W
by a fine quadrature of |E - psi| (DESHIFT-grade [DV]: high-resolution quadrature, not interval-verified).
The imaginary part of the stored chi is a null vector of the identity (W[Im chi] = 0); it is fitted separately with target 0.
Usage: python3 k2_cont.py D|Q DG [m]   -> cache/k2_<D|Q>_<DG>.npz, log to stdout
"""
import os, sys, json, pickle, time
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'deshift'))
import kern                                                  # noqa: E402  (read-only: exact R_-+ as data)

K0 = 27.48509817
CACHE = os.path.join(HERE, 'cache')


def lam_of(Dg): return np.sqrt((Dg - 2.0) * (Dg - 1.0))


class Target:
    """psi(nu) on nu >= 0 (factored evaluation near p = 1), kink slope kappa = psi'(0)."""
    def __init__(self, which, Dg):
        self.Lam = lam_of(Dg); self.which = which
        if which == 'D':
            d = json.load(open(os.path.join(ADS4, 'code', 'robust', 'functionals_noforward.json')))['NF02']
            self.r = d['r']; self.co = np.array([float(Fr(x)) for x in d['X']]) / K0      # psi = (1-p)^5 p X(p) / k0
            self.exact = [Fr(x) / Fr('27.48509817') for x in d['psi2_exact']]
        else:
            d = json.load(open(os.path.join(ADS4, 'prl', 'data', 'q3_full_delta0.2.json')))
            c = np.array([float(x) for x in d['psi2']])
            # psi2 = (1-p)^3 * cof(p): exact polynomial division (the stored psi2 has the (1-p)^3 factor up to rounding)
            q, rem = np.polynomial.polynomial.polydiv(c, np.polynomial.polynomial.polypow([1, -1], 3))
            self.r = 3; self.co_full = c; self.co = None; self.q3 = q; self.rem = rem
        self.kappa = self.dpsi0()

    def __call__(self, nu):
        nu = np.abs(np.asarray(nu, float)); p = nu / self.Lam
        if self.which == 'D':
            v = (1 - p) ** 5 * p * np.polynomial.polynomial.polyval(p, self.co)
        else:
            v = np.polynomial.polynomial.polyval(p, self.co_full)
        return np.where(p <= 1, v, 0.0)

    def dpsi0(self):
        if self.which == 'D':
            return self.co[0] / self.Lam
        return self.co_full[1] / self.Lam


class Basis:
    def __init__(self, Lam, m=2, fine=(0.15, 0.45, 40.0), coarse=(0.3, 0.9, None), s_sing=1.0, nsing=3):
        hf, sf, nf = fine; hc, sc, _ = coarse
        numax = Lam + 30.0
        self.cent = [(k * hf, sf) for k in range(0, int(round(nf / hf)) + 1)]
        k0 = int(np.floor(nf / hc)) + 1
        self.cent += [(k * hc, sc) for k in range(k0, int(np.ceil(numax / hc)) + 1)]
        self.m = m; self.ss = s_sing; self.nsing = nsing
        self.n = len(self.cent) + nsing

    def Q(self, z):
        return z * z / (z * z + 0.25) ** self.m

    def eval(self, z):
        """matrix (len(z), n) of basis functions at complex z."""
        z = np.asarray(z, complex)
        out = np.empty((len(z), self.n), complex)
        q = self.Q(z)
        for i, (c, s) in enumerate(self.cent):
            g = np.exp(-(z - c) ** 2 / (2 * s * s))
            if c > 0: g = g + np.exp(-(z + c) ** 2 / (2 * s * s))
            out[:, i] = q * g
        l = (z * z + 1) * np.log(z * z + 1)
        e = np.exp(-z * z / (2 * self.ss ** 2))
        for k in range(self.nsing):
            out[:, len(self.cent) + k] = z ** (2 * k) * l * e
        return out


def mu_grid(mumax):
    g = np.concatenate([np.geomspace(1e-5, 0.5, 300), np.arange(0.5, 40.0, 0.02), np.arange(40.0, mumax, 0.1)])
    return np.unique(g)


def quad_weights(x):
    w = np.zeros_like(x); dx = np.diff(x); w[:-1] += dx / 2; w[1:] += dx / 2
    return w


def stored(which, Dg):
    if which == 'D':
        A = np.loadtxt(os.path.join(ADS4, 'prl', 'data', 'd_abel_X_%d.txt' % Dg))
        nu, X = A[:, 0], A[:, 1] + 1j * A[:, 2]
    else:
        d = pickle.load(open(os.path.join(ADS4, 'code', 'resultq', 'cache', 'rv_k2eng_Y_%d.pkl' % Dg), 'rb'))
        nu, X = d['nu'], d['a'] * d['Y']
    Rm, Rp = kern.load_D()
    chi = X / (Rm.np(nu) + Rp.np(nu))
    return nu, chi, X


def build(which, Dg, m=2, lam_g=1e-2, part='re', chunk=2500, **kw):
    """weighted least squares by a blocked (TSQR) Householder QR of [A | b], column-scaled; memory O(n^2)."""
    import scipy.linalg as sl
    t0 = time.time()
    T = Target(which, Dg); Lam = T.Lam
    Rm, Rp = kern.load_D()
    B = Basis(Lam, m=m, **kw)
    mumax = Lam + 60.0
    mu = mu_grid(mumax); wq = quad_weights(mu)
    tgt = T(mu) if part == 're' else 0 * mu
    sw = np.sqrt(wq)
    nu_s, chi_s, X_s = stored(which, Dg)
    sel = nu_s <= Lam + 25
    st_ = 4 if Dg < 500 else 8
    nug = nu_s[sel][::st_]; cg = (chi_s.real if part == 're' else chi_s.imag)[sel][::st_]
    # beyond the stored grid the weight is pinned to 0 (the stored weight is < 1e-14 of its maximum there)
    nz = np.arange(nu_s[-1] + 0.1, Lam + 60.0, 0.1)
    nug = np.concatenate([nug, nz]); cg = np.concatenate([cg, 0 * nz])
    sc = np.max(np.abs(cg)); dg = np.sqrt(nug[1] - nug[0])

    def blocks():
        for i in range(0, len(mu), chunk):
            z = mu[i:i + chunk] - 1j
            yield 2 * (B.eval(z) * Rm.np(z)[:, None]).real * sw[i:i + chunk, None], tgt[i:i + chunk] * sw[i:i + chunk]
        for i in range(0, len(nug), chunk):
            yield B.eval(nug[i:i + chunk].astype(complex)).real * (lam_g * dg / sc), cg[i:i + chunk] * (lam_g * dg / sc)
    cn = np.zeros(B.n)
    for Ab, bb in blocks():
        cn += np.sum(Ab * Ab, axis=0)
    cn = np.sqrt(cn); cn[cn == 0] = 1
    Rr = None
    for Ab, bb in blocks():
        M = np.hstack([Ab / cn, bb[:, None]])
        if Rr is not None: M = np.vstack([Rr, M])
        Rr = sl.qr(M, mode='r', overwrite_a=True, check_finite=False)[0][:B.n + 1]
    Rn = Rr[:B.n, :B.n]; d = np.abs(np.diag(Rn))
    x = sl.solve_triangular(Rn, Rr[:B.n, B.n], check_finite=False) / cn
    rk = int(np.sum(d > 1e-15 * d.max()))
    Gx = np.concatenate([B.eval(nug[i:i + chunk].astype(complex)).real @ x for i in range(0, len(nug), chunk)])
    Ex = np.concatenate([2 * (B.eval(mu[i:i + chunk] - 1j) @ x * Rm.np(mu[i:i + chunk] - 1j)).real for i in range(0, len(mu), chunk)])
    return dict(T=T, B=B, x=x, mu=mu, E=Ex, tgt=tgt, rank=rk, nug=nug, cg=cg, Gx=Gx, t=time.time() - t0, Rm=Rm,
                diagR=(d.min(), d.max()), resid=abs(Rr[B.n, B.n]) if Rr.shape[0] > B.n else 0.0)


def eps_W(Bx, B, T, Rm, part, mumax, n=400000):
    """fine quadrature of |E - psi| on (0, mumax] (piecewise geometric/uniform grid) -- [DV], not interval-verified."""
    g = np.unique(np.concatenate([np.geomspace(1e-9, 1.0, 20000), np.linspace(1.0, mumax, n)]))
    tot = 0.0; mx = 0.0
    for i in range(0, len(g), 4000):
        mu = g[i:i + 4001]
        z = mu - 1j
        E = 2 * (B.eval(z) @ Bx * Rm.np(z)).real
        d = np.abs(E - (T(mu) if part == 're' else 0))
        tot += np.sum((d[1:] + d[:-1]) / 2 * np.diff(mu)); mx = max(mx, d.max())
    return tot, mx


def chi_eval(B, x, nu, chunk=4000):
    out = []
    nu = np.asarray(nu, float)
    for i in range(0, len(nu), chunk):
        out.append((B.eval(nu[i:i + chunk].astype(complex)) @ x).real)
    return np.concatenate(out)


def checks(which, Dg, xr, xi, B, T):
    """V3 graviton-pole identity and dk0 for chi = chi_r + i chi_i (fine real-axis quadrature), vs the stored grid weight."""
    Rm, Rp = kern.load_D()
    # composite Gauss-Legendre (16 nodes per panel of width 0.05; panel edges at 0 and Lam): smooth integrands on each panel
    from numpy.polynomial.legendre import leggauss
    xg, wg = leggauss(16)
    edges = np.unique(np.concatenate([np.arange(0, T.Lam, 0.05), [T.Lam], np.arange(T.Lam, T.Lam + 60, 0.05)]))
    lo, hi = edges[:-1], edges[1:]
    nu = ((hi - lo)[:, None] * (xg[None, :] + 1) / 2 + lo[:, None]).ravel(); dnu = ((hi - lo)[:, None] * wg[None, :] / 2).ravel()
    chr_ = chi_eval(B, xr, nu); chi_ = chi_eval(B, xi, nu)
    ch = chr_ + 1j * chi_
    Dz = lambda z: 1 / (z * z + 2.25)
    hm, hp = Rm.np(nu), Rp.np(nu)
    V3 = np.sum(dnu * ch * (hm * Dz(nu + 1j) + hp * Dz(nu - 1j))); V3t = np.sum(dnu * T(nu) * Dz(nu))
    X = ch * (hm + hp); p = nu / T.Lam
    dk = (-2 * np.sum(dnu * (X.real - T(nu)) * p ** 6) / T.Lam)
    nu_s, chi_s, X_s = stored(which, Dg)
    ds = nu_s[1] - nu_s[0]
    V3s = np.sum(ds * chi_s * (Rm.np(nu_s) * Dz(nu_s + 1j) + Rp.np(nu_s) * Dz(nu_s - 1j))); V3st = np.sum(ds * T(nu_s) * Dz(nu_s))
    dks = -2 * np.sum((X_s.real - T(nu_s)) * (nu_s / T.Lam) ** 6) * ds / T.Lam
    chs = chi_eval(B, xr, nu_s) + 1j * chi_eval(B, xi, nu_s)
    dX = np.abs(chs * (Rm.np(nu_s) + Rp.np(nu_s)) - X_s) / T(np.linspace(0, T.Lam, 100001)).max()
    return dict(V3=V3, V3t=V3t, dk=dk, V3s=V3s, V3st=V3st, dks=dks, dXmax=[(lo, hi, dX[(nu_s >= lo) & (nu_s < hi)].max())
                for lo, hi in ((0, 2), (2, 10), (10, 25), (25, T.Lam), (T.Lam, nu_s[-1] + 1))])


if __name__ == '__main__':
    which, Dg = sys.argv[1], int(sys.argv[2]); m = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    xs = {}
    for part in ['re', 'im']:
        R = build(which, Dg, m=m, part=part)
        T, B = R['T'], R['B']
        xs[part] = R['x']
        ew, mx = eps_W(R['x'], B, T, R['Rm'], part, T.Lam + 60.0)
        if part == 're': ewr = ew
        phi_norm = np.sum(np.abs(T(np.linspace(0, T.Lam, 200001)))) * T.Lam / 200000
        dch = np.max(np.abs(R['Gx'] - R['cg'])) / np.max(np.abs(R['cg']))
        print('%s Dg=%d part=%s m=%d: unknowns %d, |diag R| in [%.1e, %.1e]; eps_W = %.3e (int|psi| = %.4g; rel %.2e), max|E-psi| = %.2e; '
              'max|chi - chi_stored|/max = %.2e; kappa = %.6g [%.0fs]' % (which, Dg, part, m, B.n, R['diagR'][0], R['diagR'][1], ew, phi_norm, ew / phi_norm, mx,
                                                                      dch, T.kappa, R['t']), flush=True)
        np.savez(os.path.join(CACHE, 'k2_%s_%d_%s_m%d.npz' % (which, Dg, part, m)), x=R['x'], cent=np.array(B.cent), m=m, ss=B.ss,
                 nsing=B.nsing, eps_W=ew, maxerr=mx, intpsi=phi_norm)
    C = checks(which, Dg, xs['re'], xs['im'], B, T)
    print('%s Dg=%d checks: V3 continuum %.12e vs target %.12e (rel %.1e; bound (4/9) eps_W/V3 = %.1e) | stored grid %.12e vs %.12e (rel %.1e)' % (
        which, Dg, C['V3'].real, C['V3t'], abs(C['V3'] - C['V3t']) / C['V3t'], 4 / 9 * ewr / C['V3t'], C['V3s'].real, C['V3st'], abs(C['V3s'] - C['V3st']) / C['V3st']))
    print('%s Dg=%d checks: dk0 continuum %+.9e | stored %+.9e (rel diff %.1e); (1+dk0)^-1/2 = %.8f vs %.8f' % (
        which, Dg, C['dk'], C['dks'], abs(C['dk'] - C['dks']) / abs(C['dks']), (1 + C['dk']) ** -0.5, (1 + C['dks']) ** -0.5))
    print('%s Dg=%d checks: max|X_cont - X_stored|/max psi on ' % (which, Dg) + '; '.join('[%.0f,%.0f): %.1e' % t for t in C['dXmax']), flush=True)
