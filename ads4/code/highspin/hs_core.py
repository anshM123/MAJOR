"""
hs_core (HIGHSPIN.md): shared machinery for the high-spin near-gap states of Omega^EE at Delta_g = 1000.

* Rows: the chain's K2 rows (EE_NEARGAP / ee1000 / neargap2 rows, read-only) with DIRECT n1 points (hs_points.py, highspin/cache, and
  the older ones in neargap2/cache) injected exactly as hc_chain's HC_INJ / neargap2's ng_val do (hc_nodes.to_record; direct records
  override K2 records at the same nu), Laplace order K (KK), jackknife replica JK.  Implemented by importing neargap2/ng_val (which
  imports ng_chain = the unchanged w2/q3 chain in EE mode with the ee1000 Jacobi fix) with hc_nodes' point directory moved to
  highspin/cache.  Every module is imported read-only; hs_guard forbids writes outside highspin/cache.
* Reference-series diagnostic (refdiag): |r2/r0|, |r4/r0| of the Delta_ref = 2e4 reference record at each K2 node -- the reference's
  leading partial wave D_h(nu/m_ref) has its first zero at nu0 = x0 m_ref (HIGHSPIN.md §2), where its Laplace series stops converging.
* Curvature: kappa image of SMOOTH_BOUND §4 (sm_curvfix.curv_image, imported read-only with SM_CURV=none so it patches nothing),
  subtracted from every diagonal Gamma (the published convention: the state's base already contains the leading-Regge curvature).
* Models of Gamma(x) on a row: 'shape' (ee1000_chain.ShapeFit, copied: Gamma = c(x) D_h + s(x) x D_h'/(2J+1) + r(x), c, s, r of degree
  K, held beyond the last node), 'spline' (natural cubic, held), 'lin'.
* Quadrature of a hybrid state's correction: composite Gauss-Legendre in p on [0, p_cut] and Gauss-Kronrod-15 on [p_cut, 1] (so that
  direct nodes placed at the GK15 x-nodes of an edge window are used without any interpolation); base and N_s from ee1000_common.hyb.
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import hs_guard  # noqa: F401  (first: registers itself as ng_guard)
NG = os.path.join(HERE, '..', 'neargap2')
SF = os.path.join(HERE, '..', 'spinfunc')
HSC = os.path.join(HERE, 'cache')
sys.path.insert(0, SF); sys.path.insert(0, NG); sys.path.insert(0, HERE)
os.environ['SM_CURV'] = 'none'                 # sm_curvfix imported for curv_image only (no patch)
import numpy as np
import ng_val as V                              # -> ng_chain -> ee1000_common -> w2_chain (EE mode)
import sm_curvfix as CF
NC, HN = V.NC, V.HN
E, Q, W = NC.E, NC.Q, NC.W
HN.C = HSC                                      # direct points: highspin/cache (neargap2 points are copied there by hs_points)
HN._SC.clear()
_records0 = HN.records


def _records_rekeyed(kind, J, D, jk=None, **kw):
    """hc_nodes.records with every nu moved to round(nu/m, 6) m (m = the state's mass): p4_odd_analysis.grid / l9 rebuild nu from x
    rounded to 6 decimals and match the reference record within 1e-7 nu, so a direct node at a non-round x (the GK15 nodes) would
    otherwise be silently dropped.  The shift is |d nu| <= 5e-7 m (a phase change <= 1e-5 of the oscillation)."""
    out = _records0(kind, J, D, jk, **kw)
    m = HN.mphys(D, J)
    return {(k[0], round(round(k[1]/m, 6)*m, 9)): v for k, v in out.items()}


HN.records = _records_rekeyed
DREF = 20000.0
GKX = np.array([0.991455371120813, 0.949107912342759, 0.864864423359769, 0.741531185599394, 0.586087235467691,
                0.405845151377397, 0.207784955007898, 0.0])
GKW = np.array([0.022935322010529, 0.063092092629979, 0.104790010322250, 0.140653259715525, 0.169004726639267,
                0.190350578064785, 0.204432940075298, 0.209482141084728])
G7W = {0.949107912342759: 0.129484966168870, 0.741531185599394: 0.279705391489277, 0.405845151377397: 0.381830050505119,
       0.0: 0.417959183673469}


def gk15_nodes(a, b):
    """nodes, Kronrod weights, Gauss-7 weights (0 on the Kronrod-only nodes) on [a, b]."""
    c, h = 0.5*(a + b), 0.5*(b - a)
    xs, wk, wg = [], [], []
    for t, w in zip(GKX, GKW):
        for s in ((-1, 1) if t else (1,)):
            xs.append(c + s*h*t); wk.append(h*w); wg.append(h*G7W.get(t, 0.0))
    o = np.argsort(xs)
    return np.array(xs)[o], np.array(wk)[o], np.array(wg)[o]


# ------------------------------------------------------------------------------------------------ rows
def grid(J, K=2, inject=(), jk=None):
    """the chain grid of spin J (all its rows) at Laplace order K, with the direct points of the listed state Deltas injected (replica
    jk left out if given).  Odd: {D: {m, xs, a, C2, C3H}}; even: {D: {m, xs, a, (kind, h) ...}}."""
    NC.KK = K
    NC.PK = None if J % 2 else K
    if J % 2 == 0:
        NC.KK = 2
    V.JK[0] = jk
    if inject:
        V.INJ[J] = list(inject)
    else:
        V.INJ.pop(J, None)
    try:
        g = V.grid(J)
    finally:
        V.INJ.pop(J, None)
        NC.KK, NC.PK = 2, None
    return g


def direct_x(J, D):
    """x of the direct nodes (state AND reference present) of state D."""
    P = HN.pts('EE', J, D); m = HN.mphys(D, J)
    return np.array(sorted(round(k[1]/m, 6) for k in P if k[0] == D and (DREF, k[1]) in P))


def failed_x(J, D, rtol=0.2):
    """x of direct nodes whose RQMC failed (max relative replica error of the state or reference record > rtol, any kind): the n1
    contour fails as nu -> Delta (x >~ 0.95 on Delta ~ 1000 rows)."""
    m = HN.mphys(D, J); bad = set()
    for kind in ('EE', 'C3H'):
        for (DD, nu), r in HN.pts(kind, J, D).items():
            e = float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300)))
            if e > rtol or not np.all(np.isfinite(r['E'])):
                bad.add(round(nu/m, 6))
    return bad


def kappa(kind, h, J, m, xs):
    env = Q.P4.env if J % 2 else Q.A.env
    SIG = Q.P4.SIG if J % 2 else Q.A.SIG
    return CF.curv_image(kind, h, J, m, np.asarray(xs, float), env, SIG)


def keys_of(J):
    return ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]


def hof(key):
    return key[1] if isinstance(key, tuple) else '+-'


def kof(key):
    return key[0] if isinstance(key, tuple) else key


def subtract_kappa(g, J):
    out = {}
    for D, r in g.items():
        r = dict(r)
        for key in keys_of(J):
            if key in r:
                r[key] = np.asarray(r[key], float) - kappa(kof(key), hof(key), J, r['m'], r['xs'])
        out[D] = r
    return out


def refdiag(J, D, kind='C3H'):
    """per K2 node of row D: (x, |r2/r0|, excess) of the Delta_ref = 2e4 reference record (max over matrix entries); excess = |r2/r0|
    minus its smooth trend c/nu (c = median of |r2/r0| nu over the row).  The trend is the ordinary small-nu growth of the reference's
    Laplace terms (common to state and reference; it cancels in the ratio).  A large excess marks the node where the reference's
    leading partial wave D_h(nu/m_ref) passes through its first zero (HIGHSPIN.md §2): there the K2 series of the RATIO fails."""
    import pickle, glob
    pats = (['w1_p3_%s_J%d*.pkl', 'ee1000_p3_%s_J%d*.pkl', 'q2_p3_%s_J%d_*.pkl', 'p3_pts_%s_J%d.pkl'] if J % 2 else
            ['w1_pts_%s_J%d_k2*.pkl', 'ee1000_pts_%s_J%d_k2*.pkl', 'l9_pts_%s_J%d_k2.pkl', 'q2_pts_%s_J%d_k2*.pkl', 'u5_pts_%s_J%d_k2*.pkl'])
    R = {}
    for p in pats:
        for d in (NC.SFC, os.path.join(NG, 'cache'), HSC):
            for fn in sorted(glob.glob(os.path.join(d, p % (kind, J)))):
                R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items()})
    m = HN.mphys(D, J); out = []
    for (DD, nu), r in R.items():
        if DD != D or (DREF, nu) not in R:
            continue
        rf = R[(DREF, nu)]
        r0 = np.asarray(rf['M0']); a = np.abs(r0).max()
        out.append((nu/m, nu, float(np.abs(np.asarray(rf['M2'])).max()/a)))
    out.sort()
    if not out:
        return []
    c = float(np.median([v*nu for _, nu, v in out]))
    return [(x, v, v - c/nu) for x, nu, v in out]


def ref_zero_x(J, D, h=None):
    """x (on row D) at which the reference's leading partial wave D_h(nu/m_ref) has its first zero."""
    from scipy.optimize import brentq
    h = h or ('+-' if J % 2 else '++')
    f = lambda t: E.Dc_fixed(h, J, t)
    ts = np.linspace(1e-5, 0.3, 6000); v = f(ts)
    i = int(np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0][0])
    return brentq(f, ts[i], ts[i + 1])*HN.mphys(DREF, J)/HN.mphys(D, J)


def refbad(J, D, tol=0.1):
    """x of K2 nodes of row D whose reference series fails (excess |r2/r0| > tol in the EE/C2 or C3H record)."""
    bad = set()
    for kind in (('EE', 'C3H') if True else ()):
        for x, v, ex in refdiag(J, D, kind):
            if ex > tol:
                bad.add(round(x, 6))
    return bad


# ------------------------------------------------------------------------------------------------ shape model (ee1000_chain, copied)
def _dshape(h, J, x):
    D = E.Dc_fixed(h, J, x)
    hh = 1e-6
    Dp = (E.Dc_fixed(h, J, x + hh) - E.Dc_fixed(h, J, x - hh))/(2*hh)
    return D, x*Dp/(2*J + 1)


def _basis(h, J, x, xm, K, x0=0.0):
    D, Dp = _dshape(h, J, x)
    t = (np.clip(x, x0, xm) - x0)/(xm - x0) if x0 > 0 else np.minimum(x, xm)/xm
    P = [t**k for k in range(K + 1)]
    return np.array([p*D for p in P] + [p*Dp for p in P] + [p for p in P]).T


class ShapeFit:
    """Gamma(x) = c(x) D_h + s(x) x D_h'/(2J+1) + r(x) (polynomials of degree K in t = x/x_last, held beyond the last node)."""
    def __init__(self, J, xs, y, h, K=4, x0=0.0):
        self.J, self.K, self.h, self.x0 = J, K, h, x0
        self.xm = float(np.max(xs))
        A = _basis(h, J, np.asarray(xs, float), self.xm, K, x0)
        sc = np.abs(A).max(0); sc[sc == 0] = 1
        c, *_ = np.linalg.lstsq(A/sc, np.asarray(y, float), rcond=None)
        self.c = c/sc
        res = np.asarray(y, float) - A @ self.c
        self.rms_rel = float(np.sqrt(np.mean(res**2))/max(np.sqrt(np.mean(np.asarray(y, float)**2)), 1e-30))
        self.maxres = float(np.abs(res).max())

    def __call__(self, x):
        return _basis(self.h, self.J, np.asarray(x, float), self.xm, self.K, self.x0) @ self.c

    def parts(self, x):
        """(c-part, s-part, r-part) of the fit at x."""
        A = _basis(self.h, self.J, np.asarray(x, float), self.xm, self.K, self.x0); n = self.K + 1
        return A[:, :n] @ self.c[:n], A[:, n:2*n] @ self.c[n:2*n], A[:, 2*n:] @ self.c[2*n:]


def spline_of(xs, y):
    from scipy.interpolate import CubicSpline
    s = CubicSpline(np.asarray(xs, float), np.asarray(y, float), bc_type='natural'); xm = float(np.max(xs))
    return lambda x: s(np.minimum(np.asarray(x, float), xm))


# ------------------------------------------------------------------------------------------------ states
def state_weights(Dg, J, mfl, h, x):
    """per unit Gamma (c2) at x: the k = 2 and k = 3 density factors psi_k(p) env_k(x) [u] / m^2 * dp/dx, p = x m/Lam."""
    Lm = E.lam_of(Dg); m = mfl*Lm; u = 1/mfl**2
    F = Q.NearGap(Dg).F
    p = np.asarray(x, float)*m/Lm
    env = Q.P4.env if J % 2 else Q.A.env
    jac = m/Lm
    w2 = F.psi2(p)*env(h, x, 2)/m**2*jac
    w3 = F.psi3(p)*u*env(h, x, 3)/m**2*jac
    if J % 2 == 0 and h == 'off':
        eo2 = np.sqrt(env('++', x, 2)*env('+-', x, 2)); y = np.asarray(x)**2
        w2 = F.psi2(p)*eo2/m**2*jac; w3 = F.psi3(p)*u*eo2/(2 - y)/m**2*jac
    return w2, w3


def base_of(Dg, J, mfl, npn=256):
    hs = ['+-'] if J % 2 else ['++', '+-']
    st = {h: E.hyb(Dg, J, mfl, h, npn) for h in hs}
    return hs, np.array([st[h]['base'] for h in hs]), np.array([st[h]['Ns'] for h in hs]), st


def xgrid(Dg, mfl, xcut, n_bulk=600):
    """quadrature in x: GL on [0, xcut] (n_bulk nodes), GK15 on [xcut, x_max = Lam/m] (returned separately, with G7 weights)."""
    from numpy.polynomial.legendre import leggauss
    xmax = 1.0/mfl
    t, w = leggauss(n_bulk)
    xb = 0.5*xcut*(t + 1); wb = 0.5*xcut*w
    xe, wk, wg = gk15_nodes(xcut, xmax)
    return xb, wb, xe, wk, wg, xmax


def lam_1x1(base, Ns, d):
    return (base + d)/Ns
