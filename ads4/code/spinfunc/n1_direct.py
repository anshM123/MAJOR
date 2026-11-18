"""
n1_direct (doc 24): NON-PERTURBATIVE evaluation of the CLS-D crystal integrals (no Laplace truncation): the exact 7-variable integrand
(l1_crystal.build for the scalar, l9_points.build_spin for the graviton kinds C2 / C3H) evaluated point-wise (n1_pw) on the real
xi-contour z = z* + L xi through the scattering crystal (L: L^T(-Hess lamF)L = 1, as in l1_crystal.laplace), integrated by randomised
quasi-Monte Carlo (scrambled Sobol, Gaussian importance weight with scale s):

   E[f] := (2 pi)^{-7/2} int d^7 xi exp(lamF + G - lamF_0 - G_0) * pref  =  s^7 E_eta[ exp(lamF + G - lamF_0 - G_0)(s eta) pref e^{|eta|^2/2} ]

so that the full integral is exp(lamF_0 + G_0) (2 pi)^{7/2} det L * E[f]  (the same normalisation in which the Laplace series is
S0 + S2 + S4 (+ S6)).  Error estimate: spread over R independent scramblings (std/sqrt(R)).
Functions:
  scalar(nu, Delta, J, ...)            -> dict(E, err, lap=(S0,S2,S4), logN = lamF0 + G0 + log det L - log beta, wind)
  graviton(kind, nu, Delta, J, ...)    -> dict(E[a, c] (2x2 / 1x1 prefactor integrals), err, Mser (Laplace orders, l6 convention), ...)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.stats import qmc
from scipy.special import ndtri
import n1_pw as PW
import l1_tps as T
import l1_crystal as Cr
import l5_spin as S

CHUNK = int(os.environ.get('N1_CHUNK', '4096'))
WALL = np.exp(2.0)          # exp Re(Phi - Phi_0) above this: past a branch wall (dropped; counted)
WMAX = float(os.environ.get('N1_WMAX', '3.0'))   # per-point winding above this: dropped (counted)
RCUT = float(os.environ.get('N1_RCUT', '5.5'))


def saddle_L(nu, Delta):
    zs, a = Cr.saddle(nu, Delta)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF = Cr.build(zv, nu, Delta, 0)[0]
    _, HF = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = Cr.build(T.linsub(zs, L), nu, Delta, 0)[0]
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    return zs, L


def _etas(m, seed):
    eng = qmc.Sobol(d=7, scramble=True, seed=seed)
    u = eng.random_base2(m)
    u = np.clip(u, 1e-15, 1 - 1e-15)
    return ndtri(u)


def integrate(fun, zs, L, m=14, R=8, s=1.0, seed0=1000):
    """fun(zlist, xi) -> (lamF, G, pref[..., N]) point-wise; returns mean over R scramblings of E[f] (array over pref's leading axes),
    its standard error, and the winding diagnostic."""
    vals = []; wind = 0.0
    for r in range(R):
        eta = _etas(m, seed0 + r)
        acc = None
        for c0 in range(0, eta.shape[0], CHUNK):
            et = eta[c0:c0 + CHUNK]
            xi = s*et
            zl = PW.points(zs, L, xi)
            lamF, G, pref = fun(zl)
            ex = lamF - lamF[0] + G - G[0]
            w = np.exp(ex[1:] + 0.5*np.sum(et*et, axis=1))*s**7
            v = np.sum(pref[..., 1:]*w, axis=-1)
            acc = v if acc is None else acc + v
            wind = max(wind, PW.WIND[0])
        vals.append(acc/eta.shape[0])
    vals = np.array(vals)
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)/np.sqrt(R), wind


def scalar(nu, Delta, J, m=14, R=8, s=1.0, laplace=True):
    zs, L = saddle_L(nu, Delta)
    def fun(zl):
        lamF, G, Q = PW.Cr.build(zl, nu, Delta, J, pref=True)
        return lamF, G, Q
    E, err, wind = integrate(fun, zs, L, m, R, s)
    out = dict(E=E, err=err, wind=wind)
    zv = [T.var(i, zs[i]) for i in range(T.NV)]
    lamF, G, Q = Cr.build(T.linsub(zs, L), nu, Delta, J, pref=True)
    out['logN'] = lamF[0] + G[0] + np.log(np.linalg.det(L)) - Cr.log_beta(Delta, J)
    if laplace:
        out['lap'] = np.array(T.expand_pref(lamF, G, Q))
    return out


def graviton(kind, nu, Delta, J, m=13, R=8, s=1.0, laplace=False):
    """prefactor integrals E[a', c] for the l9 kinds C2 / C3H (conserved-structure basis of l5_spin.conserved).  The physical 2x2 is
    M = E @ inv(beta)^T times exp(logC) (l6 convention), with logC = lamF0 + G0 + log det L."""
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    zs, L = saddle_L(nu, Delta)
    def fun(zl):
        return PW.P9.build_spin(zl, nu, Delta, J, mons, sol, solt, kind)
    E, err, wind = integrate(fun, zs, L, m, R, s)
    zl = PW.points(zs, L, np.zeros((0, 7)))
    lamF, G, pref = PW.P9.build_spin(zl, nu, Delta, J, mons, sol, solt, kind)
    out = dict(E=E, err=err, wind=wind, logC=lamF[0] + G[0] + np.log(np.linalg.det(L)), pref0=pref[..., 0])
    if laplace:
        import l9_points as P9
        out['lap'] = P9.action(nu, Delta, J, kind)
    return out


if __name__ == '__main__':
    # quick test: scalar at large Delta vs Laplace
    for (nu, Dl, J) in ((400.0, 1000.0, 4), (50.0, 100.0, 4), (30.0, 100.0, 10)):
        for (m, s) in ((12, 1.0), (14, 1.0), (14, 1.2)):
            t = time.time()
            r = scalar(nu, Dl, J, m=m, R=8, s=s)
            lap = r['lap']
            print(f"nu={nu} Delta={Dl} J={J} m={m} s={s}: E = {r['E']:.10f} +- {r['err']:.1e}; Laplace S0..S4 cum "
                  f"{lap[0]:.8f} {lap[0]+lap[1]:.8f} {lap.sum():.8f}; |E/lap - 1| = {abs(r['E']/lap.sum() - 1):.2e}; "
                  f"wind {r['wind']:.2f}; {time.time() - t:.1f}s", flush=True)


# ------------------------------------------------------------------ re-centred, log-coordinate contour (doc 24.2)
LOGV = (1, 2, 3, 4)            # u_x, v_x, u_y, v_y -> z_i = c_i exp(w_i), c_i = crystal value (+-i a); others linear (z_i = w_i)


def texp(a):
    """exp of a TPS."""
    d = a.copy(); d[0] = 0.0
    out = T.const(1.0); p = T.const(1.0)
    for n in range(1, T.DEG + 1):
        p = T.mul(p, d)/n
        out = out + p
    return out*np.exp(a[0])


TANH0 = os.environ.get('N1_TANH0', '1') == '1'     # x~0 = tanh(w_0): the light-ray parameter x~0 in (-1, 1) mapped to the real line


def zmap(wv, c, pw=False):
    """w -> z (list), and the log-Jacobian: sum_i w_i (i in LOGV) [+ log sech^2 w_0 if TANH0]."""
    ex = (lambda a: np.exp(a)) if pw else texp
    zv = [c[i]*ex(wv[i]) if i in LOGV else wv[i] for i in range(7)]
    lj = sum(wv[i] for i in LOGV) + np.sum(np.log([c[i] for i in LOGV]))      # dz_i = c_i e^{w_i} dw_i (constant part included)
    if TANH0:
        e2 = ex(2*wv[0])
        if pw:
            zv[0] = (e2 - 1)/(e2 + 1)
            lj = lj + np.log(4*e2/(e2 + 1)**2)
        else:
            zv[0] = T.div(e2 - T.const(1.0), e2 + T.const(1.0))
            lj = lj + T.log(T.div(4*e2, T.mul(e2 + T.const(1.0), e2 + T.const(1.0))))
    return zv, lj


def _zzww(zv, n0=1.0, M=None):
    """Z^2 W^2 of l1_crystal.build (same geometry), for backend module M (default: TPS)."""
    M = M or Cr
    TT = M.T
    xt0, ux, vx, uy, vy, xp, yp = zv
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = M.dot(x, x); rx = TT.recip(-xsq)
    x4 = M.scal(rx, x); x2 = M.neg(x4); x1 = M.neg(y); x3 = y
    e = M.cpt(1.0, 1.0, 0.0)
    def msq(p, q):
        r = M.sub(p, q); return -M.dot(r, r)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    x15, x25, x53, x54 = M.sub(x1, e), M.sub(x2, e), M.sub(e, x3), M.sub(e, x4)
    Z = M.sub(M.scal(m25, x15), M.scal(m15, x25)); W = M.sub(M.scal(m45, x53), M.scal(m35, x54))
    return TT.mul(M.dot(Z, Z), M.dot(W, W))


def phi_w(wv, c, nu, Delta, J, n0=1.0):
    """Phi(w) = lamF + G_scal(pref form) + (J/2) log(Z^2 W^2) + log-Jacobian: the smooth, J-growing part of the scalar integrand in the
    w-variables (P_J(Z.W/|Z||W|) excluded: it stays a bounded prefactor).  TPS backend."""
    zv, lj = zmap(wv, c)
    lamF, G, Q = Cr.build(zv, nu, Delta, J, n0, pref=True)
    P = lamF + G + lj
    if J > 0:
        P = P + (J/2)*T.log(_zzww(zv, n0))
    return P


def full_saddle(nu, Delta, J, it=40, tol=1e-12, verbose=False, Jc=None):
    """Newton for grad_w Phi = 0 from the crystal (w = 0); returns (c, wc, L): L^T(-Hess_w Phi)L = 1 at wc.  Jc: spin used for the
    centring (default J)."""
    Jc = J if Jc is None else Jc
    c, _ = Cr.saddle(nu, Delta)
    wc = np.zeros(7, complex)
    for k in range(it):
        wv = [T.var(i, wc[i]) for i in range(T.NV)]
        g, H = T.grad_hess(phi_w(wv, c, nu, Delta, Jc))
        dw = np.linalg.solve(H, g)
        wc = wc - dw
        if verbose:
            print(k, np.abs(g).max(), np.abs(dw).max())
        if np.abs(dw).max() < tol:
            break
    wv = [T.var(i, wc[i]) for i in range(T.NV)]
    _, H = T.grad_hess(phi_w(wv, c, nu, Delta, Jc))
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    return c, wc, L


def integrate_w(fun, c, wc, L, m=14, R=8, s=1.0, seed0=1000, cut=None):
    """as integrate(), on the w-contour w = wc + L xi; fun(zlist) -> (lamF, G, pref); log-Jacobian added here.
    cut: optional radius in eta beyond which points are dropped (diagnostic)."""
    vals = []; wind = 0.0
    for r in range(R):
        eta = _etas(m, seed0 + r)
        acc = None
        for c0 in range(0, eta.shape[0], CHUNK):
            et = eta[c0:c0 + CHUNK]
            xi = s*et
            X = np.vstack([np.zeros((1, 7)), xi])
            Wp = wc[None, :] + X @ L.T
            PW.set_batch(Wp.shape[0])
            zl, lj = zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
            lamF, G, pref = fun(zl)
            ex = lamF - lamF[0] + G - G[0] + lj - lj[0]
            w = np.exp(ex[1:] + 0.5*np.sum(et*et, axis=1))*s**7
            if cut is not None:
                w = np.where(np.sum(xi*xi, axis=1) <= cut**2, w, 0.0)
            w = np.where(np.isfinite(w) & (np.exp(ex[1:].real) < WALL), w, 0.0)   # wall points (smooth part above its saddle value)
            v = np.sum(pref[..., 1:]*w, axis=-1)
            acc = v if acc is None else acc + v
            wind = max(wind, PW.WIND[0])
        vals.append(acc/eta.shape[0])
    vals = np.array(vals)
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)/np.sqrt(R), wind


def norm_w(fun, c, wc, L):
    """log of the normalisation exp(lamF + G + logJac)(wc) (2 pi)^{7/2} det L, and the prefactor at wc."""
    PW.set_batch(1)
    zl, lj = zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamF, G, pref = fun(zl)
    return lamF[0] + G[0] + lj[0] + 3.5*np.log(2*np.pi) + np.log(np.linalg.det(L)), pref[..., 0]


def gh_nodes(ns):
    """tensor Gauss-Hermite (probabilists') nodes/weights in 7D, ns: points per dimension (int or 7-list)."""
    ns = [ns]*7 if np.isscalar(ns) else list(ns)
    xs, ws = [], []
    for n in ns:
        x, w = np.polynomial.hermite_e.hermegauss(n)
        xs.append(x); ws.append(w/np.sqrt(2*np.pi))
    G = np.meshgrid(*xs, indexing='ij'); Wt = np.meshgrid(*ws, indexing='ij')
    X = np.stack([g.ravel() for g in G], axis=1); Wv = np.prod(np.stack([w.ravel() for w in Wt], axis=1), axis=1)
    return X, Wv


def integrate_gh(fun, c, wc, L, ns=6):
    """E[f] (same normalisation as integrate_w) by a tensor Gauss-Hermite rule on the w-contour."""
    X, Wv = gh_nodes(ns)
    acc = None; wind = 0.0
    for c0 in range(0, X.shape[0], CHUNK):
        xi = X[c0:c0 + CHUNK]
        Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), xi]) @ L.T
        PW.set_batch(Wp.shape[0])
        zl, lj = zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
        lamF, G, pref = fun(zl)
        ex = lamF - lamF[0] + G - G[0] + lj - lj[0]
        w = np.exp(ex[1:] + 0.5*np.sum(xi*xi, axis=1))*Wv[c0:c0 + CHUNK]
        v = np.sum(pref[..., 1:]*w, axis=-1)
        acc = v if acc is None else acc + v
        wind = max(wind, PW.WIND[0])
    return acc, wind


def smolyak(level, d=7, rule=lambda i: 2*i - 1):
    """Smolyak (combination technique) Gauss-Hermite rule (probabilists', N(0,1)) of the given level >= 1 in d dims:
    Q = sum_{max(1,level-d+1) <= |i|-d+1 <= level} (-1)^(level-q) C(d-1, level-q) (x) Q_{i_k}, q = |i| - d + 1.  Returns (X, W) merged."""
    from itertools import product
    from math import comb
    rules = {}
    def r1(i):
        if i not in rules:
            x, w = np.polynomial.hermite_e.hermegauss(rule(i)); rules[i] = (x, w/np.sqrt(2*np.pi))
        return rules[i]
    pts = {}
    def comps(total, k):
        if k == 1:
            yield (total,); return
        for a in range(1, total - k + 2):
            for rest in comps(total - a, k - 1):
                yield (a,) + rest
    for q in range(max(1, level - d + 1), level + 1):
        coef = (-1)**(level - q)*comb(d - 1, level - q)
        for idx in comps(q + d - 1, d):
            xs = [r1(i)[0] for i in idx]; ws = [r1(i)[1] for i in idx]
            for tup in product(*[range(len(x)) for x in xs]):
                key = tuple(round(xs[k][tup[k]], 12) for k in range(d))
                pts[key] = pts.get(key, 0.0) + coef*np.prod([ws[k][tup[k]] for k in range(d)])
    X = np.array(list(pts.keys())); W = np.array(list(pts.values()))
    return X, W


def integrate_rule(fun, c, wc, L, X, Wv, cut=None):
    """E[f] with a given rule (X: (n,7) xi-nodes, Wv weights for N(0,1)); wall/cut masks as integrate_w."""
    acc = None; wind = 0.0; nmask = 0
    for c0 in range(0, X.shape[0], CHUNK):
        xi = X[c0:c0 + CHUNK]
        Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), xi]) @ L.T
        PW.set_batch(Wp.shape[0])
        zl, lj = zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
        lamF, G, pref = fun(zl)
        ex = lamF - lamF[0] + G - G[0] + lj - lj[0]
        w = np.exp(ex[1:] + 0.5*np.sum(xi*xi, axis=1))*Wv[c0:c0 + CHUNK]
        ok = np.isfinite(w) & (np.exp(ex[1:].real) < WALL)
        if cut is not None:
            ok &= np.sum(xi*xi, axis=1) <= cut**2
        nmask += int(np.sum(~ok))
        w = np.where(ok, w, 0.0)
        v = np.sum(pref[..., 1:]*w, axis=-1)
        acc = v if acc is None else acc + v
        wind = max(wind, PW.WIND[0])
    return acc, wind, nmask


def evaluate(fun, c, wc, L, xi, nu, Delta, J):
    """point values on the w-contour: returns (val[..., n] = pref * exp(lamF + G + logJ - (same at wc)), ok-mask, Phi-excess).
    Masks: Re(Phi - Phi_0) < log WALL (Phi = scalar smooth part incl. (J/2) log Z^2W^2 and log-Jacobian) and winding < WMAX."""
    nl = len(LEVELS)
    if nl > 1:
        windacc = None
        for k, t in enumerate(LEVELS[:-1]):
            PW.path_mode('record' if k == 0 else 'replay_record')
            _evaluate1(fun, c, wc, L, t*xi, nu, Delta, J)
            windacc = PW.WINDP[1:].copy() if windacc is None else np.maximum(windacc, PW.WINDP[1:])
        PW.path_mode('replay')
        try:
            val, ok, dph = _evaluate1(fun, c, wc, L, xi, nu, Delta, J)
        finally:
            PW.path_mode('plain')
        return val, ok, dph
    return _evaluate1(fun, c, wc, L, xi, nu, Delta, J)


LEVELS = tuple(float(v) for v in os.environ.get('N1_LEVELS', '0.5,1').split(','))


def _evaluate1(fun, c, wc, L, xi, nu, Delta, J):
    Wp = wc[None, :] + np.vstack([np.zeros((1, 7)), xi]) @ L.T
    PW.set_batch(Wp.shape[0])
    zl, lj = zmap([Wp[:, i].copy() for i in range(7)], c, pw=True)
    lamS, GS, QS = PW.Cr.build(zl, nu, Delta, J, pref=True)
    ph = lamS + GS + lj
    if J > 0:
        ph = ph + (J/2)*PW.log(_zzww(zl, M=PW.Cr))
    dph = (ph - ph[0]).real[1:]
    wind1 = PW.WINDP[1:].copy()
    if fun is None:
        lamF, G, pref = lamS, GS, QS
    else:
        lamF, G, pref = fun(zl)
    wind = np.maximum(wind1, PW.WINDP[1:])
    ex = lamF - lamF[0] + G - G[0] + lj - lj[0]
    ok = (dph < np.log(WALL)) & (wind < WMAX)
    val = pref[..., 1:]*np.exp(ex[1:])
    ok &= np.all(np.isfinite(val.reshape(-1, val.shape[-1])), axis=0)
    return val, ok, dph


def integrate_rule2(fun, c, wc, L, X, Wv, nu, Delta, J, cut=None):
    """E[f] = sum_k W_k f(xi_k) e^{|xi_k|^2/2} with masks (evaluate); returns (E, n_masked, masked |weight| fraction)."""
    acc = None; nm = 0; mw = 0.0; tw = 0.0
    for c0 in range(0, X.shape[0], CHUNK):
        xi = X[c0:c0 + CHUNK]
        val, ok, _ = evaluate(fun, c, wc, L, xi, nu, Delta, J)
        if cut is not None:
            ok &= np.sum(xi*xi, axis=1) <= cut**2
        w = Wv[c0:c0 + CHUNK]*np.exp(0.5*np.sum(xi*xi, axis=1))
        nm += int(np.sum(~ok)); mw += float(np.sum(np.abs(Wv[c0:c0 + CHUNK])[~ok])); tw += float(np.sum(np.abs(Wv[c0:c0 + CHUNK])))
        v = np.sum(np.where(ok, val, 0.0)*w, axis=-1)
        acc = v if acc is None else acc + v
    return acc, nm, mw/tw


def integrate_qmc(fun, c, wc, L, nu, Delta, J, m=14, R=8, s=1.2, cut=None, seed0=1000, anti=True):
    """RQMC version of integrate_rule2 (Gaussian proposal of scale s); anti: antithetic pairs +-eta (odd terms cancel pairwise);
    returns (mean, stderr over R scramblings, n_masked)."""
    vals = []; nm = 0
    for r in range(R):
        eta = _etas(m, seed0 + r)
        if anti:
            eta = np.concatenate([eta, -eta])
        acc = None
        for c0 in range(0, eta.shape[0], CHUNK):
            et = eta[c0:c0 + CHUNK]; xi = s*et
            val, ok, _ = evaluate(fun, c, wc, L, xi, nu, Delta, J)
            if cut is not None:
                ok &= np.sum(xi*xi, axis=1) <= cut**2
            nm += int(np.sum(~ok))
            w = np.exp(0.5*np.sum(et*et, axis=1))*s**7
            v = np.sum(np.where(ok, val, 0.0)*w, axis=-1)
            acc = v if acc is None else acc + v
        vals.append(acc/eta.shape[0])
    vals = np.array(vals)
    LAST_VALS[0] = vals
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)/np.sqrt(R), nm


LAST_VALS = [None]


# ------------------------------------------------------------------ adaptive Gaussian proposal (doc 24.2)
def adapt(fun, c, wc, L, nu, Delta, J, n=8192, its=3, infl=1.2, seed=5):
    """moment-matched proposal N(mu, C) in xi for |f| (f = integrand without the Gaussian), from pilot samples; returns (mu, C)."""
    rng = np.random.default_rng(seed)
    mu = np.zeros(7); C = np.eye(7)
    for it in range(its):
        Lc = np.linalg.cholesky(C); x = mu + rng.normal(size=(n, 7)) @ Lc.T
        val, ok, _ = evaluate(fun, c, wc, L, x, nu, Delta, J)
        v = val if val.ndim == 1 else np.abs(val.reshape(-1, val.shape[-1])).sum(axis=0)
        d = x - mu
        lq = -0.5*np.einsum('ni,ij,nj->n', d, np.linalg.inv(C), d) - 0.5*np.log(np.linalg.det(C))
        a = np.where(ok, np.abs(v)*np.exp(-lq), 0.0)
        a = a/a.sum()
        mu = (a[:, None]*x).sum(0); d = x - mu
        C = infl*(a[:, None, None]*d[:, :, None]*d[:, None, :]).sum(0)
    return mu, C


def integrate_ad(fun, c, wc, L, nu, Delta, J, mu, C, m=13, R=8, seed0=2000, anti=True):
    """RQMC with the proposal N(mu, C): E[f] = (2pi)^{-7/2} int val dxi = E_q[val (2pi)^{-7/2}/q]; returns (mean, stderr, n_masked)."""
    Lc = np.linalg.cholesky(C); Ci = np.linalg.inv(C); ldet = 0.5*np.log(np.linalg.det(C))
    vals = []; nm = 0
    for r in range(R):
        eta = _etas(m, seed0 + r)
        if anti:
            eta = np.concatenate([eta, -eta])
        acc = None
        for c0 in range(0, eta.shape[0], CHUNK):
            et = eta[c0:c0 + CHUNK]; xi = mu + et @ Lc.T
            val, ok, _ = evaluate(fun, c, wc, L, xi, nu, Delta, J)
            nm += int(np.sum(~ok))
            w = np.exp(0.5*np.sum(et*et, axis=1) + ldet)
            v = np.sum(np.where(ok, val, 0.0)*w, axis=-1)
            acc = v if acc is None else acc + v
        vals.append(acc/eta.shape[0])
    vals = np.array(vals)
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)/np.sqrt(R), nm
