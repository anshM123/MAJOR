"""
n3_flow (doc 25): Lefschetz-thimble (steepest-descent) evaluation of the exact 7-variable CLS-D crystal integral at SMALL nu (5-20),
where the real-xi contour of n1_direct crosses branch walls (doc 24.2, 24.8).

Contour.  Same variables as n1_direct (w: log coordinates for u_x, v_x, u_y, v_y; x~0 = tanh w_0; others linear), saddle w_c of
   Phi = lamF + G_scal(pref form) + logJac + (J/2) log(Z^2 W^2)        (the smooth, J-growing scalar exponent; n1_direct.phi_w),
xi-coordinates w = w_c + L xi with L^T(-Hess Phi)L = 1, so Phi = Phi_c - xi.xi/2 + O(xi^3).  The contour is the (numerically integrated)
downward gradient flow of Phi out of the saddle -- the Lefschetz thimble of the crystal saddle:
   d xi/dt = -conj(dPhi/dxi),   xi(0) = r0 n,  n in S^6,   point (n, rho) <-> xi(n, t = log(rho/r0)),
which reduces to xi = rho n near the saddle (the n1 contour there).  Along the thimble Im Phi is constant and Re Phi decreases
monotonically, so the integrand never passes a branch wall: logs and fractional powers are continued along the flow lines from the
saddle (n1_pw replay mode, step by step; winding per step recorded).
Measure.  d^7xi = det[d xi/dt, V_1, ..., V_6] dt dOmega_6, with the tangent vectors V_a (V_a(0) = r0 e_a, (n, e_a) an oriented
orthonormal frame) propagated by the linearised flow dV/dt = -conj(H V) (H = Hess Phi); (xi, V) integrated jointly by RK4, so the
Jacobian is that of the contour actually used (up to the FD error of H).  Gradient and Hessian by 4-point holomorphic stencils
(f'(z) = sum_w w^-1 f(z + h w)/(4h), f'' = sum_w w^-2 f(z + h w)/(2h^2), w in {1, i, -1, -i}; O(h^4)).
Integral.  E = (2 pi)^{-7/2} int_contour d^7xi exp(Pint - Pint_c) pref  (Pint = lamF + G + logJac: n1 normalisation, so that
log Y = logN - log beta + log E as in n1_scalar_validate), = (2 pi)^{-7/2} |S^6| < int dt det(...) exp(...) pref >_n ;
the radial (t) integral by the trapezoid rule on the RK4 grid, the angular average by randomised QMC on S^6 (antithetic +-n),
error from R independent scramblings.
Functions:
  saddle(nu, Delta, J)            -> c, wc, L   (damped Newton with continuation in nu from nu >= 20 where the crystal start converges)
  flow(...)                       -> per-direction paths (xi_k, Phi_k, pref_k, det_k)
  integrate(nu, Delta, J, kind)   -> dict(E, err, ...); kind 'S' (scalar) or graviton kinds 'C2', 'C3H' (pref matrices; graviton
                                     integrand evaluated on the scalar thimble, logs continued along the stored flow lines)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import gamma, pi
from scipy.stats import qmc
from scipy.special import ndtri
import n1_direct as D
import n1_pw as PW
import l1_tps as T
import l1_crystal as Cr

H_FD = float(os.environ.get('N3_H', '0.01'))
R0 = float(os.environ.get('N3_R0', '0.02'))
DT = float(os.environ.get('N3_DT', '0.1'))
SMAX = float(os.environ.get('N3_SMAX', '20'))
WCUT = float(os.environ.get('N3_WCUT', '1.0'))          # per-evaluation winding (rad) above which a line is cut (branch safety)          # stop a flow line once Re(Pint - Pint_c) < -SMAX (and log|pref/pref_c| accounted)
TMAX = float(os.environ.get('N3_TMAX', '9.0'))
AREA6 = 16*pi**3/15                                     # |S^6|
OM = np.array([1, 1j, -1, -1j])


# ------------------------------------------------------------------ saddle with continuation in nu
def _newton(nu, Delta, J, w0, it=80, cap=0.15):
    c, _ = Cr.saddle(nu, Delta)
    wc = w0.copy()
    for k in range(it):
        wv = [T.var(i, wc[i]) for i in range(T.NV)]
        g, H = T.grad_hess(D.phi_w(wv, c, nu, Delta, J))
        dw = np.linalg.solve(H, g)
        s = 1.0 if np.abs(dw).max() <= cap else cap/np.abs(dw).max()
        wc = wc - s*dw
        if not np.all(np.isfinite(wc)):
            raise RuntimeError('newton diverged')
        if np.abs(dw).max() < 1e-12:
            break
    wv = [T.var(i, wc[i]) for i in range(T.NV)]
    g, H = T.grad_hess(D.phi_w(wv, c, nu, Delta, J))
    return c, wc, H, np.abs(g).max()


_SC = {}


def saddle(nu, Delta, J, nu_start=20.0, step=1.5):
    """saddle of Phi in w (continued from nu_start downwards if nu < nu_start); returns c, wc, L, |grad| at wc.
    Complex nu (flow function at complex nu_f, doc 25): path 20 -> Re nu + 5 (real) -> + i Im nu -> nu, avoiding the real fold nu*."""
    key = (complex(nu), float(Delta), int(J))
    if key in _SC:
        return _SC[key]
    if np.iscomplexobj(nu) and abs(np.imag(nu)) > 0:
        nr = float(np.real(nu)); ni = float(np.imag(nu)); top = max(nu_start, nr + 5.0)
        nus = list(np.arange(top, nr + 5.0 - 1e-9, -step)) + [nr + 5.0 + 1j*ni*f for f in np.linspace(0.25, 1, 4)] + \
              [nr + 5.0 - 5.0*f + 1j*ni for f in np.linspace(0.1, 1, 10)]
    else:
        nu = float(np.real(nu))
        nus = [nu] if nu >= nu_start else list(np.arange(nu_start, nu, -step)) + [nu]
    w = np.zeros(7, complex)
    for v in nus:
        c, w, H, gm = _newton(v if np.iscomplexobj(v) and np.imag(v) != 0 else float(np.real(v)), Delta, J, w)
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    _SC[key] = (c, w, L, gm)
    return _SC[key]


# ------------------------------------------------------------------ point evaluation with branch replay
class Evaluator:
    """Phi / Pint / pref on batches of xi with logs continued from per-direction references."""

    def __init__(self, nu, Delta, J, c, wc, L, kind='S', nuf=None):
        self.nu, self.Delta, self.J, self.c, self.wc, self.L, self.kind = nu, Delta, J, c, wc, L, kind
        self.nuf = nu if nuf is None else nuf             # nu of the flow function (contour); nu: of the integrand
        if kind != 'S':
            import l5_spin as S
            self.mons, self.sol, _ = S.conserved(Delta, J)
            _, self.solt, _ = S.conserved(3.0 - Delta, J)

    def _z(self, xi):
        Wp = self.wc[None, :] + xi @ self.L.T
        PW.set_batch(Wp.shape[0])
        return D.zmap([Wp[:, i].copy() for i in range(7)], self.c, pw=True)

    def phi(self, xi):
        """scalar exponent pieces: Phi (flow function), Pint (integrand exponent), Q (scalar prefactor)."""
        zl, lj = self._z(xi)
        lamS, GS, QS = PW.Cr.build(zl, self.nuf, self.Delta, self.J, pref=True)
        Pint = lamS + GS + lj
        Phi = Pint + ((self.J/2)*PW.log(D._zzww(zl, M=PW.Cr)) if self.J > 0 else 0.0)
        return Phi, Pint, QS

    def pint(self, xi):
        """the integrand exponent and prefactor at the integrand's nu (same log sequence as phi up to the last log)."""
        zl, lj = self._z(xi)
        lamS, GS, QS = PW.Cr.build(zl, self.nu, self.Delta, self.J, pref=True)
        return lamS + GS + lj, QS

    def grav(self, xi):
        zl, lj = self._z(xi)
        lamF, G, pref = PW.P9.build_spin(zl, self.nu, self.Delta, self.J, self.mons, self.sol, self.solt, self.kind)
        return lamF + G + lj, pref

    def run(self, fn, xi, refs, rep=1):
        """evaluate fn(xi) with the log references refs (list of arrays (..., Nbase)), each repeated rep times along the batch axis
        (batch = Nbase x rep, base-major).  refs None: plain mode relative to xi[0] (must be the saddle) -- or, if a branch at the saddle
        was set by set_branch (saddles away from the crystal), replay from those references.  Returns (result, new refs of all batch
        points, per-point winding)."""
        sb = getattr(self, 'saddle_refs', {}).get(fn.__name__)
        if refs is None and sb is not None:
            refs = [np.repeat(r, xi.shape[0], axis=-1) for r in sb]; rep = 1
        if refs is None:
            PW.path_mode('record')
        else:
            PW.REC[1] = [np.repeat(r, rep, axis=-1) for r in refs]
            PW.path_mode('replay_record')
        try:
            out = fn(xi)
            new = PW.REC[1]
            wind = PW.WINDP.copy()
        finally:
            PW.path_mode('plain')
        return out, new, wind


W_REF = 1e-4*np.array([0.31 + 0.17j, -0.23 + 0.41j, 0.37 - 0.29j, -0.13 - 0.47j, 0.43 + 0.11j, -0.19 + 0.07j, 0.05 - 0.21j])
# generic point next to the crystal (w = 0): at the crystal itself some log arguments are exactly on the negative real axis (symmetric
# configuration), so 'principal log there' is decided by round-off; W_REF fixes the convention (the same for every state and contour).


def set_branch(ev, w_from=None, nstep=200, fns=('phi', 'pint', 'grav')):
    """log branch at ev's saddle by continuation along the straight w-path from w_from (principal logs there: the n1 convention at the
    symmetric crystal saddle, nu >= 15) to ev.wc; stored in ev.saddle_refs (used by run(..., refs=None)).  Returns the max winding."""
    if w_from is None:                       # W_REF next to the crystal of the INTEGRAND's (real) nu, in this evaluator's w-coordinates
        cr, _ = Cr.saddle(float(np.real(ev.nu)), ev.Delta)
        w_from = W_REF.copy()
        for i in (1, 2, 3, 4):
            w_from[i] += np.log(cr[i]/ev.c[i])
    Li = np.linalg.inv(ev.L)
    ev.saddle_refs = {}
    wmax = 0.0
    for name in fns:
        if name == 'grav' and ev.kind == 'S':
            continue
        fn = getattr(ev, name)
        xs = [Li @ ((1 - s)*(w_from - ev.wc)) for s in np.linspace(0, 1, nstep + 1)]
        _, refs, w = ev.run(fn, xs[0][None, :].astype(complex), None)
        for x in xs[1:]:
            _, refs, w = ev.run(fn, x[None, :].astype(complex), refs, rep=1)
            wmax = max(wmax, float(w.max()))
        ev.saddle_refs[name] = refs
    return wmax


def _sub(refs, idx):
    return [r[..., idx] for r in refs]


# ------------------------------------------------------------------ the stencil
def _dirs():
    ds = [np.eye(7)[i] for i in range(7)]
    pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    ds += [np.eye(7)[i] + np.eye(7)[j] for (i, j) in pairs]
    return np.array(ds), pairs


DIRS, PAIRS = _dirs()                       # 28 directions
NS = 1 + 4*len(DIRS)                        # 113 stencil points (centre first)
STEN = np.vstack([np.zeros(7)] + [H_FD*w*d for d in DIRS for w in OM]).astype(complex)


def grad_hess(Fv):
    """Fv: (Nb, NS) values of Phi on the stencil -> g (Nb, 7), H (Nb, 7, 7)."""
    F = Fv[:, 1:].reshape(Fv.shape[0], len(DIRS), 4)
    d1 = (F*(1/OM)[None, None, :]).sum(-1)/(4*H_FD)
    d2 = (F*(1/OM**2)[None, None, :]).sum(-1)/(2*H_FD**2)
    g = d1[:, :7]
    Hd = d2[:, :7]
    H = np.zeros((Fv.shape[0], 7, 7), complex)
    for i in range(7):
        H[:, i, i] = Hd[:, i]
    for k, (i, j) in enumerate(PAIRS):
        hij = 0.5*(d2[:, 7 + k] - Hd[:, i] - Hd[:, j])
        H[:, i, j] = hij; H[:, j, i] = hij
    return g, H


def frames(n):
    """oriented orthonormal frames (n, e_1..e_6): returns E (Nd, 7, 7) with E[:, :, 0] = n, det = +1."""
    Nd = n.shape[0]
    E = np.zeros((Nd, 7, 7))
    rng = np.random.default_rng(12345)
    for d in range(Nd):
        M = np.c_[n[d], rng.normal(size=(7, 6))]
        Q, R = np.linalg.qr(M)
        if R[0, 0] < 0:
            Q[:, 0] *= -1
        Q[:, 0] = n[d]
        if np.linalg.det(Q) < 0:
            Q[:, 1] *= -1
        E[d] = Q
    return E


# ------------------------------------------------------------------ the flow
GSTEP = float(os.environ.get('N3_GSTEP', '0.08'))     # max step length |d xi| per RK4 substep
HSTEP = float(os.environ.get('N3_HSTEP', '0.4'))      # max h ||H|| per substep
NSUBMAX = int(os.environ.get('N3_NSUBMAX', '64'))


def flow(ev, n, dt=DT, tmax=TMAX, smax=SMAX, verbose=False):
    """flow lines from r0 n (n: (Nd, 7) unit vectors), output on the uniform grid t_k = k dt (spectral trapezoid in t), each output
    interval integrated by n_sub RK4 substeps per line (n_sub from |grad| and ||Hess|| at the start of the interval).  Returns dict
    with the per-output-step records (t, idx, xi, Pint, Phi, Q, det) and diagnostics."""
    Nd = n.shape[0]
    (Phi0, Pint0, Q0), ref0, _ = ev.run(ev.phi, np.zeros((1, 7), complex), None)
    Phic, Pintc, Qc = Phi0[0], Pint0[0], Q0[0]
    if ev.nu != ev.nuf:
        (Pi0, Qi0), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
        Pintc, Qc = Pi0[0], Qi0[0]
    refs = [np.repeat(r, Nd, axis=-1) for r in ref0]
    E = frames(n)
    xi = R0*n.astype(complex)
    V = R0*E[:, :, 1:].astype(complex)
    act = np.ones(Nd, bool); sbad = np.zeros(Nd, bool)
    out = dict(t=[], Phic=Phic, Pintc=Pintc, Qc=Qc, frames=E, nbad=0, nsub_max=0, wind=0.0, imdrift=0.0)

    def field(xs, Vs, idx, upd):
        pts = (xs[:, None, :] + STEN[None, :, :]).reshape(-1, 7)
        (Phi, Pint, Q), new, wind = ev.run(ev.phi, pts, _sub(refs, idx), rep=NS)
        Phi = Phi.reshape(-1, NS); g, H = grad_hess(Phi)
        if upd:
            cen = np.arange(0, pts.shape[0], NS)
            for a, b in zip(refs, new):
                a[..., idx] = b[..., cen]
        return (-np.conj(g), -np.conj(np.einsum('nij,njk->nik', H, Vs)), Phi[:, 0], Pint.reshape(-1, NS)[:, 0],
                Q.reshape(-1, NS)[:, 0], wind.reshape(-1, NS).max(axis=1), g, H)

    t = 0.0
    while t < tmax + 1e-9 and act.any():
        idx = np.where(act)[0]
        k1x, k1v, Phi, Pint, Q, w1, g, H = field(xi[idx], V[idx], idx, True)
        if ev.nu != ev.nuf:
            if not hasattr(ev, 'nref'):
                ev.nref = len(ev.run(ev.pint, np.zeros((1, 7), complex), None)[1])
            nref = ev.nref
            (Pint, Q), _, w1b = ev.run(ev.pint, xi[idx], [r[..., idx] for r in refs[:nref]], rep=1)
            w1 = np.maximum(w1, w1b)
        detM = np.linalg.det(np.concatenate([k1x[:, :, None], V[idx]], axis=2))
        out['t'].append(dict(t=t, idx=idx, xi=xi[idx].copy(), Pint=Pint, Phi=Phi, Q=Q, det=detM, wind=w1))
        out['wind'] = max(out['wind'], float(w1.max()))
        out['imdrift'] = max(out['imdrift'], float(np.abs((Phi - Phic).imag).max()))
        bad = ~np.isfinite(k1x).all(axis=1) | ~np.isfinite(k1v).all(axis=(1, 2)) | (w1 > WCUT) | sbad[idx]
        out['t'][-1]['ok'] = ~bad
        out['lost'] = max(out.get('lost', -1e9), float(np.max(np.where(bad, (Phi - Phic).real, -1e9))))
        stop = (Phi - Phic).real < -smax
        if verbose:
            print(f"   t={t:.2f}: active {idx.size}, Re dPhi min/med {np.min((Phi - Phic).real):.1f}/{np.median((Phi - Phic).real):.1f}, "
                  f"|Im dPhi| max {np.abs((Phi - Phic).imag).max():.1e}, wind max {w1.max():.2f}, |g| med/max {np.median(np.linalg.norm(g, axis=1)):.2f}/"
                  f"{np.max(np.linalg.norm(g, axis=1)):.2f}", flush=True)
        out['nbad'] += int(bad.sum())
        keep = ~(stop | bad)
        act[idx[~keep]] = False
        idx = idx[keep]
        if idx.size == 0:
            break
        k1x, k1v, g, H = k1x[keep], k1v[keep], g[keep], H[keep]
        gn = np.linalg.norm(g, axis=1); hn = np.linalg.norm(H.reshape(H.shape[0], -1), axis=1)
        nsub = np.clip(np.ceil(np.maximum(dt*gn/GSTEP, dt*hn/HSTEP)), 1, NSUBMAX).astype(int)
        out['nsub_max'] = max(out['nsub_max'], int(nsub.max()))
        for j in range(int(nsub.max())):
            s = nsub > j
            ii = idx[s]
            h = (dt/nsub[s])[:, None]
            x0, v0 = xi[ii], V[ii]
            if j == 0:
                a1x, a1v = k1x[s], k1v[s]
            else:
                r1 = field(x0, v0, ii, True); a1x, a1v = r1[:2]; sbad[ii] |= r1[5] > WCUT
            r2 = field(x0 + 0.5*h*a1x, v0 + 0.5*h[:, :, None]*a1v, ii, False); a2x, a2v = r2[:2]
            r3 = field(x0 + 0.5*h*a2x, v0 + 0.5*h[:, :, None]*a2v, ii, False); a3x, a3v = r3[:2]
            r4 = field(x0 + h*a3x, v0 + h[:, :, None]*a3v, ii, False); a4x, a4v = r4[:2]
            sbad[ii] |= (r2[5] > WCUT) | (r3[5] > WCUT) | (r4[5] > WCUT)
            xi[ii] = x0 + h/6*(a1x + 2*a2x + 2*a3x + a4x)
            V[ii] = v0 + h[:, :, None]/6*(a1v + 2*a2v + 2*a3v + a4v)
        t += dt
    out['dt'] = dt
    return out


def radial_scalar(fl):
    """per-direction t-integral of det * exp(Pint - Pint_c) * Q (trapezoid; the t = 0 end is O(r0^7))."""
    Nd = fl['frames'].shape[0]
    acc = np.zeros(Nd, complex)
    for rec in fl['t']:
        f = rec['det']*np.exp(rec['Pint'] - fl['Pintc'])*rec['Q']
        acc[rec['idx']] += np.where(np.isfinite(f) & rec['ok'], f, 0.0)*fl['dt']
    return acc


def radial_grav(ev, fl):
    """graviton integrand on the same flow lines: logs continued step by step along each line from the saddle."""
    Nd = fl['frames'].shape[0]
    (Pg0, pref0), ref0, _ = ev.run(ev.grav, np.zeros((1, 7), complex), None)
    Pgc = Pg0[0]
    refs = [np.repeat(r, Nd, axis=-1) for r in ref0]
    acc = None; wmax = 0.0
    for rec in fl['t']:
        idx = rec['idx']
        (Pg, pref), new, wind = ev.run(ev.grav, rec['xi'], _sub(refs, idx), rep=1)
        for a, b in zip(refs, new):
            a[..., idx] = b
        wmax = max(wmax, float(wind.max()))
        f = rec['det']*np.exp(Pg - Pgc)*pref                     # (na, nc, n)
        f = np.where(np.isfinite(f) & rec['ok'] & (wind < WCUT), f, 0.0)
        if acc is None:
            acc = np.zeros(f.shape[:-1] + (Nd,), complex)
        acc[..., idx] += f*fl['dt']
    return acc, Pgc, pref0[..., 0], wmax


def sphere_pts(m, seed):
    eng = qmc.Sobol(d=7, scramble=True, seed=seed)
    u = np.clip(eng.random_base2(m), 1e-12, 1 - 1e-12)
    g = ndtri(u)
    return g/np.linalg.norm(g, axis=1, keepdims=True)


NUF = float(os.environ.get('N3_NUF', '15'))            # flow (contour) nu: max(nu, NUF) -- avoids the soft-block pitchfork (nu* ~ 11-12)


def integrate(nu, Delta, J, kind='S', m=7, R=4, seed0=777, anti=True, verbose=False, sad=None, nuf=None):
    """E (same normalisation as n1_direct.integrate_qmc) on the thimble; kind 'S' or a graviton kind (then also the scalar on the same
    contour, key Es).  Returns dict(E, err, Es, errs, logN (scalar) / logNg (graviton), diagnostics)."""
    nuf = max(nu, NUF) if nuf is None else nuf
    c, wc, L, gm = sad or saddle(nuf, Delta, J)
    ev = Evaluator(nu, Delta, J, c, wc, L, kind, nuf=nuf)
    valsS, valsG = [], []
    diag = dict(nbad=0, wind=0.0, wind_g=0.0, steps=0)
    for r in range(R):
        n = sphere_pts(m, seed0 + r)
        if anti:
            n = np.vstack([n, -n])
        fl = flow(ev, n, verbose=verbose)
        a = radial_scalar(fl)
        valsS.append(a.mean()*AREA6/(2*pi)**3.5)
        diag['nbad'] += fl['nbad']; diag['steps'] = max(diag['steps'], len(fl['t']))
        diag['wind'] = max(diag['wind'], fl['wind']); diag['nsub'] = max(diag.get('nsub', 0), fl['nsub_max'])
        diag['imdrift'] = max(diag.get('imdrift', 0.0), fl['imdrift']); diag['lost'] = max(diag.get('lost', -1e9), fl.get('lost', -1e9))
        if kind != 'S':
            ag, Pgc, pref0, wg = radial_grav(ev, fl)
            valsG.append(ag.mean(axis=-1)*AREA6/(2*pi)**3.5)
            diag['wind_g'] = max(diag['wind_g'], wg)
    valsS = np.array(valsS)
    out = dict(Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, diag=diag, wc=wc, saddle_grad=gm)
    # normalisations: exp(Pint)(w_c) (2 pi)^{7/2} det L  (n1 convention)
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamS, GS, QS = PW.Cr.build(zl, nu, Delta, J, pref=True)
    out['logN'] = lamS[0] + GS[0] + lj[0] + 3.5*np.log(2*pi) + np.log(np.linalg.det(L))
    if kind != 'S':
        valsG = np.array(valsG)
        out.update(E=valsG.mean(axis=0), err=valsG.std(axis=0, ddof=1)/np.sqrt(R), reps=valsG)
        PW.set_batch(1)
        zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
        lamF, G, pref = PW.P9.build_spin(zl, nu, Delta, J, ev.mons, ev.sol, ev.solt, kind)
        out['logNg'] = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))           # n1_run's logN (no (2pi)^{7/2})
    return out


if __name__ == '__main__':
    nu, Dl, J = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    m = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    t0 = time.time()
    r = integrate(nu, Dl, J, 'S', m=m, R=4, verbose=True)
    print(f"nu={nu} Delta={Dl} J={J}: Es = {r['Es']:.8f} +- {r['errs']:.1e} (rel {abs(r['errs']/r['Es']):.1e}); logN {r['logN']:.6f}; "
          f"diag {r['diag']}; {time.time() - t0:.0f}s")
