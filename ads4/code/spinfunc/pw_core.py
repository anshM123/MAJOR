"""
pw_core (HWIN_THIMBLE.md): per-wave ("travelling-wave thimble") evaluation of the exact CLS-D crystal integral at large J.

Integrand (HLJ_ANALYTIC §1.2; checked point-wise, pw_tests.py): the scalar crystal integrand of l1_crystal.build is
    exp(lamF + G_0 + logJac) * P_J(t),   t = Z.W / sqrt(Z^2 W^2),
with G_0 the J = 0 slow exponent (the J-dependent 3-pt powers cancel (Z^2 W^2)^{J/2} identically), up to a J-only constant sign.
Split P_J = W_+ + W_- (pw_wave: exact Legendre-Q waves, only logarithmic singularities at theta = 0).
  mode 0  : P_J itself on the n3 contour of the J-growing scalar flow function (= the z2/fw production contour);
  mode +-1: the wave W_+- alone, on the n3_hybrid contour C(rho1, nu_f) built on the saddle of the WAVE flow function
            Phi_+- = lamF + G_0 + logJac + log W_+-(t)   (at nu_f),
            found by Newton (holomorphic-stencil derivatives) with a homotopy lambda log W_+-, lambda: 0 -> 1, from the mode-0 saddle.
sqrt(Z^2 W^2) and every log are continued (n1_pw replay) from W_REF next to the crystal (n3_flow.set_branch convention), so the wave
labels are the same along every path.  The graviton kinds (grav) use pw_grav (exact per-wave spin-J contraction).
Normalisation: as z2_points.n3_block (E: (2 pi)^{-7/2}-normalised contour average; logN = Pint(w_c) + log det L);
log Y = log(E) + logN - log beta_{Delta,J} (+ constant J-only sign/normalisation that cancels between state and reference).
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import n1_direct as D
import n1_pw as PW
import l1_tps as T
import l1_crystal as Cr
import n3_flow as F
import n3_hybrid as Hy
import pw_wave as WV


def zw_geom(zl, M=None):
    """Z.W and Z^2 W^2 of l1_crystal.build (point-wise backend)."""
    M = M or PW.Cr
    xt0, ux, vx, uy, vy, xp, yp = zl
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = M.dot(x, x); x4 = M.scal(1.0/(-xsq), x); x2 = M.neg(x4); x1 = M.neg(y); x3 = y
    e = M.cpt(1.0, 1.0, 0.0)

    def msq(p, q):
        r = M.sub(p, q); return -M.dot(r, r)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    Z = M.sub(M.scal(m25, M.sub(x1, e)), M.scal(m15, M.sub(x2, e)))
    W = M.sub(M.scal(m45, M.sub(e, x3)), M.scal(m35, M.sub(e, x4)))
    return M.dot(Z, W), M.dot(Z, Z)*M.dot(W, W)


class WEv(F.Evaluator):
    """n3 Evaluator with the mode's spin factor (scalar) or the per-wave graviton contraction (kind != 'S')."""

    def __init__(self, nu, Delta, J, c, wc, L, mode=0, kind='S', nuf=None, lam=1.0):
        super().__init__(nu, Delta, J, c, wc, L, 'S', nuf=nuf)
        self.mode, self.kind, self.lam = mode, kind, lam
        if kind != 'S':
            import l5_spin as S
            self.mons, self.sol, _ = S.conserved(Delta, J)
            _, self.solt, _ = S.conserved(3.0 - Delta, J)

    def _common(self, xi, nu):
        zl, lj = self._z(xi)
        lam, G0 = PW.Cr.build(zl, nu, self.Delta, 0)
        zw, zz = zw_geom(zl)
        t = zw*PW.powr(zz, -0.5)
        return zl, lam + G0 + lj, t

    def pint(self, xi):
        zl, P, t = self._common(xi, self.nu)
        return P, WV.wave_derivs(self.J, t, self.mode, 0)[0]

    def phi(self, xi):
        zl, P, t = self._common(xi, self.nuf)
        Q = WV.wave_derivs(self.J, t, self.mode, 0)[0]
        if self.mode == 0:
            return P, P, Q
        return P + self.lam*PW.log(Q), P, Q

    def grav(self, xi):
        import pw_grav as PG
        return PG.grav(self, xi)


# ------------------------------------------------------------------ saddles
def _gh(Fv, h):
    Fq = Fv[:, 1:].reshape(Fv.shape[0], len(F.DIRS), 4)
    d1 = (Fq*(1/F.OM)[None, None, :]).sum(-1)/(4*h)
    d2 = (Fq*(1/F.OM**2)[None, None, :]).sum(-1)/(2*h**2)
    g = d1[:, :7]; Hd = d2[:, :7]
    H = np.zeros((Fv.shape[0], 7, 7), complex)
    for i in range(7):
        H[:, i, i] = Hd[:, i]
    for k, (i, j) in enumerate(F.PAIRS):
        hij = 0.5*(d2[:, 7 + k] - Hd[:, i] - Hd[:, j]); H[:, i, j] = hij; H[:, j, i] = hij
    return g[0], H[0]


def _gradhess_w(nu, Delta, J, mode, c, w, nuf, lam, h=2e-3):
    ev = WEv(nu, Delta, J, c, w, np.eye(7, dtype=complex), mode, nuf=nuf, lam=lam)
    F.set_branch(ev, nstep=60, fns=('phi',))
    pts = F.STEN*(h/F.H_FD)
    (Phi, _, _), _, _ = ev.run(ev.phi, pts, None)
    return _gh(Phi[None, :], h)


_SAD = {}


def saddle(nuf, Delta, J, mode, verbose=False):
    """(c, wc, L) of the mode's flow function at nu_f (w-coordinates of the crystal c(nu_f))."""
    key = (round(float(nuf), 9), float(Delta), int(J), int(mode))
    if key in _SAD:
        return _SAD[key]
    c, w, L, gm = F.saddle(nuf, Delta, 0)          # = the J-growing scalar flow function of n3 (its J-terms cancel exactly)
    if mode != 0:
        w = w.copy()
        for lam in list(np.linspace(0.2, 1.0, 5)) + [1.0, 1.0, 1.0]:
            for it in range(12):
                g, H = _gradhess_w(nuf, Delta, J, mode, c, w, nuf, lam)
                dw = np.linalg.solve(H, g)
                s = 1.0 if np.abs(dw).max() <= 0.1 else 0.1/np.abs(dw).max()
                w = w - s*dw
                if np.abs(dw).max() < 1e-9:
                    break
            if verbose:
                print(f"   saddle mode {mode} lam={lam:.2f}: |g| {np.abs(g).max():.1e} |dw| {np.abs(dw).max():.1e}", flush=True)
        g, H = _gradhess_w(nuf, Delta, J, mode, c, w, nuf, 1.0)
        gm = float(np.abs(g).max())
        L = np.linalg.inv(T.csym_cholesky(-H)).T
    _SAD[key] = (c, w, L, gm)
    return _SAD[key]


NUF_LOW = float(os.environ.get('PW_NUF_LOW', '13'))


def nuf_of(nu):
    return max(nu, F.NUF) if nu >= NUF_LOW else NUF_LOW


def make_ev(nu, Delta, J, mode, kind='S', nuf=None, verbose=False):
    nuf = nuf_of(nu) if nuf is None else nuf
    c, wc, L, gm = saddle(nuf, Delta, J, mode, verbose)
    ev = WEv(nu, Delta, J, c, wc, L, mode, kind, nuf=nuf)
    wb = F.set_branch(ev, fns=('phi', 'pint') + (('grav',) if kind != 'S' else ()))
    ev.gm = gm; ev.wbranch = wb
    return ev


def block(evs, prop, m, R, seed0, rho1=3.5, gkinds=None):
    """E for the scalar evaluator(s) evs (list, same contour parameters: one per mode, each on its own contour) with the given proposal
    (mu, Ch, ldC); replicas R of 2^m antithetic Sobol points.  Returns per-ev dict(E, err, reps, logN, diag) (and graviton reps)."""
    mu, Ch, ldC = prop
    out = []
    for ev in evs:
        vals = []; gvals = {}
        diag = dict(nbad=0, n=0, nouter=0, wind=0.0, imdrift=0.0)
        gev = (gkinds or {}).get(id(ev), {})
        for r in range(R):
            eta = Hy._etas(m, seed0 + r); eta = np.concatenate([eta, -eta])
            acc = 0.0; accg = {k: None for k in gev}
            for c0 in range(0, eta.shape[0], Hy.CHUNK):
                et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
                wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
                ct = Hy.contour(ev, Z, rho1, keep_path=bool(gev))
                acc = acc + np.sum(wgt*ct['Fval'])
                d = ct['diag']; diag['nbad'] += int((~ct['ok']).sum()); diag['n'] += Z.shape[0]; diag['nouter'] += d['nouter']
                for kk in ('wind', 'imdrift'):
                    diag[kk] = max(diag[kk], d.get(kk, 0.0))
                for k, ge in gev.items():
                    Fg, Pgc, pref0, okg = grav_on_paths(ge, ct, Z.shape[0])
                    v = np.sum(wgt*Fg, axis=-1)
                    accg[k] = v if accg[k] is None else accg[k] + v
            vals.append(acc/eta.shape[0])
            for k in gev:
                gvals.setdefault(k, []).append(accg[k]/eta.shape[0])
        vals = np.array(vals)
        (Pc, _), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
        o = dict(E=vals.mean(), err=vals.std(ddof=1)/np.sqrt(R), reps=vals, logN=Pc[0] + np.log(np.linalg.det(ev.L)), diag=diag,
                 mode=ev.mode, nuf=ev.nuf, wc=ev.wc, gm=getattr(ev, 'gm', None))
        o['g'] = {}
        for k, ge in gev.items():
            vg = np.array(gvals[k])
            (Pg, _), _, _ = ge.run(ge.grav, np.zeros((1, 7), complex), None)
            o['g'][k] = dict(reps=vg, E=vg.mean(0), err=vg.std(0, ddof=1)/np.sqrt(R), logN=Pg[0] + np.log(np.linalg.det(ge.L)))
        out.append(o)
    return out


def logY(o, Delta, J):
    """log Y (z3 convention) of a list of per-mode block results summed with common normalisation; returns (value, replicas)."""
    ref = o[0]['logN']
    tot = sum(np.asarray(x['reps'])*np.exp(x['logN'] - ref) for x in o)
    return np.log(tot.mean()) + ref - Cr.log_beta(Delta, J), np.log(tot) + ref - Cr.log_beta(Delta, J)


def _gravlog(self, xi):
    import pw_grav as PG
    return PG.grav(self, xi, logs_only=True)


WEv.gravlog = _gravlog


def make_gev(ev, kind):
    """graviton evaluator on the SAME contour as the scalar evaluator ev (kind 'EE' / 'C3H', same mode)."""
    ge = WEv(ev.nu, ev.Delta, ev.J, ev.c, ev.wc, ev.L, ev.mode, kind, nuf=ev.nuf, lam=ev.lam)
    F.set_branch(ge, fns=('gravlog',))
    ge.saddle_refs['grav'] = ge.saddle_refs['gravlog']
    return ge


def grav_on_paths(ge, ct, Nn):
    """as n3_hybrid.grav_on_paths, but the expensive spin-J contraction only at the final points (the path replay uses the identical
    log sequence without the contraction, pw_grav logs_only)."""
    (Pg0, pref0), ref0, _ = ge.run(ge.grav, np.zeros((1, 7), complex), None)
    Pgc = Pg0[0]
    refs = [np.repeat(r, Nn, axis=-1) for r in ref0]
    wmax = np.zeros(Nn)
    for k, (idx, xs) in enumerate(ct['path']):
        (_, _), new, wind = ge.run(ge.gravlog, xs, [r[..., idx] for r in refs], rep=1)
        for a, b in zip(refs, new):
            a[..., idx] = b
        wmax[idx] = np.maximum(wmax[idx], wind)
    (Pg, pref), _, wind = ge.run(ge.grav, ct['xi'], refs, rep=1)
    wmax = np.maximum(wmax, wind)
    Fg = ct['detJ']*np.exp(Pg - Pgc)*pref
    ok = ct['ok'] & (wmax < Hy.WCUT) & np.all(np.isfinite(Fg.reshape(-1, Nn)), axis=0)
    return np.where(ok, Fg, 0.0), Pgc, pref0[..., 0], ok


Hy_grav_on_paths = grav_on_paths
