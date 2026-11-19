"""
n4_thimble (NEARGAP100_SMALLNU.md): crystal integrals at SMALL nu (4.5-11) BELOW the fold nu* of the crystal saddle, as a sum of the
Lefschetz thimbles of the sigma-conjugate saddle pair A+, A- of the integrand's OWN exponent Phi_nu (n3_flow / n1_direct.phi_w):

    Y(nu) = n+ Y[A+](nu) + n- Y[A-](nu)        (intersection numbers n+-, constant between Stokes points)

* Saddles: continued in nu (descending) from the symmetric crystal saddle S above the fold (n3_flow.saddle), symmetry broken just
  below nu* by a small asymmetric kick (A+: the branch with Im w_1 > Im w_3 ... labelled by sign of Im(w1 - conj w3)).
* Sheets (log branches) at A+-(nu): TRACKED in nu.  At the calibration nu_c the logs at A+- are fixed as in n3_flow.set_branch
  (W_REF next to the crystal of the real nu, straight w-path); from there every log of the scalar flow function, the scalar
  integrand and every graviton integrand is continued along the z-path traced by the saddle as nu decreases (step <= DNU, straight
  in w between consecutive saddles, RSUB sub-steps).  In this tracked convention n+- can only change at a Stokes point.
* Thimble: the n3_hybrid contour C(rho1, nu_f = nu) built on A+- (linear disk |zeta| <= rho1, then downward flow of Phi_nu itself
  in the action parametrisation), i.e. the exact thimble up to the disk; integrand on it: non-oscillating Gaussian-like.
  rho1 small (default 0.8): the A+- thimbles meet the branch locus A = 0 at |zeta| ~ 1.5 along straight rays (H0 sect. 25.2),
  flow lines cannot (Im Phi const).
* n+-: determined at nu_c (default 10.0, validated range of n3_hybrid with nu_f = 13) by comparison with the n3_hybrid value of the
  same object, and checked at several nu, Delta, J (they must be the same small integers everywhere).
Functions: track(Delta, J, nus) -> saddle records (cached in cache/n4_saddles.pkl); thimble(rec, kinds, ...) -> integrals.
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import numpy as np
from math import pi
import n1_direct as D
import n1_pw as PW
import l1_tps as T
import l1_crystal as Cr
import n3_flow as F
import n3_hybrid as Hy
import n3_saddles as SD

SPK = os.path.join(HERE, 'cache', 'n4_saddles.pkl')
DNU = float(os.environ.get('N4_DNU', '0.25'))
RSUB = int(os.environ.get('N4_RSUB', '40'))


def newton(w, nu, Delta, J, it=80, cap=0.1, tol=1e-12):
    c, _ = Cr.saddle(nu, Delta)
    for k in range(it):
        wv = [T.var(i, w[i]) for i in range(T.NV)]
        g, H = T.grad_hess(D.phi_w(wv, c, nu, Delta, J))
        dw = np.linalg.solve(H, g)
        s = 1.0 if np.abs(dw).max() <= cap else cap/np.abs(dw).max()
        w = w - s*dw
        if not np.all(np.isfinite(w)):
            raise RuntimeError('newton diverged')
        if np.abs(dw).max() < tol:
            break
    wv = [T.var(i, w[i]) for i in range(T.NV)]
    g, H = T.grad_hess(D.phi_w(wv, c, nu, Delta, J))
    return c, w, H, float(np.abs(g).max())


def wshift(w, c_old, c_new):
    """the same z-point in the w-coordinates of crystal c_new."""
    w = w.copy()
    for i in D.LOGV:
        w[i] = w[i] + np.log(c_old[i]/c_new[i])
    return w


def asym(w):
    return abs(w[1] - np.conj(w[3])) + abs(w[2] - np.conj(w[4]))


class GEval(F.Evaluator):
    """Evaluator with an optional odd integrand (kind 'odd:C2' / 'odd:C3H': p8_odd_pw.P8.build, n0 = +1)."""

    def __init__(self, nu, Delta, J, c, wc, L, kind='S', nuf=None):
        if kind.startswith('odd:'):
            super().__init__(nu, Delta, J, c, wc, L, 'S', nuf=nuf)
            self.kind = kind; self.okind = kind[4:]
            import p8_odd_pw as P8M
            self._P8 = P8M.P8
        else:
            super().__init__(nu, Delta, J, c, wc, L, kind, nuf=nuf)

    def grav(self, xi):
        if not self.kind.startswith('odd:'):
            return super().grav(xi)
        zl, lj = self._z(xi)
        lamF, G, pref = self._P8.build(zl, self.nu, self.Delta, self.J, self.okind, 1.0)[:3]
        return lamF + G + lj, pref


def _fns(kinds):
    return [('phi', 'S')] + [('grav', k) for k in kinds]


def _continue(ev, fname, w_from, refs0, nstep):
    """logs of ev.<fname> continued along the straight w-path w_from -> ev.wc (refs0: reference logs at w_from, or None = principal)."""
    fn = getattr(ev, fname)
    Li = np.linalg.inv(ev.L)
    xs = [Li @ ((1 - s)*(w_from - ev.wc)) for s in np.linspace(0, 1, nstep + 1)]
    if refs0 is None:
        _, refs, w = ev.run(fn, xs[0][None, :].astype(complex), None)
    else:
        _, refs, w = ev.run(fn, xs[0][None, :].astype(complex), refs0, rep=1)
    wmax = float(np.max(w))
    for x in xs[1:]:
        _, refs, w = ev.run(fn, x[None, :].astype(complex), refs, rep=1)
        wmax = max(wmax, float(np.max(w)))
    return refs, wmax


def track(Delta, J, nus, kinds=(), nu_c=None, kick=1e-3, verbose=False):
    """records for both branches at the descending nus (nus[0] = calibration nu_c, below the fold).  Each record: dict(nu, c, wc, L, H, refs
    {(fname, kind): refs}, Phi0, wind).  The branch at nus[0] is set from W_REF (n3_flow.set_branch convention)."""
    nus = list(nus)
    # 1. symmetric saddle S down to just above the fold, then kick
    c, w, H, gm = None, np.zeros(7, complex), None, None
    nu_top = 20.0
    path = list(np.arange(nu_top, nus[0], -1.0)) + [nus[0]]
    wS = np.zeros(7, complex); cprev = None
    for v in path:
        cn, _ = Cr.saddle(v, Delta)
        if cprev is not None:
            wS = wshift(wS, cprev, cn)
        c, wS, H, gm = newton(wS, v, Delta, J)
        wS = 0.5*(wS + SD.sigma(wS)); cprev = cn
    out = {'A+': [], 'A-': []}
    # at nus[0]: from S (sym) kick both ways
    for which, sgn in (('A+', 1), ('A-', -1)):
        w0 = wS.copy(); w0[1] += sgn*kick*(1 + 1j); w0[3] -= sgn*kick*(1 - 1j)
        c, w, H, gm = newton(w0, nus[0], Delta, J)
        if asym(w) < 1e-4:
            raise RuntimeError(f'nu={nus[0]}: Newton returned to the symmetric saddle (above the fold?)')
        out[which].append((nus[0], c, w, H, gm))
    # label: A+ := Im(w1 - conj(w3)) > 0
    d0 = out['A+'][0][2][1] - np.conj(out['A+'][0][2][3])
    if d0.imag < 0:
        out['A+'], out['A-'] = out['A-'], out['A+']
    recs = {}
    for which in ('A+', 'A-'):
        nu0, c, w, H, gm = out[which][0]
        L = np.linalg.inv(T.csym_cholesky(-H)).T
        refs = {}; wind = 0.0
        for fname, kd in _fns(kinds):
            ev = GEval(nu0, Delta, J, c, w, L, kd, nuf=nu0)
            cr = c
            wf = F.W_REF.copy()                         # W_REF next to the crystal of nu0 (= c here: w = W_REF)
            r_, wm = _continue(ev, fname, wf, None, 200)
            refs[(fname, kd)] = r_; wind = max(wind, wm)
        rec = dict(nu=nu0, c=c, wc=w, L=L, H=H, grad=gm, refs=refs, wind=wind, Phi0=_phi0(c, w, nu0, Delta, J, L, refs))
        lst = [rec]
        for v in nus[1:]:
            prev = lst[-1]
            nsub = max(1, int(np.ceil((prev['nu'] - v)/DNU)))
            for k in range(1, nsub + 1):
                vv = prev['nu'] + (v - prev['nu'])*k/nsub if k < nsub else v
                cn, _ = Cr.saddle(vv, Delta)
                wf = wshift(prev['wc'], prev['c'], cn)
                cn, wn, Hn, gmn = newton(wf.copy(), vv, Delta, J)
                Ln = np.linalg.inv(T.csym_cholesky(-Hn)).T
                refs = {}; wind = 0.0
                for fname, kd in _fns(kinds):
                    ev = GEval(vv, Delta, J, cn, wn, Ln, kd, nuf=vv)
                    r_, wm = _continue(ev, fname, wf, prev['refs'][(fname, kd)], RSUB)
                    refs[(fname, kd)] = r_; wind = max(wind, wm)
                prev = dict(nu=vv, c=cn, wc=wn, L=Ln, H=Hn, grad=gmn, refs=refs, wind=wind, Phi0=_phi0(cn, wn, vv, Delta, J, Ln, refs))
                if verbose:
                    print(f'   {which} nu={vv:.3f}: |grad| {gmn:.1e}, asym {asym(wn):.3f}, step {np.abs(wn - wf).max():.3f}, wind {wind:.2f}, '
                          f'Phi0 {prev["Phi0"]:.4f}', flush=True)
            lst.append(prev)
        recs[which] = lst
    return recs


def _phi0(c, w, nu, Delta, J, L, refs):
    ev = GEval(nu, Delta, J, c, w, L, 'S', nuf=nu)
    ev.saddle_refs = {'phi': refs[('phi', 'S')]}
    (P, _, _), _, _ = ev.run(ev.phi, np.zeros((1, 7), complex), None)
    return complex(P[0])


def evaluators(rec, Delta, J, kinds):
    """scalar evaluator (with phi and pint branches) and graviton evaluators at the saddle record rec, with the tracked sheets."""
    nu = rec['nu']
    evS = GEval(nu, Delta, J, rec['c'], rec['wc'], rec['L'], 'S', nuf=nu)
    rp = rec['refs'][('phi', 'S')]
    evS.saddle_refs = {'phi': rp}
    nref = len(evS.run(evS.pint, np.zeros((1, 7), complex), None)[1])     # (plain at saddle; only the count is used)
    evS.saddle_refs['pint'] = rp[:nref]
    evS.nref = nref
    Hy.attach_nuflow(evS)
    evs = {}
    for k in kinds:
        e = GEval(nu, Delta, J, rec['c'], rec['wc'], rec['L'], k, nuf=nu)
        e.saddle_refs = {'phi': rp, 'pint': rp[:nref], 'grav': rec['refs'][('grav', k)]}
        evs[k] = e
    return evS, evs


def adapt(evS, evs, rho1, s=1.0, m=10, infl=1.3, seed=31, its=2):
    """Gaussian proposal moment-matched to the scalar (+ graviton) integrands on the thimble."""
    mu = np.zeros(7); C = s*s*np.eye(7)
    for it in range(its):
        eta = Hy._etas(m, seed + it)
        w_, U = np.linalg.eigh(C); Ch = (U*np.sqrt(w_)) @ U.T
        Z = mu[None, :] + eta @ Ch.T
        ct = Hy.contour(evS, Z, rho1, keep_path=bool(evs))
        A_ = np.abs(ct['Fval']); a = A_/max(A_.mean(), 1e-300)
        for k, e in evs.items():
            Fg = Hy.grav_on_paths(e, ct, Z.shape[0])[0]
            B_ = np.abs(Fg).reshape(-1, Z.shape[0])
            a = a + np.sum(B_/np.maximum(B_.mean(axis=1, keepdims=True), 1e-300), axis=0)
        a = a*np.exp(0.5*np.sum(eta*eta, axis=1)); a = a/a.sum()
        mu = (a[:, None]*Z).sum(0); d = Z - mu
        C = infl*(a[:, None, None]*d[:, :, None]*d[:, None, :]).sum(0)
        w_, U = np.linalg.eigh(C); w_ = np.clip(w_, 0.3, 9.0); C = (U*w_) @ U.T
    w_, U = np.linalg.eigh(C)
    return mu, (U*np.sqrt(w_)) @ U.T, float(np.sum(np.log(w_)))


def thimble(rec, Delta, J, kinds=(), rho1=0.8, m=10, R=8, seed0=4242, prop=None, verbose=False):
    """integrals on the thimble of the saddle record rec: scalar (Es, n1 normalisation, logN incl. (2 pi)^{7/2}) and graviton kinds
    (E, logN in n1_run convention, no (2 pi)^{7/2}), on common samples."""
    evS, evs = evaluators(rec, Delta, J, kinds)
    if prop is None:
        prop = adapt(evS, evs, rho1)
    mu, Ch, ldC = prop
    valsS = []; valsG = {k: [] for k in kinds}
    diag = dict(nbad=0, nbad_g=0, wind=0.0, imdrift=0.0, sdrift=0.0, n=0, nouter=0, maxsteps=0)
    t0 = time.time()
    for r in range(R):
        eta = Hy._etas(m, seed0 + r); eta = np.concatenate([eta, -eta])
        accS = 0.0; accG = {k: None for k in kinds}
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = Hy.contour(evS, Z, rho1, keep_path=bool(kinds))
            accS = accS + np.sum(wgt*ct['Fval'])
            d = ct['diag']; diag['nbad'] += int((~ct['ok']).sum()); diag['n'] += Z.shape[0]; diag['nouter'] += d['nouter']
            for kk in ('wind', 'imdrift', 'sdrift', 'maxsteps'):
                diag[kk] = max(diag[kk], d[kk])
            diag['wind'] = max(diag['wind'], d['wind_ray'])
            for k in kinds:
                Fg, Pgc, pref0, okg = Hy.grav_on_paths(evs[k], ct, Z.shape[0])
                v = np.sum(wgt*Fg, axis=-1)
                accG[k] = v if accG[k] is None else accG[k] + v
                diag['nbad_g'] += int((~okg).sum())
        valsS.append(accS/eta.shape[0])
        for k in kinds:
            valsG[k].append(accG[k]/eta.shape[0])
        if verbose:
            print(f'   replica {r}: Es {valsS[-1]:.5e} ({time.time() - t0:.0f}s)', flush=True)
    valsS = np.array(valsS)
    (Pc, _), _, _ = evS.run(evS.pint, np.zeros((1, 7), complex), None)
    ldL = np.log(np.linalg.det(rec['L']))
    out = dict(S=dict(Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, logN=Pc[0] + 3.5*np.log(2*pi) + ldL,
                      logNs=Pc[0] + ldL), diag=diag, prop=prop, rho1=rho1, mq=m, R=R)
    for k in kinds:
        vg = np.array(valsG[k])
        (Pg, _), _, _ = evs[k].run(evs[k].grav, np.zeros((1, 7), complex), None)
        out[k] = dict(E=vg.mean(axis=0), err=vg.std(axis=0, ddof=1)/np.sqrt(R), reps=vg, logN=Pg[0] + ldL)
    return out


def save_tracks(key, recs):
    db = pickle.load(open(SPK, 'rb')) if os.path.exists(SPK) else {}
    db[key] = recs
    tmp = SPK + '.tmp%d' % os.getpid()
    pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, SPK)


def load_tracks(key):
    if not os.path.exists(SPK):
        return None
    return pickle.load(open(SPK, 'rb')).get(key)


if __name__ == '__main__':
    Dl, J = float(sys.argv[1]), int(sys.argv[2])
    nus = [float(v) for v in sys.argv[3].split(',')]
    t0 = time.time()
    recs = track(Dl, J, nus, verbose=True)
    print(f'track {time.time() - t0:.0f}s')
    for wh in ('A+', 'A-'):
        for rec in recs[wh]:
            print(wh, rec['nu'], np.round(rec['wc'], 4), rec['Phi0'], rec['wind'])
