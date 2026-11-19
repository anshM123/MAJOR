"""
n3_thimble (doc 25): the crystal integral on the Lefschetz thimble (n3_flow) with a SPARSE-GRID rule in the thimble parameters.

Parametrisation of the contour (n3_flow): zeta in R^7 -> xi(zeta) = Flow_{t}(r0 n), n = zeta/|zeta|, t = log(|zeta|/r0), the downward
gradient flow of the scalar flow function Phi_{nu_f} out of its saddle (xi = zeta for |zeta| <= r0).  Jacobian
   det(d xi/d zeta) = det[d xi/dt, V_1, ..., V_6] / |zeta|^7,     V_a(0) = r0 e_a,  dV/dt = -conj(Hess Phi V)  (RK4 with the flow).
Integral (n1 normalisation):  E = (2 pi)^{-7/2} int d^7zeta det(d xi/d zeta) exp(Pint - Pint_c) pref
                                = sum_k W_k F(zeta_k) e^{|zeta_k|^2/2}       (Smolyak Gauss-Hermite rule, n1_direct.smolyak),
optionally with a random rotation O of the rule (zeta_k -> s O zeta_k; the Gaussian weight is rotation invariant, so every rotated
rule is a valid rule; the spread over rotations is the error estimate), scale s (then W_k -> W_k s^7 e^{(s^2-1)|z|^2/2}).
Each node needs one flow line (integrated to exactly t_k with adaptive RK4 substeps); the path is stored so that the graviton integrand
(build_spin) can be continued along it (logs replayed point by point from the saddle).
Gradient / Hessian of Phi: central differences on a 57-point stencil (h = N3_H2), O(h^2).
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import pi
from scipy.stats import special_ortho_group
import n1_direct as D
import n1_pw as PW
import n3_flow as F

H2 = float(os.environ.get('N3_H2', '0.003'))
R0 = float(os.environ.get('N3_R0T', '0.05'))
HMAX = float(os.environ.get('N3_HMAX', '0.2'))       # max RK4 step in rho
GSTEP = float(os.environ.get('N3_GSTEP', '0.08'))     # max |d xi| per step
HSTEP = float(os.environ.get('N3_HSTEP', '0.4'))      # max h ||H|| per step
WCUT = float(os.environ.get('N3_WCUT', '1.0'))
NUF = float(os.environ.get('N3_NUF', '15'))
MAXSTEPS = int(os.environ.get('N3_MAXSTEPS', '400'))


def _stencil():
    P = [np.zeros(7)]
    for i in range(7):
        P += [H2*np.eye(7)[i], -H2*np.eye(7)[i]]
    pairs = [(i, j) for i in range(7) for j in range(i + 1, 7)]
    for (i, j) in pairs:
        d = H2*(np.eye(7)[i] + np.eye(7)[j])
        P += [d, -d]
    return np.array(P, complex), pairs


STEN, PAIRS = _stencil()
NS = STEN.shape[0]                                     # 57


def grad_hess(Fv):
    f0 = Fv[:, 0]
    fp = Fv[:, 1:15:2]; fm = Fv[:, 2:15:2]
    g = (fp - fm)/(2*H2)
    Hd = (fp + fm - 2*f0[:, None])/H2**2
    H = np.zeros((Fv.shape[0], 7, 7), complex)
    for i in range(7):
        H[:, i, i] = Hd[:, i]
    for k, (i, j) in enumerate(PAIRS):
        dd = (Fv[:, 15 + 2*k] + Fv[:, 16 + 2*k] - 2*f0)/H2**2           # (e_i + e_j)^T H (e_i + e_j)
        hij = 0.5*(dd - Hd[:, i] - Hd[:, j])
        H[:, i, j] = hij; H[:, j, i] = hij
    return g, H


def flow_to(ev, Z, keep_path=False):
    """ACTION-PARAMETRISED thimble: node zeta = rho n -> the point of the flow line out of r0 n at which Re(Phi_c - Phi) = rho^2/2
    exactly (flow integrated in rho: d xi/d rho = -conj(g) rho/|g|^2, dV/d rho = -conj(H V) rho/|g|^2), so that the Gaussian weight of
    the rule is exact and the rest of the integrand (Jacobian x prefactor) is smooth.  Jacobian:
         det(d xi/d zeta) = det[-conj(g), V_1..V_6] / (|g|^2 rho^5)      (V: t-flow tangent vectors; radial column along xi-dot)
    Returns dict: xi (final), Fval = det exp(Pint - Pint_c) Q, detJ, ok, path (if keep_path), diagnostics."""
    Nn = Z.shape[0]
    rho_t = np.linalg.norm(Z, axis=1); n = Z/rho_t[:, None]
    (Phi0, Pint0, Q0), ref0, _ = ev.run(ev.phi, np.zeros((1, 7), complex), None)
    Phic = Phi0[0]
    (Pi0, Qi0), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
    Pintc, Qc = Pi0[0], Qi0[0]
    refs = [np.repeat(r, Nn, axis=-1) for r in ref0]
    E = F.frames(n)
    xi = R0*n.astype(complex); V = R0*E[:, :, 1:].astype(complex)
    act = np.ones(Nn, bool); bad = np.zeros(Nn, bool); nst = np.zeros(Nn, int)
    path = [] if keep_path else None
    out = dict(wind=0.0, imdrift=0.0, maxsteps=0, sdrift=0.0)

    def field(xs, Vs, rr, idx, upd):
        pts = (xs[:, None, :] + STEN[None, :, :]).reshape(-1, 7)
        (Phi, _, _), new, wind = ev.run(ev.phi, pts, [r[..., idx] for r in refs], rep=NS)
        Phi = Phi.reshape(-1, NS); g, H = grad_hess(Phi)
        if upd:
            cen = np.arange(0, pts.shape[0], NS)
            for a, b in zip(refs, new):
                a[..., idx] = b[..., cen]
        g2 = np.sum(np.abs(g)**2, axis=1)
        fac = (rr/g2)
        return (-np.conj(g)*fac[:, None], -np.conj(np.einsum('nij,njk->nik', H, Vs))*fac[:, None, None], Phi[:, 0],
                wind.reshape(-1, NS).max(axis=1), g, H)

    # starting rho from the actual action at r0 n
    k0 = field(xi, V, np.ones(Nn), np.arange(Nn), True)
    rho = np.sqrt(np.maximum(2*(Phic - k0[2]).real, 1e-30))
    if keep_path:
        path.append((np.arange(Nn), xi.copy()))
    first = True
    while act.any():
        idx = np.where(act)[0]
        if first:
            k1x, k1v, Phi, w1, g, H = k0[0]*rho[:, None], k0[1]*rho[:, None, None], k0[2], k0[3], k0[4], k0[5]
            first = False
        else:
            k1x, k1v, Phi, w1, g, H = field(xi[idx], V[idx], rho[idx], idx, True)
        bad[idx] |= (w1 > WCUT) | ~np.isfinite(k1x).all(axis=1)
        out['imdrift'] = max(out['imdrift'], float(np.abs((Phi - Phic).imag).max()))
        out['sdrift'] = max(out['sdrift'], float(np.abs((Phic - Phi).real - 0.5*rho[idx]**2).max()))
        gn = np.linalg.norm(g, axis=1); hn = np.linalg.norm(H.reshape(H.shape[0], -1), axis=1)
        sp = np.linalg.norm(k1x, axis=1)                                    # |d xi/d rho|
        h = np.minimum.reduce([np.full(idx.size, HMAX), GSTEP/np.maximum(sp, 1e-12),
                               HSTEP/np.maximum(hn*rho[idx]/np.maximum(gn**2, 1e-30), 1e-12), rho_t[idx] - rho[idx]])
        h = np.where(np.isfinite(h) & (h > 0), h, 0.0)
        hh = h[:, None]; hv = h[:, None, None]
        x0, v0, r0_ = xi[idx], V[idx], rho[idx]
        r2 = field(x0 + 0.5*hh*k1x, v0 + 0.5*hv*k1v, r0_ + 0.5*h, idx, False)
        r3 = field(x0 + 0.5*hh*r2[0], v0 + 0.5*hv*r2[1], r0_ + 0.5*h, idx, False)
        r4 = field(x0 + hh*r3[0], v0 + hv*r3[1], r0_ + h, idx, False)
        bad[idx] |= (r2[3] > WCUT) | (r3[3] > WCUT) | (r4[3] > WCUT)
        xi[idx] = x0 + hh/6*(k1x + 2*r2[0] + 2*r3[0] + r4[0])
        V[idx] = v0 + hv/6*(k1v + 2*r2[1] + 2*r3[1] + r4[1])
        rho[idx] += h; nst[idx] += 1
        out['wind'] = max(out['wind'], float(max(w1.max(), r2[3].max(), r3[3].max(), r4[3].max())))
        if keep_path:
            path.append((idx.copy(), xi[idx].copy()))
        done = (rho[idx] >= rho_t[idx] - 1e-12) | bad[idx] | (nst[idx] >= MAXSTEPS)
        act[idx[done]] = False
    out['maxsteps'] = int(nst.max())
    idx = np.arange(Nn)
    _, _, Phi, w1, g, H = field(xi, V, rho, idx, True)
    (Pint, Q), _, w2 = ev.run(ev.pint, xi, [r for r in refs[:_nref(ev)]], rep=1)
    bad |= (w1 > WCUT) | (w2 > WCUT) | ((nst >= MAXSTEPS) & (rho < rho_t - 1e-12))
    g2 = np.sum(np.abs(g)**2, axis=1)
    detJ = np.linalg.det(np.concatenate([-np.conj(g)[:, :, None], V], axis=2))/(g2*rho_t**5)
    Fval = detJ*np.exp(Pint - Pintc)*Q
    ok = ~bad & np.isfinite(Fval)
    out['sdrift'] = max(out['sdrift'], float(np.abs((Phic - Phi).real - 0.5*rho_t**2)[ok].max()) if ok.any() else 0.0)
    out.update(xi=xi, Fval=np.where(ok, Fval, 0.0), detJ=detJ, Pint=Pint, Q=Q, Phi=Phi, Phic=Phic, Pintc=Pintc, Qc=Qc, ok=ok,
               lost_mag=float(np.max(np.where(~ok, -0.5*rho_t**2, -1e9))) if (~ok).any() else -1e9, nsteps=nst, path=path)
    return out


def _nref(ev):
    if not hasattr(ev, 'nref'):
        ev.nref = len(ev.run(ev.pint, np.zeros((1, 7), complex), None)[1])
    return ev.nref


def grav_on_paths(ev, fl, Nn):
    """graviton integrand at the nodes, logs continued along the stored paths (pref matrices; returns F (na, nc, Nn) with the same
    det factor and the graviton normalisation exp(-Pg_c))."""
    (Pg0, pref0), ref0, _ = ev.run(ev.grav, np.zeros((1, 7), complex), None)
    Pgc = Pg0[0]
    refs = [np.repeat(r, Nn, axis=-1) for r in ref0]
    wmax = np.zeros(Nn)
    last = None
    for (idx, xs) in fl['path']:
        (Pg, pref), new, wind = ev.run(ev.grav, xs, [r[..., idx] for r in refs], rep=1)
        for a, b in zip(refs, new):
            a[..., idx] = b
        wmax[idx] = np.maximum(wmax[idx], wind)
    # final point = last path entry of every line (each line's last stored xi is its final xi)
    (Pg, pref), _, wind = ev.run(ev.grav, fl['xi'], refs, rep=1)
    wmax = np.maximum(wmax, wind)
    ok = fl['ok'] & (wmax < WCUT)
    Fg = fl['detJ']*np.exp(Pg - Pgc)*pref
    Fg = np.where(ok & np.all(np.isfinite(Fg.reshape(-1, Nn)), axis=0), Fg, 0.0)
    return Fg, Pgc, pref0[..., 0], ok


def rule(level, rot=None):
    X, W = D.smolyak(level)
    if rot is not None:
        X = X @ rot.T
    return X, W*np.exp(0.5*np.sum(X*X, axis=1))          # W_k e^{|eta_k|^2/2}


def _quadfeat(Z):
    """features of the quadratic model c + a.z - z.B.z/2 (B symmetric): [1, z_i, -z_i z_j/2 (i<=j, off-diag doubled)]."""
    cols = [np.ones(Z.shape[0])] + [Z[:, i] for i in range(7)]
    for i in range(7):
        for j in range(i, 7):
            cols.append(-0.5*Z[:, i]*Z[:, j]*(1.0 if i == j else 2.0))
    return np.array(cols).T


def adapt(ev, level=3, cap=4.0):
    """pilot on the Smolyak level-`level` nodes: fit log|h(zeta)|, h = F e^{|zeta|^2/2}, to c + a.z - z.B.z/2; returns mu, Chalf, logdetC
    of the adapted Gaussian N(mu, C), C = (1 + B)^{-1}, mu = C a (eigenvalues of C clipped to [1/cap, cap], |mu| <= 2.5)."""
    X, W = D.smolyak(level)
    r = np.linalg.norm(X, axis=1); nz = r > 1e-12
    fl = flow_to(ev, X[nz])
    h = fl['Fval']*np.exp(0.5*r[nz]**2)/fl['Qc']
    ok = fl['ok'] & (np.abs(h) > 0)
    A = _quadfeat(X[nz][ok]); y = np.log(np.abs(h[ok]))
    A = np.vstack([A, _quadfeat(np.zeros((1, 7)))]); y = np.r_[y, 0.0]
    sol = np.linalg.lstsq(A, y, rcond=None)[0]
    a = sol[1:8]; B = np.zeros((7, 7)); k = 8
    for i in range(7):
        for j in range(i, 7):
            B[i, j] = B[j, i] = sol[k]; k += 1
    resid = y - A @ sol
    M = np.eye(7) + B
    w, U = np.linalg.eigh(M)
    w = np.clip(w, 1/cap, cap)
    C = (U/w) @ U.T
    mu = C @ a
    if np.linalg.norm(mu) > 2.5:
        mu *= 2.5/np.linalg.norm(mu)
    Ch = (U/np.sqrt(w)) @ U.T
    return dict(mu=mu, Ch=Ch, logdetC=-np.sum(np.log(w)), a=a, B=B, fit_rms=float(np.sqrt(np.mean(resid**2))), npilot=int(ok.sum()))


def integrate(nu, Delta, J, kind='S', level=5, rots=0, nuf=None, adapt_level=3, verbose=False, seed=99, ad=None):
    """E (n1 normalisation) by a Smolyak rule on the thimble, in the adapted variable eta: zeta = mu + C^{1/2} eta.
    rots = 0: the unrotated rule; rots = R >= 2: R random rotations of the eta-rule (mean and spread)."""
    nuf = max(nu, NUF) if nuf is None else nuf
    c, wc, L, gm = F.saddle(nuf, Delta, J)
    ev = F.Evaluator(nu, Delta, J, c, wc, L, kind, nuf=nuf)
    if ad is None:
        ad = adapt(ev, adapt_level) if adapt_level else dict(mu=np.zeros(7), Ch=np.eye(7), logdetC=0.0)
    rng = np.random.default_rng(seed)
    Os = [None] if rots == 0 else [special_ortho_group.rvs(7, random_state=rng) for _ in range(rots)]
    valsS, valsG = [], []; diag = dict(nbad=0, lost=-1e9, wind=0.0, imdrift=0.0, maxsteps=0, nodes=0, sdrift=0.0,
                                       fit_rms=ad.get('fit_rms'), mu=float(np.linalg.norm(ad['mu'])))
    for O in Os:
        X, Wx = rule(level, O)
        Z = ad['mu'][None, :] + X @ ad['Ch'].T
        wk = Wx*np.exp(0.5*ad['logdetC'])
        r = np.linalg.norm(Z, axis=1)
        cen = r < 1e-9
        fl = flow_to(ev, Z[~cen], keep_path=(kind != 'S'))
        Fs = np.zeros(Z.shape[0], complex); Fs[~cen] = fl['Fval']; Fs[cen] = fl['Qc']
        valsS.append(np.sum(wk*Fs))
        diag['nbad'] += int((~fl['ok']).sum()); diag['lost'] = max(diag['lost'], fl['lost_mag']); diag['wind'] = max(diag['wind'], fl['wind'])
        diag['imdrift'] = max(diag['imdrift'], fl['imdrift']); diag['maxsteps'] = max(diag['maxsteps'], fl['maxsteps'])
        diag['sdrift'] = max(diag['sdrift'], fl['sdrift']); diag['nodes'] += Z.shape[0]
        if kind != 'S':
            Fg, Pgc, pref0, okg = grav_on_paths(ev, fl, int((~cen).sum()))
            G = np.zeros(Fg.shape[:-1] + (Z.shape[0],), complex); G[..., ~cen] = Fg; G[..., cen] = pref0[..., None]
            valsG.append(np.sum(wk*G, axis=-1))
            diag['nbad_g'] = diag.get('nbad_g', 0) + int((~okg).sum())
    valsS = np.array(valsS)
    out = dict(Es=valsS.mean(), errs=(valsS.std(ddof=1)/np.sqrt(len(valsS)) if len(valsS) > 1 else np.nan), repsS=valsS, diag=diag,
               wc=wc, nuf=nuf, ad=ad)
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamS, GS, QS = PW.Cr.build(zl, nu, Delta, J, pref=True)
    out['logN'] = lamS[0] + GS[0] + lj[0] + 3.5*np.log(2*pi) + np.log(np.linalg.det(L))
    if kind != 'S':
        valsG = np.array(valsG)
        out.update(E=valsG.mean(axis=0), err=(valsG.std(axis=0, ddof=1)/np.sqrt(len(valsG)) if len(valsG) > 1 else np.nan*valsG[0]),
                   reps=valsG)
        PW.set_batch(1)
        zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
        lamF, G, pref = PW.P9.build_spin(zl, nu, Delta, J, ev.mons, ev.sol, ev.solt, kind)
        out['logNg'] = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))
    return out


if __name__ == '__main__':
    nu, Dl, J = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    for lev in [int(a) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 else (3, 4, 5):
        t0 = time.time()
        r = integrate(nu, Dl, J, 'S', level=lev)
        print(f"nu={nu} Delta={Dl} J={J} level {lev}: Es = {r['Es']:.8e}; diag {r['diag']}; {time.time() - t0:.0f}s", flush=True)
