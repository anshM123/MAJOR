"""
n3_hybrid (doc 25): exact crystal integral on a WALL-FREE contour at small nu, by randomised QMC in the contour parameter zeta in R^7.

Contour C(rho1, nu_f) (any choice is a valid deformation of the integration cycle; the integral does not depend on it -- used as a check):
  * inner disk  |zeta| <= rho1 :  xi = zeta (the n1 linear contour through the saddle of Phi_{nu_f}), logs continued along the ray from
                                  the saddle in steps of <= RSTEP (n1 used two steps; that is what produced its "walls" at nu ~ 15);
  * outer part  |zeta| >  rho1 :  xi = point of the DOWNWARD gradient-flow line of Phi_{nu_f} started at rho1 n (n = zeta/|zeta|) at which
                                  Re Phi = Re Phi(rho1 n) - (|zeta|^2 - rho1^2)/2 (action parametrisation, flow integrated in rho).
  Along flow lines Im Phi is constant, so they cannot reach the branch locus of A^{i nu/2} (A = 0, where Im Phi -> infinity) and the
  integrand stays single valued; Re Phi decreases, so the contour leaves through the decay regions.
  Jacobian: 1 inside; det[-conj(g), V_1..V_6]/(|g|^2 |zeta|^5) outside (V_a(rho1) = rho1 e_a, dV/drho = -conj(H V) rho/|g|^2).
  Phi_{nu_f}: the scalar flow function of n3_flow (lamF + G + logJac + (J/2) log Z^2W^2) at nu_f; the integrand is evaluated at nu.
Integral (n1 normalisation): E = (2 pi)^{-7/2} int d^7zeta det exp(Pint_nu - Pint_nu(saddle)) pref = E_eta[ F(s eta) s^7 e^{|eta|^2/2} ],
scrambled Sobol (antithetic), error from R scramblings.  Saddle: n3_flow.saddle (nu_f; symmetric crystal branch) or a given (c, wc, L).
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import pi
from scipy.stats import qmc
from scipy.special import ndtri
import n1_direct as D
import n1_pw as PW
import n3_flow as F
import n3_thimble as TH

RSTEP = float(os.environ.get('N3_RSTEP', '0.25'))
WCUT = TH.WCUT
CHUNK = int(os.environ.get('N3_CHUNK', '2048'))
GSTRIDE = int(os.environ.get('N3_GSTRIDE', '1'))       # graviton log replay on every GSTRIDE-th stored path point (winding checked)


def _etas(m, seed):
    eng = qmc.Sobol(d=7, scramble=True, seed=seed)
    u = np.clip(eng.random_base2(m), 1e-15, 1 - 1e-15)
    return ndtri(u)


def contour(ev, Z, rho1, keep_path=False):
    """points of C(rho1, nu_f) for the parameters Z (Nn, 7).  Returns dict(xi, detJ, Pint, Q, ok, path, diagnostics)."""
    Nn = Z.shape[0]
    rho_t = np.linalg.norm(Z, axis=1)
    n = Z/np.maximum(rho_t, 1e-300)[:, None]
    outer = rho_t > rho1
    Zin = np.where(outer[:, None], rho1*n, Z)
    (Phi0, _, _), ref0, _ = ev.run(ev.phi, np.zeros((1, 7), complex), None)
    Phic = Phi0[0]
    (Pi0, Qi0), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
    Pintc, Qc = Pi0[0], Qi0[0]
    refs = [np.repeat(r, Nn, axis=-1) for r in ref0]
    K = max(1, int(np.ceil(rho1/RSTEP)))
    path = [] if keep_path else None
    wmax = np.zeros(Nn)
    for k in range(1, K + 1):
        pts = (k/K)*Zin.astype(complex)
        (_, _, _), new, wind = ev.run(ev.phi, pts, refs, rep=1)
        refs = new
        wmax = np.maximum(wmax, wind)
        if keep_path:
            path.append((np.arange(Nn), pts.copy()))
    xi = Zin.astype(complex).copy()
    detJ = np.ones(Nn, complex)
    Jm = np.repeat(np.eye(7, dtype=complex)[None], Nn, axis=0)
    bad = wmax > WCUT
    diag = dict(wind_ray=float(wmax.max()), nouter=int(outer.sum()), imdrift=0.0, sdrift=0.0, maxsteps=0)
    if outer.any():
        io = np.where(outer)[0]
        sub_refs = [r[..., io] for r in refs]
        fo = _flow_rho(ev, xi[io], n[io], rho1, rho_t[io], sub_refs, Phic, keep_path)
        xi[io] = fo['xi']; detJ[io] = fo['detJ']; Jm[io] = fo['Jm']; bad[io] |= fo['bad']
        for a, b in zip(refs, fo['refs']):
            a[..., io] = b
        if keep_path:
            for (ii, xs) in fo['path']:
                path.append((io[ii], xs))
        for kk in ('imdrift', 'sdrift', 'maxsteps', 'wind'):
            diag[kk] = fo[kk]
    if TFLOW > 0 and getattr(ev, 'evnu', None) is not None:
        gf = _gflow(ev.evnu, xi, Jm, refs, TFLOW, keep_path)
        xi = gf['xi']; detJ = gf['detJ']; bad |= gf['bad']; refs = gf['refs']
        if keep_path:
            path += gf['path']
        diag['gflow_wind'] = gf['wind']; diag['gflow_steps'] = gf['maxsteps']
    (Pint, Q), _, w2 = ev.run(ev.pint, xi, refs[:TH._nref(ev)], rep=1)
    bad |= w2 > WCUT
    Fval = detJ*np.exp(Pint - Pintc)*Q
    ok = ~bad & np.isfinite(Fval)
    return dict(xi=xi, detJ=detJ, Pint=Pint, Q=Q, Pintc=Pintc, Qc=Qc, Fval=np.where(ok, Fval, 0.0), ok=ok, path=path, diag=diag,
                refs=refs)


def _flow_rho(ev, xi, n, rho1, rho_t, refs, Phic, keep_path):
    Nn = xi.shape[0]
    E = F.frames(n)
    xi = xi.copy(); V = rho1*E[:, :, 1:].astype(complex)
    rho = np.full(Nn, rho1, float)
    act = np.ones(Nn, bool); bad = np.zeros(Nn, bool); nst = np.zeros(Nn, int)
    out = dict(imdrift=0.0, sdrift=0.0, wind=0.0, maxsteps=0)
    path = [] if keep_path else None
    STEN, NS = TH.STEN, TH.NS
    Phi_start = [None]

    def field(xs, Vs, rr, idx, upd):
        pts = (xs[:, None, :] + STEN[None, :, :]).reshape(-1, 7)
        (Phi, _, _), new, wind = ev.run(ev.phi, pts, [r[..., idx] for r in refs], rep=NS)
        Phi = Phi.reshape(-1, NS); g, H = TH.grad_hess(Phi)
        if upd:
            cen = np.arange(0, pts.shape[0], NS)
            for a, b in zip(refs, new):
                a[..., idx] = b[..., cen]
        g2 = np.sum(np.abs(g)**2, axis=1); fac = rr/g2
        return (-np.conj(g)*fac[:, None], -np.conj(np.einsum('nij,njk->nik', H, Vs))*fac[:, None, None], Phi[:, 0],
                wind.reshape(-1, NS).max(axis=1), g, H)

    while act.any():
        idx = np.where(act)[0]
        k1x, k1v, Phi, w1, g, H = field(xi[idx], V[idx], rho[idx], idx, True)
        if Phi_start[0] is None:
            Phi_start[0] = Phi.copy()                          # all lines active at the first step
        bad[idx] |= (w1 > WCUT) | ~np.isfinite(k1x).all(axis=1)
        out['imdrift'] = max(out['imdrift'], float(np.abs((Phi - Phi_start[0][idx]).imag).max()))
        out['sdrift'] = max(out['sdrift'], float(np.abs((Phi_start[0][idx] - Phi).real - 0.5*(rho[idx]**2 - rho1**2)).max()))
        gn = np.linalg.norm(g, axis=1); hn = np.linalg.norm(H.reshape(H.shape[0], -1), axis=1)
        sp = np.linalg.norm(k1x, axis=1)
        h = np.minimum.reduce([np.full(idx.size, TH.HMAX), TH.GSTEP/np.maximum(sp, 1e-12),
                               TH.HSTEP/np.maximum(hn*rho[idx]/np.maximum(gn**2, 1e-30), 1e-12), rho_t[idx] - rho[idx]])
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
        done = (rho[idx] >= rho_t[idx] - 1e-12) | bad[idx] | (nst[idx] >= TH.MAXSTEPS)
        act[idx[done]] = False
    idx = np.arange(Nn)
    _, _, Phi, w1, g, H = field(xi, V, rho, idx, True)
    bad |= (w1 > WCUT) | ((nst >= TH.MAXSTEPS) & (rho < rho_t - 1e-12))
    g2 = np.sum(np.abs(g)**2, axis=1)
    Mf = np.concatenate([(-np.conj(g)*(rho_t/g2)[:, None])[:, :, None], V/rho_t[:, None, None]], axis=2)   # columns in the frame basis
    Jm = np.einsum('nij,nkj->nik', Mf, E)                                                             # d xi/d zeta (Cartesian)
    detJ = np.linalg.det(Mf)
    out.update(xi=xi, detJ=detJ, Jm=Jm, bad=bad, refs=refs, path=path, maxsteps=int(nst.max()))
    return out


TFLOW = float(os.environ.get('N3_TFLOW', '0'))         # generalised-thimble time: flow C(rho1, nu_f) by the INTEGRAND's own Phi_nu


def _gflow(ev, xi, U, refs, T, keep_path, hmax=0.1):
    """downward gradient flow of Phi_nu (ev: Evaluator with nuf = nu, same coordinates) for time T of all points, with the tangent
    matrix U (d xi/d zeta) propagated (dU/dt = -conj(H U)); returns xi_T, det U_T, refs, path.  The flowed cycle is equivalent to the
    original one (flow lines avoid poles and cannot reach A = 0)."""
    Nn = xi.shape[0]
    xi = xi.copy(); U = U.copy(); t = np.zeros(Nn); bad = np.zeros(Nn, bool); act = np.ones(Nn, bool); nst = np.zeros(Nn, int)
    STEN, NS = TH.STEN, TH.NS
    path = [] if keep_path else None
    wmax = 0.0

    def field(xs, Us, idx, upd):
        pts = (xs[:, None, :] + STEN[None, :, :]).reshape(-1, 7)
        (Phi, _, _), new, wind = ev.run(ev.phi, pts, [r[..., idx] for r in refs], rep=NS)
        Phi = Phi.reshape(-1, NS); g, H = TH.grad_hess(Phi)
        if upd:
            cen = np.arange(0, pts.shape[0], NS)
            for a_, b_ in zip(refs, new):
                a_[..., idx] = b_[..., cen]
        return -np.conj(g), -np.conj(np.einsum('nij,njk->nik', H, Us)), wind.reshape(-1, NS).max(axis=1), g, H

    while act.any():
        idx = np.where(act)[0]
        k1x, k1u, w1, g, H = field(xi[idx], U[idx], idx, True)
        bad[idx] |= (w1 > WCUT) | ~np.isfinite(k1x).all(axis=1)
        gn = np.linalg.norm(g, axis=1); hn = np.linalg.norm(H.reshape(H.shape[0], -1), axis=1)
        h = np.minimum.reduce([np.full(idx.size, hmax), TH.GSTEP/np.maximum(gn, 1e-12), TH.HSTEP/np.maximum(hn, 1e-12), T - t[idx]])
        h = np.where(np.isfinite(h) & (h > 0), h, 0.0)
        hh = h[:, None]; hv = h[:, None, None]
        x0, u0 = xi[idx], U[idx]
        r2 = field(x0 + 0.5*hh*k1x, u0 + 0.5*hv*k1u, idx, False)
        r3 = field(x0 + 0.5*hh*r2[0], u0 + 0.5*hv*r2[1], idx, False)
        r4 = field(x0 + hh*r3[0], u0 + hv*r3[1], idx, False)
        bad[idx] |= (r2[2] > WCUT) | (r3[2] > WCUT) | (r4[2] > WCUT)
        wmax = max(wmax, float(max(w1.max(), r2[2].max(), r3[2].max(), r4[2].max())))
        xi[idx] = x0 + hh/6*(k1x + 2*r2[0] + 2*r3[0] + r4[0])
        U[idx] = u0 + hv/6*(k1u + 2*r2[1] + 2*r3[1] + r4[1])
        t[idx] += h; nst[idx] += 1
        if keep_path:
            path.append((idx.copy(), xi[idx].copy()))
        done = (t[idx] >= T - 1e-12) | bad[idx] | (nst[idx] >= TH.MAXSTEPS)
        act[idx[done]] = False
    idx = np.arange(Nn)
    _, _, w1, _, _ = field(xi, U, idx, True)
    bad |= (w1 > WCUT) | ((nst >= TH.MAXSTEPS) & (t < T - 1e-12))
    return dict(xi=xi, detJ=np.linalg.det(U), bad=bad, refs=refs, path=path, wind=wmax, maxsteps=int(nst.max()))


def grav_on_paths(ev, ct, Nn):
    """graviton integrand (build_spin, kind ev.kind) at the contour points, logs replayed along the stored paths."""
    (Pg0, pref0), ref0, _ = ev.run(ev.grav, np.zeros((1, 7), complex), None)
    Pgc = Pg0[0]
    refs = [np.repeat(r, Nn, axis=-1) for r in ref0]
    wmax = np.zeros(Nn)
    for k, (idx, xs) in enumerate(ct['path']):
        if GSTRIDE > 1 and (k % GSTRIDE) != GSTRIDE - 1 and k != len(ct['path']) - 1:
            continue
        (_, _), new, wind = ev.run(ev.grav, xs, [r[..., idx] for r in refs], rep=1)
        for a, b in zip(refs, new):
            a[..., idx] = b
        wmax[idx] = np.maximum(wmax[idx], wind)
    (Pg, pref), _, wind = ev.run(ev.grav, ct['xi'], refs, rep=1)
    wmax = np.maximum(wmax, wind)
    Fg = ct['detJ']*np.exp(Pg - Pgc)*pref
    ok = ct['ok'] & (wmax < WCUT) & np.all(np.isfinite(Fg.reshape(-1, Nn)), axis=0)
    return np.where(ok, Fg, 0.0), Pgc, pref0[..., 0], ok


def attach_nuflow(ev):
    """for TFLOW > 0 and nu != nu_f: the evaluator of Phi_nu in ev's coordinates (same log sequence, same branch references)."""
    if TFLOW > 0 and ev.nu != ev.nuf:
        e2 = F.Evaluator(ev.nu, ev.Delta, ev.J, ev.c, ev.wc, ev.L, 'S', nuf=ev.nu)
        if hasattr(ev, 'saddle_refs'):
            e2.saddle_refs = ev.saddle_refs
        ev.evnu = e2
    else:
        ev.evnu = None


def adapt(ev, rho1, s=1.5, m=11, infl=1.3, seed=31, its=2):
    """moment-matched Gaussian proposal N(mu, C) in zeta for |F| (pilot RQMC with N(mu, C) starting from N(0, s^2)); returns (mu, Ch, logdetC)."""
    mu = np.zeros(7); C = s*s*np.eye(7)
    for it in range(its):
        eta = _etas(m, seed + it)
        w_, U = np.linalg.eigh(C); Ch = (U*np.sqrt(w_)) @ U.T
        Z = mu[None, :] + eta @ Ch.T
        ct = contour(ev, Z, rho1)
        a = np.abs(ct['Fval'])*np.exp(0.5*np.sum(eta*eta, axis=1))
        a = a/a.sum()
        mu = (a[:, None]*Z).sum(0); d = Z - mu
        C = infl*(a[:, None, None]*d[:, :, None]*d[:, None, :]).sum(0)
        w_, U = np.linalg.eigh(C); w_ = np.clip(w_, 0.3, 16.0); C = (U*w_) @ U.T
    w_, U = np.linalg.eigh(C)
    return mu, (U*np.sqrt(w_)) @ U.T, float(np.sum(np.log(w_)))


def integrate(nu, Delta, J, kind='S', rho1=3.5, m=11, R=8, s=1.0, nuf=None, sad=None, seed0=4242, anti=True, verbose=False,
              branch_from=None, adaptive=False, prop=None):
    """E (n1 normalisation) on C(rho1, nu_f).  sad: optional (c, wc, L) (e.g. an asymmetric saddle).  Returns dict(Es, errs[, E, err],
    logN[, logNg], reps, diagnostics)."""
    nuf = max(nu, F.NUF) if nuf is None else nuf
    if sad is None:
        c, wc, L, _ = F.saddle(nuf, Delta, J)
    else:
        c, wc, L = sad[:3]
    ev = F.Evaluator(nu, Delta, J, c, wc, L, kind, nuf=nuf)
    if branch_from is not None:
        F.set_branch(ev, None if isinstance(branch_from, str) else branch_from)    # 'ref': W_REF next to the real-nu crystal
    attach_nuflow(ev)
    if prop is None and adaptive:
        prop = adapt(ev, rho1)
    mu, Ch, ldC = prop if prop is not None else (np.zeros(7), s*np.eye(7), 14*np.log(s))
    valsS, valsG = [], []
    diag = dict(nbad=0, wind=0.0, imdrift=0.0, sdrift=0.0, maxsteps=0, nouter=0, n=0, badw=0.0, totw=0.0, nbad_g=0)
    t0 = time.time()
    for r in range(R):
        eta = _etas(m, seed0 + r)
        if anti:
            eta = np.concatenate([eta, -eta])
        accS = 0.0; accG = None
        for c0 in range(0, eta.shape[0], CHUNK):
            et = eta[c0:c0 + CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = contour(ev, Z, rho1, keep_path=(kind != 'S'))
            accS = accS + np.sum(wgt*ct['Fval'])
            d = ct['diag']
            diag['nbad'] += int((~ct['ok']).sum()); diag['n'] += Z.shape[0]; diag['nouter'] += d['nouter']
            for kk in ('wind', 'imdrift', 'sdrift', 'maxsteps'):
                diag[kk] = max(diag[kk], d[kk])
            diag['wind'] = max(diag['wind'], d['wind_ray'])
            if kind != 'S':
                Fg, Pgc, pref0, okg = grav_on_paths(ev, ct, Z.shape[0])
                v = np.sum(wgt*Fg, axis=-1)
                accG = v if accG is None else accG + v
                diag['nbad_g'] += int((~okg).sum())
        valsS.append(accS/eta.shape[0])
        if kind != 'S':
            valsG.append(accG/eta.shape[0])
        if verbose:
            print(f"   replica {r}: Es {valsS[-1]:.6e}  ({time.time() - t0:.0f}s)", flush=True)
    valsS = np.array(valsS)
    out = dict(Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, diag=diag, wc=wc, nuf=nuf, rho1=rho1, prop=(mu, Ch, ldC))
    (Pc, _), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)          # (continued branch if set_branch was used)
    out['logN'] = Pc[0] + 3.5*np.log(2*pi) + np.log(np.linalg.det(L))
    if kind != 'S':
        valsG = np.array(valsG)
        out.update(E=valsG.mean(axis=0), err=valsG.std(axis=0, ddof=1)/np.sqrt(R), reps=valsG)
        (Pg, _), _, _ = ev.run(ev.grav, np.zeros((1, 7), complex), None)
        out['logNg'] = Pg[0] + np.log(np.linalg.det(L))
    return out


if __name__ == '__main__':
    nu, Dl, J = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    rho1 = float(sys.argv[4]) if len(sys.argv) > 4 else 3.5
    m = int(sys.argv[5]) if len(sys.argv) > 5 else 10
    nuf = float(sys.argv[6]) if len(sys.argv) > 6 else None
    s = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
    t0 = time.time()
    r = integrate(nu, Dl, J, 'S', rho1=rho1, m=m, R=8, nuf=nuf, s=s, verbose=False)
    print(f"nu={nu} Delta={Dl} J={J} rho1={rho1} nuf={r['nuf']} s={s}: Es = {r['Es']:.6e} +- {r['errs']:.1e} (rel {abs(r['errs']/r['Es']):.1e}); "
          f"logN {r['logN']:.6f}; diag {r['diag']}; {time.time() - t0:.0f}s", flush=True)
