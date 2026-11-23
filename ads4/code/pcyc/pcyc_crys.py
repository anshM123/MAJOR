"""
pcyc_crys (PCYC_CLOSE.md section 3.2): the crystal-cycle integrals of the scalar and of the EE (k = 4, J = 0) integrand at one (Delta, nu),
with the EE integrand taken as Q x scalar POINT-WISE (RD_CLOSE section 1: Phi_EE = a Q Phi_scal exactly, Q single valued), so that no EE
branch bookkeeping enters:

    <Q>_crys = int_crys Q Phi_scal / int_crys Phi_scal,      compared with   <Q>_D' = I_E/I_S from pcyc_dquad.

P_cyc (crystal = c(nu) D' as twisted cycles, c the same for every form on the local system) implies <Q>_crys = <Q>_D' EXACTLY, at every
Delta and nu, with no normalisation, branch or Delta-only constant (they cancel in each ratio).

Crystal cycle = the project's: below the fold the n4 Lefschetz thimbles A+ + A- (n+ = n- = 1; saddles tracked from nu_c = 10,
contour n3_hybrid C(rho1 = 0.8, nu_f = nu)); above it the n3 hybrid contour C(3.5, nu_f) of the symmetric saddle S.  The contour points
and the scalar integrand on them come from the project's code (n4_thimble / n3_hybrid / n3_flow, imported read-only, no cache written);
the raw point z and the sheet of y0 at each contour point are read off by replaying the scalar integrand's own branch references, and Q
is evaluated there with rd_core.Qratio (RD_CLOSE, read-only).  Check: |exp(Pint)| recomputed from (z, y0) with rd_core agrees.
Proposal: Gaussian in the contour parameter, moment-matched to |F| (1 + |Q|/<|Q|>) (scalar and EE weights), re-used for all replicas.
Usage: python3 pcyc_crys.py Delta nu [m R]  -> cache/cr_D<Delta>_nu<nu>.pkl
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
SPF = os.path.join(HERE, '..', 'spinfunc'); RDD = os.path.join(HERE, '..', 'rd')
sys.path.insert(0, RDD); sys.path.insert(0, SPF)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from math import pi
import n4_thimble as N4
import n3_hybrid as Hy
import n3_flow as F
import n3_thimble as TH
import n1_pw as PW
import rd_core as C
np.seterr(all='ignore')

NUC = float(os.environ.get('PC_NUC', '10.0'))
NUSPLIT = float(os.environ.get('PC_NUSPLIT', '10.5'))
SEED = int(os.environ.get('PC_SEED', '4242'))
STORE = os.environ.get('PC_STORE', '1') == '1'          # keep per-point data (control variates, pcyc_cv.py)


def zy0(ev, xi, refs):
    """raw z (7, N), y0 (N), Pint and logJac(w) at contour points xi, the sheet of y0 replayed with the scalar integrand's references."""
    def fn(x):
        zl, lj = ev._z(x)
        out = PW.Cr.build(zl, ev.nu, ev.Delta, ev.J, parts=True)
        return np.asarray(out[2]['y0']), [np.asarray(v) for v in zl], out[0] + out[1] + lj, lj
    nref = TH._nref(ev)
    (y0, zl, Pint, lj), _, wind = ev.run(fn, xi, refs[:nref], rep=1)
    return np.array(zl), y0, Pint, lj, wind


def evaluate(ev, Z, rho1, c):
    ct = Hy.contour(ev, Z, rho1, keep_path=False)
    z, y0, Pint, lj, wind = zy0(ev, ct['xi'], ct['refs'])
    g = C.geometry([z[i] for i in range(7)], y0val=y0)
    Q = C.Qratio(g, ev.Delta)
    ok = ct['ok'] & np.isfinite(Q)
    # consistency of (z, y0-sheet) with the replayed integrand:  Re Pint - Re lj - sum Re(e_k) log|f_k| = -(nu/2) arg_branch(A),
    # so  k = [that/(-nu/2) - Arg A]/(2 pi)  must be an integer (a wrong y0 sheet changes |f_k| and A).
    Fl = C.Flow(ev.nu, ev.Delta, c)
    vals = {k: v for k, v in g.items() if not isinstance(v, list)}
    lA = Fl.logA(vals)
    Lrd = np.real(Fl.logabs_parts(vals))
    kk = ((np.real(Pint) - np.real(lj) - Lrd)/(-ev.nu/2) - np.angle(np.exp(lA)))/(2*pi)
    dk = np.abs(kk - np.round(kk))
    Fv = np.where(ok, ct['Fval'], 0.0)
    W = ev.wc[None, :] + ct['xi'] @ ev.L.T
    return dict(F=Fv, Q=np.where(ok, Q, 0.0), ok=ok, dk=np.where(ok, dk, 0.0), wind=wind, Qcond=C.Q_condition(z), diag=ct['diag'],
                W=W, y0=y0)


def adapt(ev, rho1, c, s=1.0, m=10, infl=1.3, seed=31, its=2):
    mu = np.zeros(7); Cv = s*s*np.eye(7)
    for it in range(its):
        eta = Hy._etas(m, seed + it)
        w_, U = np.linalg.eigh(Cv); Ch = (U*np.sqrt(w_)) @ U.T
        Zp = mu[None, :] + eta @ Ch.T
        o = evaluate(ev, Zp, rho1, c)
        A = np.abs(o['F']); B = np.abs(o['F']*o['Q'])
        a = A/max(A.mean(), 1e-300) + B/max(B.mean(), 1e-300)
        a = a*np.exp(0.5*np.sum(eta*eta, axis=1)); a = a/a.sum()
        mu = (a[:, None]*Zp).sum(0); d = Zp - mu
        Cv = infl*(a[:, None, None]*d[:, :, None]*d[:, None, :]).sum(0)
        w_, U = np.linalg.eigh(Cv); w_ = np.clip(w_, 0.3, 9.0 if rho1 < 2 else 16.0); Cv = (U*w_) @ U.T
    w_, U = np.linalg.eigh(Cv)
    return mu, (U*np.sqrt(w_)) @ U.T, float(np.sum(np.log(w_)))


def integrate_thimble(ev, rho1, c, m, R, prop):
    mu, Ch, ldC = prop
    pts = []
    S0 = []; S1 = []; diag = dict(nbad=0, n=0, maxdk=0.0, wind=0.0, maxQcond=0.0)
    for r in range(R):
        eta = Hy._etas(m, SEED + r); eta = np.concatenate([eta, -eta])
        a0 = 0.0; a1 = 0.0
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Zp = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            o = evaluate(ev, Zp, rho1, c)
            a0 = a0 + np.sum(wgt*o['F']); a1 = a1 + np.sum(wgt*o['F']*o['Q'])
            if STORE:
                pts.append(dict(r=r, wgt=wgt, F=o['F'], Q=o['Q'], W=o['W'], y0=o['y0'], ok=o['ok']))
            diag['nbad'] += int((~o['ok']).sum()); diag['n'] += Zp.shape[0]
            diag['wind'] = max(diag['wind'], float(np.max(o['wind'])))
            wq = np.abs(o['F'])*wgt; wq = wq/max(wq.sum(), 1e-300)
            diag['maxQcond'] = max(diag['maxQcond'], float(np.sum(wq*o['Qcond'])))
            wf = np.abs(o['F'])*wgt
            diag['maxdk'] = max(diag['maxdk'], float(np.max(o['dk'][wf > 1e-8*wf.max()])) if wf.max() > 0 else 0.0)
        S0.append(a0/eta.shape[0]); S1.append(a1/eta.shape[0])
    return np.array(S0), np.array(S1), diag, pts


def run(Delta, nu, m=10, R=8, verbose=True):
    t0 = time.time()
    out = dict(Delta=Delta, nu=nu, m=m, R=R, th={})
    if nu < NUSPLIT:
        nus = [NUC] + ([nu] if nu < NUC - 1e-9 else [])
        recs = N4.track(Delta, 0, nus, kinds=())
        for wh in ('A+', 'A-'):
            rec = recs[wh][-1]
            evS, _ = N4.evaluators(rec, Delta, 0, ())
            rho1 = 0.8
            prop = adapt(evS, rho1, rec['c'])
            S0, S1, diag, pts = integrate_thimble(evS, rho1, rec['c'], m, R, prop)
            (Pc, _), _, _ = evS.run(evS.pint, np.zeros((1, 7), complex), None)
            logN = Pc[0] + 3.5*np.log(2*pi) + np.log(np.linalg.det(rec['L']))
            out['th'][wh] = dict(S0=S0, S1=S1, logN=logN, diag=diag, wc=rec['wc'], prop=prop, c=rec['c'], L=rec['L'], pts=pts)
            if verbose:
                print(f'  {wh}: <Q> = {np.mean(S1)/np.mean(S0):.6e}  rel err S0 {np.std(S0)/np.sqrt(R)/abs(np.mean(S0)):.1e}  diag {diag}  '
                      f'({time.time() - t0:.0f}s)', flush=True)
    else:
        nuf = max(nu, 15.0) if nu >= 13 else 13.0
        c, wc, L, _ = F.saddle(nuf, Delta, 0)
        ev = F.Evaluator(nu, Delta, 0, c, wc, L, 'S', nuf=nuf)
        Hy.attach_nuflow(ev)
        rho1 = 3.5
        prop = adapt(ev, rho1, c, s=1.0)
        S0, S1, diag, pts = integrate_thimble(ev, rho1, c, m, R, prop)
        (Pc, _), _, _ = ev.run(ev.pint, np.zeros((1, 7), complex), None)
        logN = Pc[0] + 3.5*np.log(2*pi) + np.log(np.linalg.det(L))
        out['th']['S'] = dict(S0=S0, S1=S1, logN=logN, diag=diag, wc=wc, prop=prop, nuf=nuf, c=c, L=L, pts=pts)
    # combine (replica-wise)
    keys = list(out['th'].keys())
    l0 = out['th'][keys[0]]['logN']
    S0 = sum(out['th'][k]['S0']*np.exp(out['th'][k]['logN'] - l0) for k in keys)
    S1 = sum(out['th'][k]['S1']*np.exp(out['th'][k]['logN'] - l0) for k in keys)
    out.update(S0=S0, S1=S1, logN0=l0, Qcrys=np.mean(S1)/np.mean(S0), time=time.time() - t0)
    # jackknife error of the ratio
    jk = np.array([(S1.sum() - S1[i])/(S0.sum() - S0[i]) for i in range(R)])
    out['Qerr'] = float(np.sqrt((R - 1)/R*np.sum(np.abs(jk - jk.mean())**2)))
    fn = os.path.join(HERE, 'cache', 'cr_D%g_nu%g%s.pkl' % (Delta, nu, os.environ.get('PC_TAG', '')))
    pickle.dump(out, open(fn, 'wb'))
    if verbose:
        print(f'crystal D={Delta:g} nu={nu:g}: <Q>_crys = {out["Qcrys"]:.8e} +- {out["Qerr"]:.1e} (rel {out["Qerr"]/abs(out["Qcrys"]):.1e}); '
              f'{out["time"]:.0f}s -> {fn}', flush=True)
    return out


if __name__ == '__main__':
    Dl, nu = float(sys.argv[1]), float(sys.argv[2])
    m = int(sys.argv[3]) if len(sys.argv) > 3 else 10
    R = int(sys.argv[4]) if len(sys.argv) > 4 else 8
    run(Dl, nu, m, R)
