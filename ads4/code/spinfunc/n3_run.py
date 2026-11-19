"""
n3_run (doc 25): graviton C2 and C3H crystal integrals at SMALL nu on the wall-free hybrid contour of n3_hybrid (nu_f = max(nu, 15),
rho1 = 3.5, adaptive Gaussian proposal, RQMC), both kinds and the scalar on the SAME samples, at the l9 Delta = 100 row nodes
x in N3_XS (default 0.15, 0.1, 0.05; nu = x m(100, J)) and their heavy references (Delta_ref = 2e4, same nu).
Stored in n1_run's format (cache/n3_pts_<kind>_J<J>.pkl, key (Delta, round(nu, 9)); ljc = True; logN = log[exp(Pint_g)(w_c) det L];
E, err, reps (per replica), Es, errs, repsS (scalar on the same samples), diagnostics), so that n1_analysis.grid_direct can use them.
Usage: OMP_NUM_THREADS=1 nohup python3 n3_run.py J [J ...]   (resumable; env N3_XS, N3_M (log2 pts/replica, 11), N3_R (8), N3_DS)
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np
from math import pi
import n3_hybrid as Hy
import n3_flow as F

XS = [float(v) for v in os.environ.get('N3_XS', '0.15,0.1,0.05').split(',')]
MQ = int(os.environ.get('N3_M', '11'))
RQ = int(os.environ.get('N3_R', '8'))
DS = [float(v) for v in os.environ.get('N3_DS', '100,20000').split(',')]
KINDS = ('C2', 'C3H')


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def fn(kind, J):
    return os.path.join(HERE, 'cache', 'n3_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


ADAPT_G = os.environ.get('N3_ADAPT_G', '1') == '1'
NUF_LOW = float(os.environ.get('N3_NUF_LOW', '13'))     # contour nu_f for nu < NUF_LOW (above the soft-block fold nu* ~ 11-12)


def adapt_g(evS, evs, rho1, s=1.5, m=11, infl=1.3, seed=31, its=2):
    """Gaussian proposal moment-matched to the GRAVITON integrands (sum over kinds and elements of |F_ij|/<|F_ij|>)."""
    mu = np.zeros(7); C = s*s*np.eye(7)
    for it in range(its):
        eta = Hy._etas(m, seed + it)
        w_, U = np.linalg.eigh(C); Ch = (U*np.sqrt(w_)) @ U.T
        Z = mu[None, :] + eta @ Ch.T
        ct = Hy.contour(evS, Z, rho1, keep_path=True)
        a = np.zeros(Z.shape[0])
        for k, e in evs.items():
            Fg = Hy.grav_on_paths(e, ct, Z.shape[0])[0]
            A_ = np.abs(Fg).reshape(-1, Z.shape[0])
            a += np.sum(A_/np.maximum(A_.mean(axis=1, keepdims=True), 1e-300), axis=0)
        a = a*np.exp(0.5*np.sum(eta*eta, axis=1)); a = a/a.sum()
        mu = (a[:, None]*Z).sum(0); d = Z - mu
        C = infl*(a[:, None, None]*d[:, :, None]*d[:, None, :]).sum(0)
        w_, U = np.linalg.eigh(C); w_ = np.clip(w_, 0.3, 16.0); C = (U*w_) @ U.T
    w_, U = np.linalg.eigh(C)
    return mu, (U*np.sqrt(w_)) @ U.T, float(np.sum(np.log(w_)))


def point(J, Delta, nu, m=MQ, R=RQ, rho1=3.5, seed0=4242):
    nuf = max(nu, F.NUF) if nu >= NUF_LOW else NUF_LOW
    c, wc, L, _ = F.saddle(nuf, Delta, J)
    evS = F.Evaluator(nu, Delta, J, c, wc, L, 'S', nuf=nuf); F.set_branch(evS); Hy.attach_nuflow(evS)
    evs = {k: F.Evaluator(nu, Delta, J, c, wc, L, k, nuf=nuf) for k in KINDS}
    for e in evs.values():
        F.set_branch(e)
    mu, Ch, ldC = adapt_g(evS, evs, rho1) if ADAPT_G else Hy.adapt(evS, rho1)
    valsS = []; valsG = {k: [] for k in KINDS}
    diag = dict(nbad=0, nbad_g=0, wind=0.0, imdrift=0.0, sdrift=0.0, n=0, nouter=0)
    for r in range(R):
        eta = Hy._etas(m, seed0 + r); eta = np.concatenate([eta, -eta])
        accS = 0.0; accG = {k: None for k in KINDS}
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = Hy.contour(evS, Z, rho1, keep_path=True)
            accS = accS + np.sum(wgt*ct['Fval'])
            d = ct['diag']; diag['nbad'] += int((~ct['ok']).sum()); diag['n'] += Z.shape[0]; diag['nouter'] += d['nouter']
            for kk in ('wind', 'imdrift', 'sdrift'):
                diag[kk] = max(diag[kk], d[kk])
            for k in KINDS:
                Fg, Pgc, pref0, okg = Hy.grav_on_paths(evs[k], ct, Z.shape[0])
                v = np.sum(wgt*Fg, axis=-1)
                accG[k] = v if accG[k] is None else accG[k] + v
                diag['nbad_g'] += int((~okg).sum())
        valsS.append(accS/eta.shape[0])
        for k in KINDS:
            valsG[k].append(accG[k]/eta.shape[0])
    valsS = np.array(valsS)
    (Pc, _), _, _ = evS.run(evS.pint, np.zeros((1, 7), complex), None)
    logNs = Pc[0] + np.log(np.linalg.det(L))                    # scalar, same convention as the graviton logN (no (2pi)^{7/2})
    out = {}
    for k in KINDS:
        vg = np.array(valsG[k])
        (Pg, _), _, _ = evs[k].run(evs[k].grav, np.zeros((1, 7), complex), None)
        out[k] = dict(ljc=True, E=vg.mean(axis=0), err=vg.std(axis=0, ddof=1)/np.sqrt(R), reps=vg, logN=Pg[0] + np.log(np.linalg.det(L)),
                      Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, logNs=logNs, nu=nu, Delta=Delta, J=J,
                      m=mphys(Delta, J), nuf=nuf, rho1=rho1, mq=m, R=R, diag=dict(diag), method='n3_hybrid', tflow=Hy.TFLOW)
    return out


if __name__ == '__main__':
    for J in [int(a) for a in sys.argv[1:]]:
        for x in XS:
            for Dl in DS:
                nu = x*mphys(100.0, J)
                key = (Dl, round(nu, 9))
                if all(key in load(k, J) for k in KINDS):
                    continue
                t = time.time()
                res = point(J, Dl, nu)
                for k in KINDS:
                    db = load(k, J); db[key] = res[k]
                    tmp = fn(k, J) + '.tmp%d' % os.getpid()
                    pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, fn(k, J))
                r0 = res['C2']
                rel = {k: float(np.max(np.abs(res[k]['err']))/np.max(np.abs(res[k]['E']))) for k in KINDS}
                print(f"J={J} Delta={Dl:.0f} x={x} nu={nu:.4f}: max rel err {rel}, scalar rel err {abs(r0['errs']/r0['Es']):.1e}, "
                      f"diag {r0['diag']}, {time.time() - t:.0f}s", flush=True)
