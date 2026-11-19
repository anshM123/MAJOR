"""
n3_run_odd (doc 25.7): the PARITY-ODD (odd J, +-) integrands of p8_odd_pw (P8.build, kinds C2 / C3H, n0 = +1; imported read-only)
on the small-nu wall-free contour of n3_hybrid (same rule as n3_run: nu_f = max(nu, 15) for nu >= 13, nu_f = 13 below; rho1 = 3.5;
proposal adapted to the odd integrands), at the Delta = 100 row nodes x in N3_XS (nu = x m(100, J)) and the references Delta = 2e4.
Stored in p8's format (cache/n3_odd_pts_<kind>_J<J>.pkl, key (Delta, round(nu, 9)): E (1x1), err, reps, logN, Es, errs, repsS, ljc).
Phase/sign conventions of the shadow coefficient (PARITY_ODD_ROUTEA.md sect. 12) are NOT touched here: logN is the integral's own
normalisation, exactly as p8 stores it.
Usage: OMP_NUM_THREADS=1 nohup python3 n3_run_odd.py J [J ...]    (env N3_XS default 0.15,0.1; N3_M 12; N3_R 8)
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
warnings.filterwarnings('ignore')
import numpy as np
import n3_hybrid as Hy
import n3_flow as F
import n3_run as NR
import n1_pw as PW
import p8_odd_pw as P8M                       # read-only use of the odd integrand

XS = [float(v) for v in os.environ.get('N3_XS', '0.15,0.1').split(',')]
MQ = int(os.environ.get('N3_M', '12'))
RQ = int(os.environ.get('N3_R', '8'))
KINDS = ('C2', 'C3H')


class OddEval(F.Evaluator):
    def __init__(self, nu, Delta, J, c, wc, L, okind, nuf):
        super().__init__(nu, Delta, J, c, wc, L, 'S', nuf=nuf)
        self.kind = 'odd'; self.okind = okind

    def grav(self, xi):
        zl, lj = self._z(xi)
        lamF, G, pref = P8M.P8.build(zl, self.nu, self.Delta, self.J, self.okind, 1.0)[:3]
        return lamF + G + lj, pref


def fn(kind, J):
    return os.path.join(HERE, 'cache', 'n3_odd_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def point(J, Delta, nu, m=MQ, R=RQ, rho1=3.5, seed0=4242):
    nuf = max(nu, F.NUF) if nu >= NR.NUF_LOW else NR.NUF_LOW
    c, wc, L, _ = F.saddle(nuf, Delta, J)
    evS = F.Evaluator(nu, Delta, J, c, wc, L, 'S', nuf=nuf); F.set_branch(evS)
    evs = {k: OddEval(nu, Delta, J, c, wc, L, k, nuf) for k in KINDS}
    for e in evs.values():
        F.set_branch(e)
    mu, Ch, ldC = NR.adapt_g(evS, evs, rho1)
    valsS = []; valsG = {k: [] for k in KINDS}; nbad = 0; nbg = 0
    for r in range(R):
        eta = Hy._etas(m, seed0 + r); eta = np.concatenate([eta, -eta])
        accS = 0.0; accG = {k: None for k in KINDS}
        for c0 in range(0, eta.shape[0], Hy.CHUNK):
            et = eta[c0:c0 + Hy.CHUNK]; Z = mu[None, :] + et @ Ch.T
            wgt = np.exp(0.5*ldC + 0.5*np.sum(et*et, axis=1))
            ct = Hy.contour(evS, Z, rho1, keep_path=True)
            accS = accS + np.sum(wgt*ct['Fval']); nbad += int((~ct['ok']).sum())
            for k in KINDS:
                Fg, _, _, okg = Hy.grav_on_paths(evs[k], ct, Z.shape[0])
                v = np.sum(wgt*Fg, axis=-1); accG[k] = v if accG[k] is None else accG[k] + v
                nbg += int((~okg).sum())
        valsS.append(accS/eta.shape[0])
        for k in KINDS:
            valsG[k].append(accG[k]/eta.shape[0])
    valsS = np.array(valsS)
    out = {}
    for k in KINDS:
        vg = np.array(valsG[k])
        (Pg, _), _, _ = evs[k].run(evs[k].grav, np.zeros((1, 7), complex), None)
        out[k] = dict(ljc=True, E=vg.mean(axis=0), err=vg.std(axis=0, ddof=1)/np.sqrt(R), reps=vg, logN=Pg[0] + np.log(np.linalg.det(L)),
                      Es=valsS.mean(), errs=valsS.std(ddof=1)/np.sqrt(R), repsS=valsS, nu=nu, Delta=Delta, J=J, m=NR.mphys(Delta, J),
                      nuf=nuf, rho1=rho1, mq=m, R=R, nbad=nbad, nbad_g=nbg, method='n3_hybrid (odd, P8.build)')
    return out


if __name__ == '__main__':
    for J in [int(a) for a in sys.argv[1:]]:
        for x in XS:
            for Dl in (100.0, 20000.0):
                nu = x*NR.mphys(100.0, J)
                key = (Dl, round(nu, 9))
                if all(key in load(k, J) for k in KINDS):
                    continue
                t = time.time()
                res = point(J, Dl, nu)
                for k in KINDS:
                    db = load(k, J); db[key] = res[k]
                    tmp = fn(k, J) + '.tmp%d' % os.getpid()
                    pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, fn(k, J))
                rel = {k: float(np.max(np.abs(res[k]['err']))/np.max(np.abs(res[k]['E']))) for k in KINDS}
                print(f"J={J} Delta={Dl:.0f} x={x} nu={nu:.4f}: rel err {rel}, scalar {abs(res['C2']['errs']/res['C2']['Es']):.1e}, "
                      f"nbad {res['C2']['nbad']} nbad_g {res['C2']['nbad_g']}, {time.time() - t:.0f}s", flush=True)
