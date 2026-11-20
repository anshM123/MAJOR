"""
hc_points (HYPOTHESES_CLOSE.md): direct all-orders crystal points for the hypothesis-closing runs, state and heavy reference
(Delta_ref = 2e4, same nu) on COMMON RANDOM NUMBERS -- the z2_points / fw_points method unchanged (n3 hybrid contour, nu_f = max(nu, 15)
above 13, nu_f = 13 below; proposal adapted on the state and re-used for the reference, same Sobol seeds; Q_J overflow rescaling QSC;
scalar Q_J by the Legendre recurrence (even J); spin-J Fourier floor fw_frame.patch_S).  NOT applied: fw_frame.patch_conserved_clean
(only needed at J >= 150; it rotates the SVD basis of the conserved <TTO> null space, and the points here must be in the SAME structure
basis as the existing K2 rows they are injected into).
  even J: kinds EE (w0 integrand, C^'_EE) and C3H (o_H q) -- n1_pw / l9_points integrands (z2);
  HC_METHOD=n1: the n1 real re-centred contour (n1_run.point / p8_odd_pw.point integrands), same Sobol seeds for state and reference
  (used at large nu, where n3 is poor), with the same three fixes.
  odd  J: kinds EE (w0.make_build_odd_EE on p8's point-wise copy of p3_odd_action) and C3H (p8 integrand) -- n3_run_odd.OddEval
          evaluators inside z2's CRN block (the scalar evaluator is unchanged).
Caches (one file per STATE, never shared between states; the reference points of a state live in that state's file):
  cache/hc_pts_<kind>_J<J>_D<Dstate>.pkl  keyed (Delta, round(nu, 9)) in z2 / n1_run format; raw: cache/hc_raw_J<J>_D<Dstate>.pkl.
Resumable (finished nodes skipped).
Usage: OMP_NUM_THREADS=1 [HC_METHOD=n3|n1 HC_M=10 HC_MREF=11 HC_ADAPT=1 HC_R=8] python3 hc_points.py J Dstate x1,x2,...      (x = nu/m(Dstate, J))
       OMP_NUM_THREADS=1 python3 hc_points.py J Dstate nu=n1,n2,...                     (explicit nu)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ['Z2_M'] = os.environ.get('HC_M', '10')
os.environ['Z2_R'] = os.environ.get('HC_R', '8')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import z2_points as Z2
import fw_frame as FW

C = os.environ.get('HC_CACHE', Z2.C)          # HC_CACHE: test runs elsewhere (scratch)
STATE = [None]


def fn(kind, J):
    return os.path.join(C, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, STATE[0]))


def rawfn(J):
    return os.path.join(C, 'hc_raw_J%d_D%g.pkl' % (J, STATE[0]))


Z2.fn = fn
_save0 = Z2._save


def _save(path, key, val):
    if os.path.basename(path).startswith('z2') and path.endswith('_raw.pkl'):
        path = rawfn(key[2])
    return _save0(path, key, val)


Z2._save = _save


def setup(odd=False):
    import n1_pw as PW
    FW.patch_qrec_pw(PW, HERE)          # QSC + recurrence (z2.patch_qscale then returns early)
    Z2.patch()
    FW.install_pw(PW)                   # Fourier floor on the point-wise spin engine
    import l1_crystal as Cr
    FW.patch_qrec_tps(Cr, HERE)
    import l5_spin as S
    FW.patch_S(S)
    if odd:
        import y5_n4 as Y5
        Y5.patch()                      # p8_odd_pw.P8.build: kind 'EE' -> w0.make_build_odd_EE
        import n3_flow as F
        import n3_run_odd as N3O
        Orig = F.Evaluator

        def Ev(nu, Delta, J, c, wc, L, kind, nuf=None):
            if kind == 'S':
                return Orig(nu, Delta, J, c, wc, L, 'S', nuf=nuf)
            return N3O.OddEval(nu, Delta, J, c, wc, L, kind, nuf)
        F.Evaluator = Ev


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))



# ------------------------------------------------------------------------------------------------ adaptive proposal (copy of n1_direct.integrate_ad
# that also returns the replicas)
def integrate_ad(fun, c, wc, L, nu, Delta, J, mu, Cv, m=13, R=8, seed0=2000, anti=True):
    import n1_direct as D
    Lc = np.linalg.cholesky(Cv); ldet = 0.5*np.log(np.linalg.det(Cv))
    vals = []; nm = 0
    for r in range(R):
        eta = D._etas(m, seed0 + r)
        if anti:
            eta = np.concatenate([eta, -eta])
        acc = None
        for c0 in range(0, eta.shape[0], D.CHUNK):
            et = eta[c0:c0 + D.CHUNK]; xi = mu + et @ Lc.T
            val, ok, _ = D.evaluate(fun, c, wc, L, xi, nu, Delta, J)
            nm += int(np.sum(~ok))
            w = np.exp(0.5*np.sum(et*et, axis=1) + ldet)
            v = np.sum(np.where(ok, val, 0.0)*w, axis=-1)
            acc = v if acc is None else acc + v
        vals.append(acc/eta.shape[0])
    vals = np.array(vals)
    return vals.mean(axis=0), vals.std(axis=0, ddof=1)/np.sqrt(R), nm, vals

# ------------------------------------------------------------------------------------------------ n1 (real re-centred contour), shared seeds
def n1_point(kind, J, Delta, nu, mq, R):
    """n1_run.point (even: PW.P9.build_spin; odd: p8's P8.build) with the Q_J rescaling (set_qsc at the n1 saddle), Fourier floor,
    Q_J recurrence; Sobol seeds 1000 + r for every call (state and reference share them when mq is equal)."""
    import n1_direct as D
    import n1_pw as PW
    import l5_spin as S
    from n1_run import _both
    c, wc, L = D.full_saddle(nu, Delta, J)
    lq = Z2.set_qsc(Delta, J, c, wc, L)
    if J % 2 == 0:
        mons, sol, _ = S.conserved(Delta, J)
        _, solt, _ = S.conserved(3.0 - Delta, J)
        gfun = lambda zl: PW.P9.build_spin(zl, nu, Delta, J, mons, sol, solt, kind)
        na = sol.shape[1]
    else:
        import p8_odd_pw as P8M
        gfun = lambda zl: P8M.P8.build(zl, nu, Delta, J, kind, 1.0)[:3]
        na = 1
    both = lambda zl: _both(gfun, zl, nu, Delta, J)
    if ADAPT:       # moment-matched Gaussian proposal (n1_direct.adapt, pilot seed 5), RQMC seeds 2000 + r (shared by state and reference)
        mu, Cv = D.adapt(both, c, wc, L, nu, Delta, J, n=8192)
        E, err, nm, reps = integrate_ad(both, c, wc, L, nu, Delta, J, mu, Cv, m=mq, R=R, seed0=2000)
    else:
        E, err, nm = D.integrate_qmc(both, c, wc, L, nu, Delta, J, mq, R, 1.2)
        reps = D.LAST_VALS[0]
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamF, G, _ = gfun(zl)
    logN = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))
    return dict(ljc=True, E=E[:na], err=err[:na], reps=reps[:, :na], Es=E[na, 0], errs=err[na, 0], repsS=reps[:, na, 0], logN=logN,
                lqsc=lq, logNs=None, nmask=nm, nu=nu, Delta=Delta, J=J, m=mphys(Delta, J), mq=mq, R=R,
                method='n1 shared-seed (hc)' + (' adaptive' if ADAPT else ''))


ADAPT = os.environ.get('HC_ADAPT', '1') == '1'
MQREF = int(os.environ.get('HC_MREF', '11'))


def run_n1(J, nus, D, mq, R):
    kinds = Z2.KINDS
    for nu in nus:
        for DD in (D, 20000.0):
            mq_ = mq if DD < 1e4 else MQREF
            key = (DD, round(nu, 9))
            for k in kinds:
                if key in Z2.load(k, J):
                    continue
                t = time.time()
                r = n1_point(k, J, DD, nu, mq_, R)
                Z2._save(fn(k, J), key, r)
                rel = float(np.max(np.abs(r['err']))/np.max(np.abs(r['E'])))
                print(f"  n1 {k} J={J} Delta={DD:g} nu={nu:.4f} x={nu/mphys(DD, J):.4f}: max rel err {rel:.2e}, nmask {r['nmask']}, "
                      f"{time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2])
    STATE[0] = D
    m = mphys(D, J)
    a = sys.argv[3]
    nus = [float(v) for v in a[3:].split(',')] if a.startswith('nu=') else [float(v)*m for v in a.split(',')]
    meth = os.environ.get('HC_METHOD', 'n3')
    setup(odd=(J % 2 == 1))
    t0 = time.time()
    print(f"# hc_points [{meth}] J={J} Delta={D:g} m={m:.3f}: nu = {[round(v, 4) for v in nus]} (x = {[round(v/m, 5) for v in nus]}); "
          f"M={Z2.MQ} R={Z2.RQ} kinds {Z2.KINDS}", flush=True)
    if meth == 'n1':
        run_n1(J, nus, D, Z2.MQ, Z2.RQ)
    else:
        Z2.run(J, nus, [D, 20000.0], Z2.MQ, Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
