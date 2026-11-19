"""
r1_run (LARGEJ_OFFDIAG.md): crystal points for the FIXED-b-hat ladder at Delta_g = 300 (Lam = 298.5):
    (Delta, J) = (336, 18) [exists: q2 tag _h + q4],  (672, 36),  (1344, 72)        b-hat = (2J+1) Lam / m ~ 33 for all three
and, optionally, any other (Delta, J).  Same engines as q2_run (K2 Laplace, l9_points) and q4_direct (non-perturbative RQMC, n1_run.point);
project modules are imported read-only, only their cache file names are redirected:
    K2 state/ref points   cache/r1_pts_<kind>_J<J>_k2.pkl      shadow   cache/r1_shadow_k2.pkl
    reference basis (C2 at Delta_ref, x_ref = 0.3 ... 0.8)       same files
    direct (RQMC)         cache/r1_n1_<kind>_J<J>.pkl  key (Delta, nu)
Nodes are given in nu (not x) so that the ladder shares the nu-grid of the (336, 18) row: nu_i = x_i * m(336, 18).
Usage:  OMP_NUM_THREADS=1 python3 r1_run.py k2 J Delta nu1,nu2,...        (C2 + C3H, state + ref)
        OMP_NUM_THREADS=1 python3 r1_run.py refbasis J
        OMP_NUM_THREADS=1 python3 r1_run.py direct kind J Delta nu1,nu2,...  (state + ref)
"""
import os, sys, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    P._fn = lambda kind, J: os.path.join(C, 'r1_pts_%s_J%d_k2.pkl' % (kind, J))
    sfn = os.path.join(C, 'r1_shadow_k2.pkl')

    def shadow(Delta, J):
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        key = (float(Delta), int(J))
        if key not in db:
            # reuse any shadow entry already computed by q2 (same engine, same seed)
            import glob
            for fn in glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl')):
                d = pickle.load(open(fn, 'rb'))
                if key in d:
                    db[key] = d[key]; break
        if key not in db:
            t = time.time()
            r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
            db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
            db[key] = (r['b'], r['logK'], r['resid'])
            print(f"shadow {key}: {time.time() - t:.0f}s, resid {r['resid']}", flush=True)
        pickle.dump(db, open(sfn + '.tmp%d' % os.getpid(), 'wb')); os.replace(sfn + '.tmp%d' % os.getpid(), sfn)
        return db[key]
    P.shadow = shadow
    return P


def dfn(kind, J):
    return os.path.join(C, 'r1_n1_%s_J%d.pkl' % (kind, J))


def dload(kind, J):
    f = dfn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def point_g(kind, J, Delta, nu):
    """n1_run.point WITHOUT the stacked scalar row (its relative exponent overflows at the heavy reference for J >= 36, which masked
    every sample: r1_t2.log); identical graviton integrand, contour, samples and normalisations; Es := nan."""
    import n1_run as R
    D, PW, S = R.D, R.PW, R.S
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    c, wc, L = D.full_saddle(nu, Delta, J)
    gfun = lambda zl: PW.P9.build_spin(zl, nu, Delta, J, mons, sol, solt, kind)
    E, err, nm = D.integrate_qmc(gfun, c, wc, L, nu, Delta, J, R.MQ if Delta < 1e4 else R.MQREF, R.RQ, 1.2)
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamF, G, _ = gfun(zl)
    logN = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))
    zs, Lc = D.saddle_L(nu, Delta)
    PW.set_batch(1)
    lamF0, G0, _ = gfun([np.array([zs[i]]) for i in range(7)])
    logC9 = lamF0[0] + G0[0] + np.log(np.linalg.det(Lc))
    na = sol.shape[1]
    return dict(ljc=True, E=E[:na], err=err[:na], Es=np.nan, errs=np.nan, logN=logN, logC9=logC9, nmask=nm, wc=wc, m=R.mphys(Delta, J),
                reps=D.LAST_VALS[0][:, :na])


def direct(kind, J, Delta, nus, Ds=None):
    import n1_run as R
    R.point = point_g
    for D in (Ds or (Delta, DREF)):
        for nu in nus:
            key = (D, round(nu, 9))
            if key in dload(kind, J):
                continue
            t = time.time()
            r = R.point(kind, J, D, nu)
            db = dload(kind, J); db[key] = r
            pickle.dump(db, open(dfn(kind, J) + '.tmp%d' % os.getpid(), 'wb')); os.replace(dfn(kind, J) + '.tmp%d' % os.getpid(), dfn(kind, J))
            rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
            print(f"direct {kind} J={J} Delta={D:.0f} nu={nu:.3f} (x={nu/mphys(D, J):.4f}): max rel err {rel:.1e}, nmask {r['nmask']}, "
                  f"{time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'k2':
        P = setup_even()
        J = int(sys.argv[2]); D = float(sys.argv[3]); nus = [float(v) for v in sys.argv[4].split(',')]
        for nu in nus:
            t0 = time.time()
            for kind in ('C2', 'C3H'):
                P.get(kind, J, D, [nu], verbose=False)
                P.get(kind, J, DREF, [nu], verbose=False)
            print(f"k2 J={J} Delta={D:.0f} nu={nu:.3f} (x={nu/mphys(D, J):.4f}): {time.time() - t0:.0f}s", flush=True)
    elif mode == 'k2c2':                       # C2 state only (basis points, x >= 0.15)
        P = setup_even()
        J = int(sys.argv[2]); D = float(sys.argv[3]); nus = [float(v) for v in sys.argv[4].split(',')]
        for nu in nus:
            t0 = time.time()
            P.get('C2', J, D, [nu], verbose=False)
            print(f"k2c2 J={J} Delta={D:.0f} nu={nu:.3f} (x={nu/mphys(D, J):.4f}): {time.time() - t0:.0f}s", flush=True)
    elif mode == 'refbasis':
        P = setup_even()
        J = int(sys.argv[2]); mr = mphys(DREF, J)
        for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
            t0 = time.time()
            P.get('C2', J, DREF, [x*mr], verbose=False)
            print(f"refbasis J={J} x_ref={x}: {time.time() - t0:.0f}s", flush=True)
    elif mode == 'direct':
        direct(sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), [float(v) for v in sys.argv[5].split(',')])
