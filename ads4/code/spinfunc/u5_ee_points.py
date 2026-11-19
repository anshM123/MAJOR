"""
u5 (EE_FUNCTIONAL.md §5): K2 crystal points of C^_EE = Psi_2[EE] (j' = 0, f = v'-u', celestial j' = 0) with the l9/l6 engine, i.e. the
same integrand as l9_points kind 'C2' with the photon factor
    q = (2 H20 H40 + V042^2 H24)/(X20 X40)       replaced by      E = V042^2 (H24 - 2 V_{2,04} V_{4,02})/(X20 X40)
(u1: covariant form of V0^4 (z2'.z4')^2; E = q exactly at the crystal point, so M0 must agree with C2's M0 -- the built-in check).
l9_points is imported read-only; only its build_spin is wrapped (monkey-patch inside this process) and its cache name redirected to
cache/u5_pts_EE_J<J>_k2<tag>.pkl.  Keys (Delta, nu) are copied from the q^2 records (q2_pts_C2_J<J>_k2*.pkl / l8 records) so that the q3
chain can use EE rows in place of C2 rows one-for-one.
Usage: OMP_NUM_THREADS=1 U5_TAG=a python3 u5_ee_points.py J "D1,D2" [stride]     (Delta = 20000 refs are done with their states)
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l5_spin as S
import l9_points as P

TAG = os.environ.get('U5_TAG', '')
C = os.path.join(HERE, 'cache')
_orig_build = P.build_spin


def build_spin_EE(zv, nu, Delta, J, mons, sol, solt, kind):
    if kind != 'EE':
        return _orig_build(zv, nu, Delta, J, mons, sol, solt, kind)
    # patch the q line through a local copy of l9_points.build_spin: replace q by E
    import l1_crystal as Cr
    import l6_spin_action as A6
    lamF, G0, parts = Cr.build(zv, nu, Delta, 0, n0=1.0, parts=True)
    g = parts['g']
    Pg = A6.geom(zv, 1.0)
    x0, x1, x2, x3, x4, e, z0 = (Pg[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    msq = A6.msq
    m12, m34, m13, m24, m14, m23 = msq(x1, x2), msq(x3, x4), msq(x1, x3), msq(x2, x4), msq(x1, x4), msq(x2, x3)
    m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
    v042 = S.Vfun(x0, z0, x4, x2)
    up = T.div(T.mul(m12, m34), T.mul(m13, m24)); vp = T.div(T.mul(m14, m23), T.mul(m13, m24))
    lf = T.log(vp - up)
    Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] + lf - 4*T.log(v042)
          - ((10 - J)/2)*T.log(m12) - (J/2)*(T.log(m15) + T.log(m25))
          - ((7 - J)/2)*T.log(m34) - ((3 + J)/2)*(T.log(m35) + T.log(m45)))
    W2, WJ = S.null_circle(2), S.null_circle(J)
    k2, kJ = S.kappa(2), S.kappa(J)
    x10, x30 = Cr.sub(x1, x0), Cr.sub(x3, x0)
    F = S.tto_samples(x1, x10, x2, W2, e, WJ, J, mons, sol)
    Gt = S.tto_samples(x3, x30, x4, W2, e, WJ, J, mons, solt)
    F = S.fourier(S.fourier(F, 1, 2), 2, J); Gt = S.fourier(S.fourier(Gt, 1, 2), 2, J)
    H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
    H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
    V204 = S.lin_samples([S.Vfun(x2, S.BASIS[m], x0, x4) for m in range(3)], W2)
    V402 = S.lin_samples([S.Vfun(x4, S.BASIS[m], x0, x2) for m in range(3)], W2)
    inv = T.recip(T.mul(S.Xfun(x2, x0), S.Xfun(x4, x0)))
    Eb = S.bm(S.bm(S.bm(v042, v042)[None, None], H24 - 2*S.bm(V204[:, None], V402[None, :])), inv)
    D = S.bm(Eb, Eb)
    D = S.fourier(S.fourier(D, 0, 2), 1, 2)
    na, nc = F.shape[0], Gt.shape[0]
    pref = np.zeros((na, nc, T.N), complex)
    for a in range(na):
        Em = np.einsum('i,ijkn->jkn', 1/k2, S.bm(D[:, :, None, :], F[a][::-1][:, None, :, :]))
        for c in range(nc):
            Gr = Gt[c][::-1, ::-1]
            pref[a, c] = np.einsum('j,k,jkn->n', 1/k2, 1/kJ, S.bm(Em, Gr))
    return lamF, Gs, pref


P.build_spin = build_spin_EE
P._fn = lambda kind, J: os.path.join(C, 'u5_pts_%s_J%d_k2%s.pkl' % (kind, J, TAG))
_sfn = os.path.join(C, 'u5_shadow_k2%s.pkl' % TAG)


def _shadow(Delta, J):
    # shadow matrices are structure-independent: reuse the q2/l9 K2 caches when present, else compute and store in u5's own cache
    for fn in sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) + [os.path.join(C, 'l9_shadow_k2.pkl'), _sfn]:
        if os.path.exists(fn):
            db = pickle.load(open(fn, 'rb'))
            if (float(Delta), int(J)) in db:
                return db[(float(Delta), int(J))]
    import l5_shadow_matrix as SM
    r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
    db = pickle.load(open(_sfn, 'rb')) if os.path.exists(_sfn) else {}
    db[(float(Delta), int(J))] = (r['b'], r['logK'], r['resid'])
    pickle.dump(db, open(_sfn + '.tmp', 'wb')); os.replace(_sfn + '.tmp', _sfn)
    return db[(float(Delta), int(J))]


P.shadow = _shadow


def q2_keys(J):
    R = {}
    for fn in sorted(glob.glob(os.path.join(C, 'q2_pts_C2_J%d_k2*.pkl' % J))):
        R.update(pickle.load(open(fn, 'rb')))
    if not R:
        import l8_analysis as L8
        R = L8.records(J)
    return R


def load_all(J):
    R = {}
    for fn in sorted(glob.glob(os.path.join(C, 'u5_pts_EE_J%d_k2*.pkl' % J))):
        R.update(pickle.load(open(fn, 'rb')))
    return R


if __name__ == '__main__' and sys.argv[1] != 'diag':
    J = int(sys.argv[1]); Ds = [float(v) for v in sys.argv[2].split(',')]
    stride = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    R = q2_keys(J)
    for D in Ds:
        nus = sorted(k[1] for k in R if k[0] == D)[::stride]
        t0 = time.time()
        for nu in nus:
            P.get('EE', J, D, [nu], verbose=False)
            P.get('EE', J, 20000.0, [nu], verbose=False)
        Rn = load_all(J)
        # raw M0 differ by nu-only factors (saddle curvature; cancel against the Delta_ref row at the same nu): compare state/ref
        def nrm(RR, nu):
            a, b = RR[(D, round(nu, 9))]['M0'], RR[(20000.0, round(nu, 9))]['M0']
            return a/np.trace(b)
        dev = max(np.max(np.abs(nrm(Rn, nu) - nrm(R, nu)))/np.max(np.abs(nrm(R, nu))) for nu in nus if (20000.0, round(nu, 9)) in R)
        print(f"done EE J={J} Delta={D:.0f}: {len(nus)} nu (+refs) in {time.time() - t0:.0f}s; max |M0/tr M0_ref (EE - q2)| rel = {dev:.1e}",
              flush=True)


def diag_scaling(Js=(4, 8, 18), Ds=(400.0, 1600.0, 6400.0), xs=(0.15, 0.3)):
    """M0 ratio EE/q2 (raw, same (Delta, nu); logC, logK are structure-independent) vs J and m: tests whether EE = q^2 at the crystal
    holds for the engine's full (null-circle) polarisation contraction in the limit m -> infinity at fixed J."""
    out = []
    for J in Js:
        for D in Ds:
            m = P.mphys(D, J)
            for x in xs:
                t0 = time.time()
                rq = P.action(x*m, D, J, 'C2'); re = P.action(x*m, D, J, 'EE')
                r0 = re['M'][:, :, 0]/rq['M'][:, :, 0]
                k1 = (re['M'][:, :, 0] + re['M'][:, :, 1])/(rq['M'][:, :, 0] + rq['M'][:, :, 1])
                s = (f"J={J:2d} Delta={D:6.0f} m={m:8.1f} x={x}: J^2/m={J*J/m:.3f}  M0 ratio EE/q2 = " +
                     " ".join(f"{v.real:+.4f}{v.imag:+.4f}i" for v in r0.flatten()) +
                     "   (M0+M2) ratio = " + " ".join(f"{v.real:+.4f}{v.imag:+.4f}i" for v in k1.flatten()) + f"  [{time.time() - t0:.0f}s]")
                print(s, flush=True); out.append(s)
    open(os.path.join(C, 'u5_diag_scaling.txt'), 'w').write("\n".join(out) + "\n")


if __name__ == '__main__' and sys.argv[1] == 'diag':
    diag_scaling()
