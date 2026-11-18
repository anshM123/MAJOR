"""
n1_run (doc 24): non-perturbative (direct, RQMC) graviton crystal integrals for the l9 kinds C2 / C3H at the l9 grid points
(Delta = 100 row, x = nu/m in N1_XS) and their heavy references (Delta_ref = 2e4, same nu), on the re-centred log-coordinate contour of
n1_direct (path-continued logs, wall/winding masks).  Stored per point (resumable; cache/n1_pts_<kind>_J<J>.pkl, key (Delta, nu)):
  E[a', c]   : prefactor integral in the conserved-structure basis (same convention as l6/l9 'Ser' up to the normalisation below)
  err[a', c] : RQMC standard error (R scramblings)
  logN       : log[ exp(lamF + G + logJac)(w_c) det L_w ]           (direct normalisation)
  logC9      : log[ exp(lamF + G)(z*) det L_crystal ]                (l9 normalisation, logC of l9_points.action)
  so that M (l9 units, all orders) = exp(logN - logC9) E @ inv(b)^T,   b from the K2 shadow cache (cache/l9_shadow_k2.pkl).
  nmask, wind diagnostics; the scalar crystal (pref form, same contour and samples) as Es (for graviton/scalar ratios).
Usage: OMP_NUM_THREADS=1 python3 n1_run.py kind J [J ...]     env N1_XS (default 0.25,0.3,0.4,0.5,0.6,0.7,0.8), N1_M (log2 pts, 13),
       N1_R (scramblings, 8), N1_DS (default 100,20000)
"""
import os, sys, pickle, time, warnings
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
warnings.filterwarnings('ignore')
import numpy as np
import n1_direct as D
import n1_pw as PW
import l5_spin as S
import l1_tps as T
import l1_crystal as Cr

XS = [float(v) for v in os.environ.get('N1_XS', '0.25,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]
MQ = int(os.environ.get('N1_M', '14'))
MQREF = int(os.environ.get('N1_MREF', '12'))
RQ = int(os.environ.get('N1_R', '8'))
DS = [float(v) for v in os.environ.get('N1_DS', '100,20000').split(',')]
DROW = float(os.environ.get('N1_DROW', '100'))          # the l9 row whose x-grid defines nu = x m(DROW, J)


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def fn(kind, J):
    return os.path.join(HERE, 'cache', 'n1_pts_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def point(kind, J, Delta, nu):
    mons, sol, _ = S.conserved(Delta, J)
    _, solt, _ = S.conserved(3.0 - Delta, J)
    c, wc, L = D.full_saddle(nu, Delta, J)
    gfun = lambda zl: PW.P9.build_spin(zl, nu, Delta, J, mons, sol, solt, kind)
    sfun = lambda zl: (lambda r: (r[0], r[1], r[2][None, None, :]))(PW.Cr.build(zl, nu, Delta, J, pref=True))
    both = lambda zl: _both(gfun, zl, nu, Delta, J)
    E, err, nm = D.integrate_qmc(both, c, wc, L, nu, Delta, J, MQ if Delta < 1e4 else MQREF, RQ, 1.2)
    # normalisations
    PW.set_batch(1)
    zl, lj = D.zmap([np.array([wc[i]]) for i in range(7)], c, pw=True)
    lamF, G, _ = gfun(zl)
    logN = lamF[0] + G[0] + lj[0] + np.log(np.linalg.det(L))
    zs, Lc = D.saddle_L(nu, Delta)
    PW.set_batch(1)
    lamF0, G0, _ = gfun([np.array([zs[i]]) for i in range(7)])
    logC9 = lamF0[0] + G0[0] + np.log(np.linalg.det(Lc))
    na = sol.shape[1]
    return dict(ljc=True, E=E[:na], err=err[:na], Es=E[na, 0], errs=err[na, 0], logN=logN, logC9=logC9, nmask=nm, wc=wc, m=mphys(Delta, J),
                reps=D.LAST_VALS[0])


def _both(gfun, zl, nu, Delta, J):
    """graviton prefactor (na x nc) stacked with the scalar Q_J (as an extra row), with the graviton exponent; the scalar row is
    multiplied by exp(scalar exponent - graviton exponent) so that both are integrated on the same samples."""
    lamF, G, pref = gfun(zl)
    lamS, GS, QS = PW.Cr.build(zl, nu, Delta, J, pref=True)
    rel = np.exp((lamS - lamS[0] + GS - GS[0]) - (lamF - lamF[0] + G - G[0]))
    na, nc = pref.shape[:2]
    out = np.zeros((na + 1, nc) + pref.shape[2:], complex)
    out[:na] = pref
    out[na, 0] = QS*rel*np.exp((lamS[0] + GS[0]) - (lamF[0] + G[0]))
    return lamF, G, out


if __name__ == '__main__':
    kind = sys.argv[1]
    for J in [int(a) for a in sys.argv[2:]]:
        for Dl in DS:
            for x in XS:
                nu = x*mphys(DROW, J)
                key = (Dl, round(nu, 9))
                if key in load(kind, J):
                    continue
                t = time.time()
                r = point(kind, J, Dl, nu)
                db = load(kind, J); db[key] = r
                tmp = fn(kind, J) + '.tmp%d' % os.getpid()
                pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, fn(kind, J))
                rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
                print(f"{kind} J={J} Delta={Dl:.0f} x={x} nu={nu:.4f}: max rel err {rel:.1e}, scalar rel err {abs(r['errs']/r['Es']):.1e}, "
                      f"nmask {r['nmask']}, {time.time() - t:.0f}s", flush=True)
