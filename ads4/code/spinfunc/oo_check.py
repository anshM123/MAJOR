"""
oo_check (HYPOTHESES_CLOSE2.md): validation of the direct OO points (oo_points.py) node by node:
  (a) scalar check against the exact Mellin object (z3.scalar_check: d_D - d_ref must vanish; Delta c2 with jackknife error) -- the
      contour test of every node (the scalar is integrated on the same samples as the 3x3 graviton action);
  (b) relative RQMC errors of the OO element (state, reference, CRN ratio; dir and scalar control variate);
  (c) where a K2 OO node exists at the same nu: direct/K2 of the state/reference ratio of the OO element (normalisation and
      shadow-basis consistency; K2 is expected to fail at small nu, cf. the EE ++ at x = 0.05);
  (d) the order-0 coupling test of the OO element, C3H/C2 x (2 - y) on the direct points (state/ref ratio; = the ++ value ~2 at LO).
Usage: python3 oo_check.py J:D [J:D ...]  -> cache/oo_check.txt (appended per call: OO_CHECK_TAG)
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import z3_nodes as Z3

C = os.path.join(HERE, 'cache')
DREF = 20000.0
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def load(path):
    return pickle.load(open(path, 'rb')) if os.path.exists(path) else {}


def k2oo(kind, J):
    import glob
    R = {}
    for f in [os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind, J))] + sorted(glob.glob(os.path.join(C, 'oo_p3_%s_J%d_D*.pkl' % (kind, J)))):
        R.update(load(f))
    return {(k[0], round(k[1], 6)): v for k, v in R.items() if k[2] == 1.0}


if __name__ == '__main__':
    import glob
    for spec in sys.argv[1:]:
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        m = Z3.mphys(D, J)
        P = {kd: load(os.path.join(C, 'oo_pts_%s_J%d_D%g.pkl' % (kd, J, D))) for kd in ('EE', 'C3H')}
        K2 = {kd: k2oo(kd, J) for kd in ('EE', 'C3H')}
        hc = os.path.join(C, 'hc_exact.pkl'); db = load(hc)
        Z3.EXPK = hc if ((round(D, 9), J) in db and (round(DREF, 9), J) in db) else os.path.join(C, 'oo_exact.pkl')
        pr(f"## J={J} Delta={D:g} (m = {m:.3f}); exact scalar cache {os.path.basename(Z3.EXPK)}")
        nus = sorted(set(k[1] for k in P['EE'] if k[0] == D and (DREF, k[1]) in P['EE']))
        for nu in nus:
            dd = {DD: P['EE'][(DD, nu)] for DD in (D, DREF)}
            try:
                o = Z3.scalar_check(J, nu, dd, D)
                sc = f"scalar d-d_ref Re {o['d'].real:+.4f}+-{o['err'].real:.4f}; Dc2 {o['dc2']:+7.0f}+-{o['ec2']:5.0f} ({abs(o['dc2'])/max(o['ec2'], 1e-9):.1f}s)"
            except Exception as ex:
                o = None; sc = f"scalar check failed ({ex})"
            s = f"  nu={nu:7.3f} x={nu/m:.4f}: {sc}"
            rat = {}
            for kd in ('EE', 'C3H'):
                st, rf = P[kd].get((D, nu)), P[kd].get((DREF, nu))
                if st is None or rf is None:
                    continue
                R = len(st['reps'])
                f1 = lambda s_: st['reps'][s_].mean(0)[-1, -1]/rf['reps'][s_].mean(0)[-1, -1]
                v0, _, e = Z3._jk(f1, R)
                f2 = lambda s_: (st['reps'][s_].mean(0)[-1, -1]/st['repsS'][s_].mean())/(rf['reps'][s_].mean(0)[-1, -1]/rf['repsS'][s_].mean())
                _, _, e2 = Z3._jk(f2, R)
                rel_st = abs(st['err'][-1, -1]/st['E'][-1, -1]); rel_rf = abs(rf['err'][-1, -1]/rf['E'][-1, -1])
                # full ratio in K2 units: exp(logN_st - logN_ref) E_st/E_ref vs exp(logC_st - logC_ref) M_st/M_ref x (b_st/b_ref)
                full = np.exp(st['logN'] - rf['logN'])*v0
                rat[kd] = full
                s += f" | {kd}: OO rel err st {rel_st:.3f} ref {rel_rf:.3f} ratio(dir) {e/abs(v0):.3f} (cv) {e2/abs(v0):.3f}"
                k2s, k2r = K2[kd].get((D, round(nu, 6))), K2[kd].get((DREF, round(nu, 6)))
                if k2s is not None and k2r is not None:
                    bs = k2s['M0'][-1, -1]; br = k2r['M0'][-1, -1]
                    tk = np.exp(k2s['logC'] - k2r['logC'])*sum(k2s[o_][-1, -1] for o_ in ('M0', 'M2', 'M4'))/sum(k2r[o_][-1, -1] for o_ in ('M0', 'M2', 'M4'))
                    sh = {}
                    for DD, rr in ((D, k2s), (DREF, k2r)):
                        import oo_chain as S
                        (b, logK, _), _ = S.odd_shadow(DD, J, rr['shadow_resid'])
                        sh[DD] = np.asarray(b)[-1, -1]
                    dirK = full/(sh[D]/sh[DREF])
                    s += f" direct/K2 {complex(dirK/tk).real:+.4f}{complex(dirK/tk).imag:+.4f}i"
            if 'EE' in rat and 'C3H' in rat:
                y = (nu/m)**2
                s += f" | LO coupling C3H/C2 (2-y) (state/ref, m^2/m_r^2 scaled) {complex(rat['C3H']/rat['EE']*(m*m/Z3.mphys(DREF, J)**2)*(2 - y)).real:+.3f}"
            pr(s)
    tag = os.environ.get('OO_CHECK_TAG', '')
    with open(os.path.join(C, 'oo_check%s.txt' % tag), 'a') as f:
        f.write("\n".join(lines) + "\n")
