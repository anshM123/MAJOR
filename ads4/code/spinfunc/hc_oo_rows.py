"""
hc_oo_rows (HYPOTHESES_CLOSE.md, H_OO): K2 (L1_DEG = 6) crystal rows of the parity-odd structure at EVEN J -- the OO block
(X^O_same, 1x1 after the S^0 sum cancels the even-odd interference; PARITY_ODD_ROUTEA §3) -- for Omega^EE:
  kind 'EE'  : w0_ee_common.make_build_odd_EE on p3_odd_action (C^'_EE integrand, full [even(2), odd(1)] structure set at even J);
  kind 'C3H' : p3_odd_action's own C3H (o_H q) integrand (C^_3 unchanged in Omega^EE).
State rows Delta (default 400, 1600: the EE rows of J <= 10 in the w2 chain) at x in HC_XS, each with the heavy reference
Delta_ref = 2e4 at the same nu; n0 = +1 (the n0 = -1 point is identical on the EE/OO entries, PARITY_ODD_ROUTEA §3 [DV]).
Caches (new names only): cache/hc_p3_<kind>_J<J>.pkl (key (Delta, round(nu, 9), n0), p3 format: M0, M2, M4 3x3, logC, logK, ...);
shadows: looked up in the project's p2 shadow caches, new ones -> cache/hc_p2shadow_<kind>_J<J>.pkl (deterministic: seed 7).
Usage: OMP_NUM_THREADS=1 [HC_XS=...] [HC_DS=400,1600] python3 hc_oo_rows.py kind J [J ...]
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
XS19 = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]
XS = [float(v) for v in os.environ.get('HC_XS', ','.join(map(str, XS19))).split(',')]
DS = [float(v) for v in os.environ.get('HC_DS', '400,1600').split(',')]


def setup(kind, J):
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)            # kind 'EE' -> E (x) E; other kinds unchanged
    P._fn = lambda kind_, J_: os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind_, J_))
    sfn = os.path.join(C, 'hc_p2shadow_%s_J%d.pkl' % (kind, J))
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J_):
        key = (float(Delta), int(J_), P2.T.DEG)
        for fn in ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'hc_p2shadow*.pkl')))):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = orig_all(float(Delta), int(J_))
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P2.shadow_odd = shadow_odd
    return P


if __name__ == '__main__':
    kind = sys.argv[1]
    for J in [int(a) for a in sys.argv[2:]]:
        assert J % 2 == 0 and J >= 4
        P = setup(kind, J)
        for D in DS:
            t0 = time.time(); n = 0
            m = P.mphys(D, J)
            for x in XS:
                for DD in (D, DREF):
                    key = (float(DD), round(float(x*m), 9), 1.0)
                    if key not in P.load(kind, J):
                        P.get(kind, J, DD, [x*m], verbose=False); n += 1
            print(f"done OO {kind} J={J} Delta={D:g} m={m:.1f}: {n} new points in {time.time() - t0:.0f}s", flush=True)
