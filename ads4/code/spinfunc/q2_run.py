"""
q2_run (S1 fix, S1_S4_FIXES.md): production driver for the K2 (L1_DEG = 6) crystal points of C^_2 (q^2) and C^_3 (o_H q) at J = 11-20,
on the l9 grid (x = 0.05 ... 0.8) of the rows Delta = 200, 400, 800, 1600, each with its heavy reference Delta_ref = 2e4 at the SAME nu.
  even J : l9_points.action / scalar6 (2x2 parity-even conserved basis; identical integrand to l9 C2 / C3H)
  odd  J : p3_odd_action.action (1x1 parity-odd conserved structure, +- only; identical integrand to p3 C2 / C3H)
The project modules are imported read-only; ONLY their cache file names are redirected to q2_* files (per-process tag Q2_TAG, so that two
processes never write the same pickle):
  cache/q2_pts_<kind>_J<J>_k2<tag>.pkl, cache/q2_shadow_k2<tag>.pkl  (even);  cache/q2_p3_<kind>_J<J><tag>.pkl, cache/q2_p2shadow<tag>.pkl (odd)
Usage: OMP_NUM_THREADS=1 Q2_TAG=a python3 q2_run.py even|odd "12,14" "400,800" [x-list]
"""
import os, sys, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

TAG = os.environ.get('Q2_TAG', '')
XS = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
DREF = 20000.0
C = os.path.join(HERE, 'cache')


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    P._fn = lambda kind, J: os.path.join(C, 'q2_pts_%s_J%d_k2%s.pkl' % (kind, J, TAG))
    sfn = os.path.join(C, 'q2_shadow_k2%s.pkl' % TAG)

    def shadow(Delta, J):
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        key = (float(Delta), int(J))
        if key not in db:
            r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
            db[key] = (r['b'], r['logK'], r['resid'])
            pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P.shadow = shadow
    return P


def setup_odd():
    import p2_shadow_odd as P2
    P2._DB = os.path.join(C, 'q2_p2shadow%s.pkl' % TAG)
    import p3_odd_action as P
    P._fn = lambda kind, J: os.path.join(C, 'q2_p3_%s_J%d%s.pkl' % (kind, J, TAG))
    return P


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


if __name__ == '__main__':
    par = sys.argv[1]
    Js = [int(a) for a in sys.argv[2].split(',')]
    Ds = [float(a) for a in sys.argv[3].split(',')]
    xs = [float(a) for a in sys.argv[4].split(',')] if len(sys.argv) > 4 else XS
    P = setup_even() if par == 'even' else setup_odd()
    for D in Ds:
        for J in Js:
            t0 = time.time()
            m = mphys(D, J)
            for x in xs:
                for kind in ('C2', 'C3H'):
                    P.get(kind, J, D, [x*m], verbose=False)
                    P.get(kind, J, DREF, [x*m], verbose=False)
            print(f"done {par} J={J} Delta={D:.0f} ({len(xs)} x, 2 kinds, state + ref): {time.time() - t0:.0f}s", flush=True)
