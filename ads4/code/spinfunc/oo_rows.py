"""
oo_rows (HYPOTHESES_CLOSE2.md, H_OO at J >= 12): K2 (L1_DEG = 6) crystal rows of the full 3x3 [even(2), odd(1)] action at EVEN J for
the OO block (entry [-1, -1]) -- hc_oo_rows unchanged (kind 'EE' = w0 make_build_odd_EE on p3_odd_action; 'C3H' = p3's o_H q), but
ONE CACHE FILE PER STATE ROW: cache/oo_p3_<kind>_J<J>_D<D>.pkl (key (Delta, round(nu, 9), n0 = 1); the row's state points and its own
references Delta_ref = 2e4 at the same nu).  Nodes already present in the previous agent's hc_p3_<kind>_J<J>.pkl (J = 18: rows 336, 400,
11 nodes) are not recomputed; oo_chain merges both (the reference keys are per-nu and distinct between rows -- asserted there).
Usage: OMP_NUM_THREADS=1 [OO_XS=...] python3 oo_rows.py kind J D [D ...]
"""
import os, sys, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import hc_oo_rows as HR          # XS19 only

C = os.path.join(HERE, 'cache')
DREF = 20000.0
XS19 = HR.XS19
XS = [float(v) for v in os.environ.get('OO_XS', ','.join(map(str, XS19))).split(',')]


def setup(kind, J):
    """hc_oo_rows.setup with the new-shadow cache renamed cache/oo_p2shadow_<kind>_J<J>.pkl (lookup order: project p2 shadow caches,
    hc_p2shadow*, oo_p2shadow*; deterministic shadows, seed 7)."""
    import glob
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    if not getattr(P, '_oo', False):
        P.build = W0.make_build_odd_EE(P)
        P._oo = True
    sfn = os.path.join(C, 'oo_p2shadow_%s_J%d.pkl' % (kind, J))
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J_):
        key = (float(Delta), int(J_), P2.T.DEG)
        for fn in ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'hc_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'oo_p2shadow*.pkl')))):
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


def ofile(kind, J, D):
    return os.path.join(C, 'oo_p3_%s_J%d_D%g.pkl' % (kind, J, D))


if __name__ == '__main__':
    kind = sys.argv[1]; J = int(sys.argv[2]); Ds = [float(v) for v in sys.argv[3:]]
    assert J % 2 == 0 and J >= 4
    P = setup(kind, J)
    old_f = os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind, J))
    old = pickle.load(open(old_f, 'rb')) if os.path.exists(old_f) else {}
    for D in Ds:
        P._fn = lambda kind_, J_, D=D: ofile(kind_, J_, D)
        t0 = time.time(); n = 0
        m = P.mphys(D, J)
        for x in XS:
            nu = x*m
            if (float(D), round(float(nu), 9), 1.0) in old and (DREF, round(float(nu), 9), 1.0) in old:
                continue
            for DD in (D, DREF):
                key = (float(DD), round(float(nu), 9), 1.0)
                if key not in P.load(kind, J):
                    t = time.time()
                    P.get(kind, J, DD, [nu], verbose=False); n += 1
                    print(f"  OO {kind} J={J} Delta={DD:g} (row {D:g}) x={x} nu={nu:.4f}: {time.time() - t:.0f}s", flush=True)
        print(f"done OO {kind} J={J} Delta={D:g} m={m:.1f}: {n} new points in {time.time() - t0:.0f}s", flush=True)
