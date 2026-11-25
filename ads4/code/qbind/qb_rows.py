"""
qb_rows (QBIND_COROLLARY.md §1; copy of neargap2/ng_rows.py, cache redirected to qbind/cache): K2 crystal rows (L1_DEG = 6; state + heavy reference Delta_ref = 2e4 at the same nu) at the OWN mass of a
near-gap state, e.g. the binding state (J = 5, +-, parity-odd, twist edge Delta = 1005 at Delta_g = 1000), on a dense x grid.
It is w1_rows.py (EE_NEARGAP.md), imported read-only and run unchanged, with every cache path redirected to neargap2/cache:
  odd : cache/w1_p3_<kind>_J<J>_ng.pkl          (p3_odd_action engine, 1x1 parity-odd +-)
  even: cache/w1_pts_<kind>_J<J>_k2_ng.pkl      (l9_points engine, w0 EE integrand, 2x2)
  shadows: looked up in the project caches (read-only) and in neargap2/cache; new ones -> cache/ng_shadow*.pkl
Usage: OMP_NUM_THREADS=1 python3 ng_rows.py odd|even KINDS J Delta x1,x2,...|dense
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import qb_guard  # noqa: F401
import warnings
warnings.filterwarnings('ignore')
import numpy as np
os.chdir(SF)
import w1_rows as W1

NGC = os.path.join(HERE, 'cache')
NG2 = os.path.join(HERE, '..', 'neargap2', 'cache')       # read-only
SFC = os.path.join(SF, 'cache')
DREF = 20000.0
DENSE = [round(0.025*i, 4) for i in range(1, 37)]          # x = 0.025 ... 0.9


def setup_odd():
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    P._fn = lambda kind, J: os.path.join(NGC, 'w1_p3_%s_J%d_qb.pkl' % (kind, J))
    P2._DB = os.path.join(NGC, 'qb_p2shadow_db.pkl')
    sfn = os.path.join(NGC, 'qb_p2shadow.pkl')
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(SFC, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(SFC, 'w1_p2shadow*.pkl'))) + [os.path.join(NG2, 'ng_p2shadow.pkl'), sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = orig_all(float(Delta), int(J))
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P2.shadow_odd = shadow_odd
    return P


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(NGC, 'w1_pts_%s_J%d_k2_qb.pkl' % (kind, J))
    sfn = os.path.join(NGC, 'qb_shadow_k2.pkl')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(SFC, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(SFC, 'u5_shadow_k2.pkl'), os.path.join(SFC, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(SFC, 'w1_shadow_k2*.pkl'))) + [os.path.join(NG2, 'ng_shadow_k2.pkl'), sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P.shadow = shadow
    # l9_points may keep other caches; make sure every module-level cache path it uses points to NGC
    return P


def have(R, D, nu):
    return W1.have(R, D, nu)


if __name__ == '__main__':
    par, kinds, J, D = sys.argv[1], sys.argv[2].split(','), int(sys.argv[3]), float(sys.argv[4])
    xs = DENSE if sys.argv[5] == 'dense' else [float(v) for v in sys.argv[5].split(',')]
    P = setup_even() if par == 'even' else setup_odd()
    m = W1.mphys(D, J)
    t0 = time.time(); n = 0
    for kind in kinds:
        for x in xs:
            for DD in (D, DREF):
                fn = P._fn(kind, J)
                R = {(k[0], k[1]): v for k, v in (pickle.load(open(fn, 'rb')).items() if os.path.exists(fn) else [])}
                if not have(R, DD, x*m):
                    t1 = time.time()
                    P.get(kind, J, DD, [x*m], verbose=False); n += 1
                    print(f"  {kind} J={J} D={DD:g} x={x:.4f} nu={x*m:.3f}: {time.time() - t1:.1f}s", flush=True)
    print(f"done {par} {','.join(kinds)} J={J} Delta={D:.6g} m={m:.3f} J^2/m={J*J/m:.3f}: {n} new points in {time.time() - t0:.0f}s",
          flush=True)
