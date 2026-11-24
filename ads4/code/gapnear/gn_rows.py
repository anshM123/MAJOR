"""
gn_rows (GAP_NEARGAP.md): own-mass K2 crystal rows (L1_DEG = 6; state + heavy reference Delta_ref = 2e4 at the same nu) of a near-gap
state at Delta_g = 2000 / 5000 (Delta = Delta_g + J; light J = 0, 2 at Delta_L = 0.99 Delta_g), on the dense x grid 0.025 ... 0.9.
Odd J = gapscale/gs_rows.py (= neargap2/ng_rows.py's odd branch); even J = neargap2/ng_rows.py's even branch; both re-typed with every
cache path redirected to gapnear/cache/rows/.  Engines (spinfunc/w1_rows, p3_odd_action, p2_shadow_odd, l9_points, l5_shadow_matrix,
w0_ee_common) are imported read-only and run unchanged.  Shadows are looked up read-only in spinfunc/cache, neargap2/cache,
gapscale/cache/rows; new ones go to gapnear/cache/rows.  gn_guard refuses any write outside gapnear/cache.
Usage: OMP_NUM_THREADS=1 python3 gn_rows.py odd|even KINDS J Delta x1,x2,...|dense
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gn_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc')
NGC = os.path.join(HERE, '..', 'neargap2', 'cache')
GSR = os.path.join(HERE, '..', 'gapscale', 'cache', 'rows')
ROWS = os.path.join(HERE, 'cache', 'rows')
os.makedirs(ROWS, exist_ok=True)
sys.path.insert(0, SF)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
os.chdir(SF)
import w1_rows as W1

SFC = os.path.join(SF, 'cache')
DREF = 20000.0
DENSE = [round(0.025*i, 4) for i in range(1, 37)]          # x = 0.025 ... 0.9


def setup_odd():
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    P._fn = lambda kind, J: os.path.join(ROWS, 'w1_p3_%s_J%d_gn.pkl' % (kind, J))
    P2._DB = os.path.join(ROWS, 'gn_p2shadow_db.pkl')
    sfn = os.path.join(ROWS, 'gn_p2shadow.pkl')
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(SFC, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(SFC, 'w1_p2shadow*.pkl'))) +
                   [os.path.join(NGC, 'ng_p2shadow.pkl'), os.path.join(GSR, 'gs_p2shadow.pkl'), sfn]):
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
    P._fn = lambda kind, J: os.path.join(ROWS, 'w1_pts_%s_J%d_k2_gn.pkl' % (kind, J))
    sfn = os.path.join(ROWS, 'gn_shadow_k2.pkl')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(SFC, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(SFC, 'u5_shadow_k2.pkl'), os.path.join(SFC, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(SFC, 'w1_shadow_k2*.pkl'))) + [os.path.join(NGC, 'ng_shadow_k2.pkl'), sfn]):
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
    return P


if __name__ == '__main__':
    par, kinds, J, D = sys.argv[1], sys.argv[2].split(','), int(sys.argv[3]), float(sys.argv[4])
    xs = DENSE if sys.argv[5] == 'dense' else [float(v) for v in sys.argv[5].split(',')]
    assert (par == 'odd') == (J % 2 == 1)
    P = setup_even() if par == 'even' else setup_odd()
    m = W1.mphys(D, J)
    t0 = time.time(); n = 0
    for kind in kinds:
        for x in xs:
            for DD in (D, DREF):
                fn = P._fn(kind, J)
                R = {(k[0], k[1]): v for k, v in (pickle.load(open(fn, 'rb')).items() if os.path.exists(fn) else [])}
                if not W1.have(R, DD, x*m):
                    t1 = time.time()
                    P.get(kind, J, DD, [x*m], verbose=False); n += 1
                    print(f"  {kind} J={J} D={DD:g} x={x:.4f} nu={x*m:.3f}: {time.time() - t1:.1f}s", flush=True)
    print(f"done {par} {','.join(kinds)} J={J} Delta={D:.6g} m={m:.3f} J^2/m={J*J/m:.3f}: {n} new points in {time.time() - t0:.0f}s",
          flush=True)
