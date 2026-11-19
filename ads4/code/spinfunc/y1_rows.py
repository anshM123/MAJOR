"""
y1 (EE_DG100.md): K2 (L1_DEG = 6) crystal points of C^_EE at Delta_g = 100 -- the rows the Delta_g = 300/1000 sweep (w1, other agent)
does not produce: the Delta = 100 row (even J = 0 ... 10, odd J = 5, 7, 9), the even J = 10 Delta = 200 row, and the odd extra rows
m = 1.05/1.1/1.2/1.3 Lam(100) (p8_k2 keys, non-even-integer Delta).  State + heavy reference (Delta_ref = 2e4, same nu), x nodes = those
of the existing C3H (o_H q) row, so that EE rows replace the q^2 C2 rows one for one.
Integrand: w0_ee_common (E = V042^2 (H24 - 2 V204 V402)/(X20 X40) in place of q; imported read-only).  Project modules read-only; only
cache names are redirected:  cache/y1_pts_EE_J<J>_k2.pkl (even, l9_points format), cache/y1_p3_EE_J<J>.pkl (odd, p3 format);
new shadow matrices (if any) -> cache/y1_shadow_k2.pkl / cache/y1_p2shadow.pkl (existing shadow caches are only read).
Usage: OMP_NUM_THREADS=1 python3 y1_rows.py even J-list D-list     |   odd J-list D-list|xrows
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
LAM100 = 98.49873095629202


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def delta_of(mfl, J):
    M = mfl*LAM100; b = 2*J + 1
    return round((-b + np.sqrt(b*b + 4*M*M))/2 + J + 2, 4)


def _pk(pats, odd=False):
    R = {}
    for p in pats:
        for fn in sorted(glob.glob(os.path.join(C, p))):
            d = pickle.load(open(fn, 'rb'))
            R.update({(k[0], k[1]): v for k, v in d.items() if (not odd or k[2] == 1.0)} if odd else d)
    return R


def have(R, D, nu):
    return any(k[0] == float(D) and abs(k[1] - nu) < 1e-7*nu for k in R)


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(C, 'y1_pts_%s_J%d_k2.pkl' % (kind, J))
    sfn = os.path.join(C, 'y1_shadow_k2.pkl')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) + [sfn]):
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


def setup_odd():
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    P._fn = lambda kind, J: os.path.join(C, 'y1_p3_%s_J%d.pkl' % (kind, J))
    sfn = os.path.join(C, 'y1_p2shadow.pkl')
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + [sfn]):
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


def even_nodes(J, D):
    """x nodes of the existing C3H (o_H q) row at Delta (l9 / q2 caches) that also have a C3H reference."""
    Rk = _pk(['l9_pts_C3H_J%d_k2.pkl' % J, 'l9_pts_C3H_J%d.pkl' % J, 'q2_pts_C3H_J%d_k2*.pkl' % J])
    m = mphys(D, J)
    xs = sorted(set(round(k[1]/m, 6) for k in Rk if k[0] == float(D)))
    return [x for x in xs if have(Rk, DREF, x*m)]


def odd_nodes(J, D):
    Rk = _pk(['p3_pts_C3H_J%d.pkl' % J, 'p8_k2_C3H_J%d.pkl' % J, 'q2_p3_C3H_J%d_*.pkl' % J], odd=True)
    m = mphys(D, J)
    xs = sorted(set(round(k[1]/m, 6) for k in Rk if k[0] == float(D)))
    return [x for x in xs if have(Rk, DREF, x*m)]


def xrow_deltas(J):
    Rk = _pk(['p8_k2_C3H_J%d.pkl' % J], odd=True)
    return sorted(set(k[0] for k in Rk if 100.0 < k[0] < 150.0))


if __name__ == '__main__':
    par = sys.argv[1]
    Js = [int(a) for a in sys.argv[2].split(',')]
    if par == 'evenx':                 # even extra rows m = mfl Lam(100) (non-even-integer Delta), kinds EE, C3H, C2 (q^2), fixed x list
        P = setup_even()
        XSX = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
        for J in Js:
            for mfl in [float(v) for v in sys.argv[3].split(',')]:
                D = delta_of(mfl, J); m = mphys(D, J); t0 = time.time(); n = 0
                for kind in ('EE', 'C3H', 'C2'):
                    R = _pk(['y1_pts_%s_J%d_k2.pkl' % (kind, J)])
                    for x in XSX:
                        for DD in (D, DREF):
                            if not have(R, DD, x*m):
                                P.get(kind, J, DD, [x*m], verbose=False); n += 1
                print(f"done evenx J={J} mfl={mfl} Delta={D:.4f} m={m:.2f}: {n} new points in {time.time() - t0:.0f}s", flush=True)
        sys.exit(0)
    P = setup_even() if par == 'even' else setup_odd()
    for J in Js:
        Ds = xrow_deltas(J) if sys.argv[3] == 'xrows' else [float(a) for a in sys.argv[3].split(',')]
        for D in Ds:
            t0 = time.time(); n = 0
            m = mphys(D, J)
            xs = even_nodes(J, D) if par == 'even' else odd_nodes(J, D)
            R = (_pk(['y1_pts_EE_J%d_k2.pkl' % J, 'w1_pts_EE_J%d_k2*.pkl' % J, 'u5_pts_EE_J%d_k2*.pkl' % J]) if par == 'even' else
                 _pk(['y1_p3_EE_J%d.pkl' % J, 'w1_p3_EE_J%d*.pkl' % J], odd=True))      # skip points the 300/1000 sweep already has
            for x in xs:
                for DD in (D, DREF):
                    if par == 'even':
                        R.update(_pk(['w1_pts_EE_J%d_k2*.pkl' % J]))
                    if not have(R, DD, x*m):
                        P.get('EE', J, DD, [x*m], verbose=False); n += 1
            print(f"done {par} EE J={J} Delta={D:.6g} m={m:.2f} x={xs}: {n} new points in {time.time() - t0:.0f}s", flush=True)
