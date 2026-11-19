"""
w1 (EE_NEARGAP.md): production driver for K2 (L1_DEG = 6) crystal points, state + heavy reference (Delta_ref = 2e4, same nu):
  even J: kinds EE (w0 integrand), C2 (q^2), C3H (o_H q)  -- l9_points engine;  'refb' = six q^2 AND EE points at the reference
          itself (x_ref = 0.3 ... 0.8; needed by l8_analysis.basis for the reference helicity basis)
  odd  J: kinds EE, C2, C3H -- p3_odd_action engine (1x1 parity-odd +-; n0 = +1)
Project modules imported read-only; only cache names are redirected:
  cache/w1_pts_<kind>_J<J>_k2<tag>.pkl (even), cache/w1_p3_<kind>_J<J><tag>.pkl (odd),
  shadows: looked up in every existing K2 shadow cache (l9/q2/u5/r1/w1 resp. p2/q2_p2), new ones -> cache/w1_shadow_k2<tag>.pkl /
  cache/w1_p2shadow<tag>.pkl.
Usage: OMP_NUM_THREADS=1 W1_TAG=a python3 w1_rows.py even|odd KINDS J-list D-list [x-list]
       OMP_NUM_THREADS=1 W1_TAG=a python3 w1_rows.py refb J-list
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

TAG = os.environ.get('W1_TAG', '')
C = os.path.join(HERE, 'cache')
DREF = 20000.0
XS19 = [0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8]


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(C, 'w1_pts_%s_J%d_k2%s.pkl' % (kind, J, TAG))
    sfn = os.path.join(C, 'w1_shadow_k2%s.pkl' % TAG)

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl')))):
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
    P._fn = lambda kind, J: os.path.join(C, 'w1_p3_%s_J%d%s.pkl' % (kind, J, TAG))
    sfn = os.path.join(C, 'w1_p2shadow%s.pkl' % TAG)
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl')))):
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


def existing_even(kind, J):
    """all existing records of (kind, J) in the project caches (to skip recomputation): q2, l9, u5 (EE), w1."""
    R = {}
    pats = {'EE': ['u5_pts_EE_J%d_k2*.pkl', 'w1_pts_EE_J%d_k2*.pkl'],
            'C2': ['l9_pts_C2_J%d_k2.pkl', 'q2_pts_C2_J%d_k2*.pkl', 'w1_pts_C2_J%d_k2*.pkl'],
            'C3H': ['l9_pts_C3H_J%d_k2.pkl', 'q2_pts_C3H_J%d_k2*.pkl', 'w1_pts_C3H_J%d_k2*.pkl']}[kind]
    for p in pats:
        for fn in sorted(glob.glob(os.path.join(C, p % J))):
            R.update(pickle.load(open(fn, 'rb')))
    return R


def have(R, D, nu):
    if (float(D), round(float(nu), 9)) in R:
        return True
    return any(k[0] == float(D) and abs(k[1] - nu) < 1e-7*nu for k in R)


if __name__ == '__main__':
    par = sys.argv[1]
    if par == 'refb':
        P = setup_even()
        for J in [int(a) for a in sys.argv[2].split(',')]:
            t0 = time.time()
            mr = mphys(DREF, J)
            for kind in ('C2', 'EE'):
                R = existing_even(kind, J)
                nus = [x*mr for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8) if not have(R, DREF, x*mr)]
                if nus:
                    P.get(kind, J, DREF, nus, verbose=False)
            print(f"done refb J={J}: {time.time() - t0:.0f}s", flush=True)
        sys.exit(0)
    kinds = sys.argv[2].split(',')
    Js = [int(a) for a in sys.argv[3].split(',')]
    Ds = [float(a) for a in sys.argv[4].split(',')]
    xarg = sys.argv[5] if len(sys.argv) > 5 else ''
    xs0 = [float(a) for a in xarg.split(',')] if xarg not in ('', 'keys') else XS19
    P = setup_even() if par == 'even' else setup_odd()
    for J in Js:
        for D in Ds:
            t0 = time.time(); n = 0
            m = mphys(D, J)
            xs = xs0
            if xarg == 'keys':                      # the x nodes of the existing C3H row (so that EE rows match q^2 keys)
                if par == 'even':
                    Rk = existing_even('C3H', J)
                else:
                    Rk = {}
                    for p in ('p3_pts_C3H_J%d.pkl', 'q2_p3_C3H_J%d_*.pkl'):
                        for fn in sorted(glob.glob(os.path.join(C, p % J))):
                            Rk.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0})
                xs = sorted(set(round(k[1]/m, 6) for k in Rk if k[0] == float(D)))
                xs = [x for x in xs if have(Rk, DREF, x*m)]
            for kind in kinds:
                if par == 'even':
                    R = existing_even(kind, J)
                else:
                    R = {}
                    for p in ('p3_pts_%s_J%d.pkl', 'q2_p3_%s_J%d_*.pkl', 'w1_p3_%s_J%d*.pkl'):
                        for fn in sorted(glob.glob(os.path.join(C, p % (kind, J)))):
                            R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0})
                for x in xs:
                    for DD in (D, DREF):
                        if not have(R, DD, x*m):
                            P.get(kind, J, DD, [x*m], verbose=False); n += 1
            print(f"done {par} {','.join(kinds)} J={J} Delta={D:.4g} m={m:.1f} J^2/m={J*J/m:.2f}: {n} new points in "
                  f"{time.time() - t0:.0f}s", flush=True)
