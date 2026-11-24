"""
edger2/e2_k2.py (EDGE_R2.md §2): K2 (L1_DEG = 6) crystal records and <TTO> shadow matrices for the even-J helicity bases, with
spinfunc/w1_rows.py's even setup (l9_points engine, w0 EE integrand) imported read-only; ONLY the cache paths are redirected to
edger2/cache/k2/ (w1_pts_<kind>_J<J>_k2_e2.pkl, w1_shadow_k2_e2.pkl, w1_p2shadow_e2.pkl: the chain's own file-name patterns, so
that e2_gamma can point the chain's loaders at this directory); existing project shadow caches are read first.
  refb  J-list         : the six EE (and C2) reference points x_ref = 0.3 ... 0.8 at Delta_ref = 2e4 (l8_analysis.basis of the reference)
  state J D x-list     : EE K2 points of the state row (helicity basis of the state, x >= 0.15)
  shadow J D           : the even <TTO> shadow b(D, J) (l5_shadow_matrix, xi = 0.5, 6 triples, seed 7)
  oshadow J D          : the odd shadow (p2_shadow_odd.shadow_all, 8 triples, seed 7), key (D, J, 6)
Usage: python3 e2_k2.py refb 34,38  |  python3 e2_k2.py state 34 1034 0.2,0.3,...  |  python3 e2_k2.py shadow 34 1034
"""
import e2_guard  # noqa: F401  (first)
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF)
CACHE = os.path.join(HERE, 'cache', 'k2')
os.makedirs(CACHE, exist_ok=True)
DREF = 20000.0
import numpy as np


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def setup():
    cwd = os.getcwd(); os.chdir(SF)
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    os.chdir(cwd)
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(CACHE, 'w1_pts_%s_J%d_k2_e2.pkl' % (kind, J))
    sfn = os.path.join(CACHE, 'w1_shadow_k2_e2.pkl')
    C = os.path.join(SF, 'cache')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) + [sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        t = time.time()
        r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp%d' % os.getpid(), 'wb')); os.replace(sfn + '.tmp%d' % os.getpid(), sfn)
        print(f"shadow {key}: {time.time() - t:.1f}s, resid {np.max(np.abs(r['resid'])):.2e}", flush=True)
        return db[key]
    P.shadow = shadow
    P.load = lambda kind, J: pickle.load(open(P._fn(kind, J), 'rb')) if os.path.exists(P._fn(kind, J)) else {}
    return P


if __name__ == '__main__':
    P = setup()
    if sys.argv[1] == 'refb':
        for J in [int(a) for a in sys.argv[2].split(',')]:
            t0 = time.time(); mr = mphys(DREF, J)
            kinds = sys.argv[3].split(',') if len(sys.argv) > 3 else ['EE']
            for kind in kinds:
                P.get(kind, J, DREF, [x*mr for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8)])
            print(f"done refb J={J}: {time.time() - t0:.0f}s", flush=True)
    elif sys.argv[1] == 'shadow':
        P.shadow(float(sys.argv[3]), int(sys.argv[2]))
    elif sys.argv[1] == 'oshadow':
        cwd = os.getcwd(); os.chdir(SF)
        import p2_shadow_odd as P2
        os.chdir(cwd)
        J, D = int(sys.argv[2]), float(sys.argv[3])
        sfn = os.path.join(CACHE, 'w1_p2shadow_e2.pkl')
        key = (D, J, P2.T.DEG)
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        if key not in db:
            t = time.time(); r = P2.shadow_all(D, J)
            db[key] = (r['b'], r['logK'], r['resid'])
            pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
            print(f"odd shadow {key}: {time.time() - t:.1f}s resid {r['resid']}", flush=True)
    else:
        J, D = int(sys.argv[2]), float(sys.argv[3]); m = mphys(D, J)
        P.get('EE', J, D, [float(x)*m for x in sys.argv[4].split(',')])
