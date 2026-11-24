"""
b5_k2 (BATCH5_FOLD.md §1): the K2 (L1_DEG = 6) EE STATE records that the even chain's l8 helicity basis needs on an own row
Delta = 1000 + J that has no K2 row (even J = 18 ... 32): l8_analysis.basis(J, D, R) diagonalises the order-0 matrices M0 of >= 4
K2 records with x >= 0.15 at that Delta; direct records carry no 'slog'/'sS' and are not used for it (w2_chain._basis reads the K2
loaders only).  This is ngfill/nf_rows.py's even set-up (= neargap2/ng_rows = spinfunc/w1_rows, l9_points engine with the w0 EE
integrand, imported read-only) with the cache moved to b5fold/cache (w1_pts_EE_J<J>_k2_b5.pkl) and the shadow read from
cache/b5_shadow_k2.pkl (b5_shadow.py).  STATE points only (no Delta_ref = 2e4 record at these nu): they are NOT row nodes (a node needs
state and reference of both kinds); every node of the row stays a direct n1 node (EDGE_R2's e2_k2 'state' mode does the same).
Usage: python3 b5_k2.py J x1,x2,...   -> cache/w1_pts_EE_J<J>_k2_b5.pkl
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b5_guard  # noqa: F401  (first)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF)
os.chdir(SF)
import w1_rows as W1
B5C = os.path.join(HERE, 'cache')
SFC = os.path.join(SF, 'cache')


def setup_even():
    import l9_points as P
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(B5C, 'w1_pts_%s_J%d_k2_b5.pkl' % (kind, J))
    sfn = os.path.join(B5C, 'b5_shadow_k2.pkl')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(SFC, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(SFC, 'u5_shadow_k2.pkl'), os.path.join(SFC, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(SFC, 'w1_shadow_k2*.pkl'))) + [sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        raise KeyError(('no shadow (run b5_shadow.py first)', key))
    P.shadow = shadow
    return P


if __name__ == '__main__':
    J = int(sys.argv[1]); D = 1000.0 + J
    xs = [float(v) for v in sys.argv[2].split(',')]
    P = setup_even()
    m = W1.mphys(D, J)
    t0 = time.time()
    for x in xs:
        fn = P._fn('EE', J)
        R = {(k[0], k[1]): v for k, v in (pickle.load(open(fn, 'rb')).items() if os.path.exists(fn) else [])}
        if not W1.have(R, D, x*m):
            t1 = time.time(); P.get('EE', J, D, [x*m], verbose=False)
            print(f"  EE J={J} D={D:g} x={x:.4f} nu={x*m:.3f}: {time.time() - t1:.1f}s", flush=True)
    print(f"done EE state K2 points J={J} Delta={D:g} m={m:.3f}: {time.time() - t0:.0f}s", flush=True)
