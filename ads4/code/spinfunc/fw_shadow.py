"""
fw_shadow (FOLD_WINDOW_LARGEJ.md §1): robust spinning shadow matrix b(Delta, J) at large J.
l5_shadow_matrix.shadow_matrix fits S_c(t) = sum_b b_cb r_b(t) over 6 random complex null polarisation triples t by an UNSCALED lstsq.
At J >~ 56 the triples' magnitudes differ by many orders (V^J factors), one triple dominates and b comes out near-singular and seed-
dependent: (1064, 64) seed 7: det b = 8e-18; seed 11: det 0.57 (the cached ee1000 b for J = 56, 64 and r1's J = 72 are of the first kind).
Fix: row equilibration (each triple scaled by its rhs norm) and column equilibration, ntrip = 12; checked seed/xi-independent to ~1e-4
(J = 64) [NE].  Same Laplace order (K2, L1_DEG = 6) and logK as l5.  Output: cache/fw_shadowk2.pkl {(Delta, J): (b, logK, resid)}.
Usage: OMP_NUM_THREADS=1 python3 fw_shadow.py D1:J1,D2:J2,...  [--check]
"""
import os, sys, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l5_shadow_matrix as SM
import fw_frame as FW
FW.patch_conserved_clean()          # J >= 150

FN = os.path.join(HERE, 'cache', 'fw_shadowk2.pkl')


def robust(D, J, xi=0.5, seed=7, ntrip=12):
    rng = np.random.default_rng(seed)
    trips = [(SM._rnull(rng), SM._rnull(rng), SM._rnull(rng)) for _ in range(ntrip)]
    logI0, Sv = SM.lhs_series(D, J, xi, trips)
    den, r = SM.rhs(D, J, xi, trips)
    rs = np.abs(r).max(axis=1, keepdims=True)
    r2 = r/rs; S2 = Sv[:, :, 2]/rs
    cs = np.abs(r2).max(axis=0)
    b, *_ = np.linalg.lstsq(r2/cs, S2, rcond=None)
    b = b/cs[:, None]
    res = np.linalg.norm(r2 @ b - S2)/np.linalg.norm(S2)
    return b.T, logI0 - den, res


def get(D, J):
    db = pickle.load(open(FN, 'rb')) if os.path.exists(FN) else {}
    key = (float(D), int(J))
    if key not in db:
        b, lK, res = robust(D, J)
        db = pickle.load(open(FN, 'rb')) if os.path.exists(FN) else {}
        db[key] = (b, lK, [res, res, res])
        pickle.dump(db, open(FN + '.tmp%d' % os.getpid(), 'wb')); os.replace(FN + '.tmp%d' % os.getpid(), FN)
    return db[key]


if __name__ == '__main__':
    for spec in sys.argv[1].split(','):
        D, J = float(spec.split(':')[0]), int(spec.split(':')[1])
        t = time.time()
        if '--check' in sys.argv:
            b0, l0, r0 = robust(D, J)
            b1, l1, r1 = robust(D, J, 0.4, 11)
            print(f"check Delta={D:g} J={J}: seed/xi dev {np.abs(np.exp(l1 - l0)*b1 - b0).max()/np.abs(b0).max():.1e}, resid {r0:.1e}/{r1:.1e}, "
                  f"det {abs(np.linalg.det(b0)):.2e}", flush=True)
            continue
        b, lK, res = get(D, J)
        print(f"Delta={D:g} J={J}: b = {np.round(b, 5).ravel()}, det {abs(np.linalg.det(b)):.2e}, resid {res[0]:.1e}, {time.time() - t:.0f}s", flush=True)
