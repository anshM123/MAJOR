"""
ee1000_direct (EE_DG1000_LARGEJ.md): direct (RQMC, n1 contour; n1_run.point unchanged) crystal values of C^_EE (w3_direct.setup_pw) at
Delta_g = 1000 near-gap rows, to validate the K2 rows at large J^2/m and large x.  State and reference (Delta_ref = 2e4, same nu).
Cache: cache/ee1000_n1_<kind>_J<J>.pkl keyed (Delta, nu rounded 1e-9).  Resumable.
Usage: OMP_NUM_THREADS=1 python3 ee1000_direct.py J Delta x1,x2,... [EE|C3H]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
C = os.path.join(HERE, 'cache')


def dfn(J, kind='EE'):
    return os.path.join(C, 'ee1000_n1_%s_J%d.pkl' % (kind, J))


def dload(J, kind='EE'):
    f = dfn(J, kind)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


if __name__ == '__main__':
    import w3_direct as W3
    W3.setup_pw()
    import n1_run as R
    J = int(sys.argv[1]); Delta = float(sys.argv[2]); xs = [float(v) for v in sys.argv[3].split(',')]
    kind = sys.argv[4] if len(sys.argv) > 4 else 'EE'
    m = R.mphys(Delta, J)
    for x in xs:
        for D in (Delta, 20000.0):
            nu = x*m
            key = (D, round(nu, 9))
            if key in dload(J, kind):
                continue
            t = time.time()
            r = R.point(kind, J, D, nu)
            db = dload(J, kind); db[key] = r
            pickle.dump(db, open(dfn(J, kind) + '.tmp', 'wb')); os.replace(dfn(J, kind) + '.tmp', dfn(J, kind))
            rel = np.max(np.abs(r['err']))/np.max(np.abs(r['E']))
            print(f"direct {kind} J={J} Delta={D:.0f} x={x} nu={nu:.3f}: max rel err {rel:.1e}, nmask {r['nmask']}, {time.time() - t:.0f}s",
                  flush=True)
