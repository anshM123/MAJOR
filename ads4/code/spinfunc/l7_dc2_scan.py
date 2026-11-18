"""
l7 (M6 step 2, items 3-5 data): scan of the graviton C^_2 crystal action (l6) against the scalar k = 2 action (l1) on a grid of
physical x = nu/m (m^2 = alpha(alpha + 2J + 1), alpha = Delta - J - 2), J (even, parity-even sector) and Delta.  Resumable: one pickle per J,
keyed (Delta, x).  Analysis in l7_dc2_analysis.py.
Usage: python3 l7_dc2_scan.py J1 [J2 ...]   (env L7_DELTAS, L7_XS to override the grids)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l6_spin_action as A

DELTAS = [float(v) for v in os.environ.get('L7_DELTAS', '100,200,400,800,1600,3200').split(',')]
XS = [float(v) for v in os.environ.get('L7_XS', '0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def run(J):
    fn = os.path.join(HERE, 'cache', 'l7_J%d.pkl' % J)
    res = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    for Dl in DELTAS:
        for x in XS:
            if (Dl, x) in res:
                continue
            t = time.time()
            m = mphys(Dl, J)
            r = A.action(x*m, Dl, J)
            res[(Dl, x)] = dict(M=r['M'], logC=r['logC'], logK=r['logK'], sl=r['scal_log'], sS=r['scal_S'], m=m,
                                shadow_resid=r['shadow_resid'])
            pickle.dump(res, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
            print(f"J={J} Delta={Dl} x={x}: {time.time() - t:.1f}s", flush=True)


if __name__ == '__main__':
    for J in [int(a) for a in sys.argv[1:]]:
        run(J)
