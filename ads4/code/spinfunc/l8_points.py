"""
l8_points (M6 round 4): cached evaluation of the graviton C^_2 crystal action (l6_spin_action.action) at arbitrary (J, Delta, nu), with the
Laplace orders kept SEPARATE (needed for the cumulant / log form, doc 22.2).  One pickle per J: cache/l8_pts_J<J>.pkl, keyed (Delta, nu)
(nu rounded to 1e-9).  Stored: M0, M2, M4 (order eps^0, eps^2, eps^4 contributions of the 2x2 / 1x1 action, structure basis of
l5_spin.conserved(Delta, J)), logC, logK, scalar (logY0, S0, S2, S4), m.
Usage: python3 l8_points.py J Delta nu1 [nu2 ...]       (or import and call get(J, Delta, nus))
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l6_spin_action as A
import l1_crystal as Cr


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def _fn(J):
    return os.path.join(HERE, 'cache', 'l8_pts_J%d.pkl' % J)


def load(J):
    fn = _fn(J)
    return pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}


def get(J, Delta, nus, verbose=True):
    res = load(J)
    out = []
    for nu in nus:
        key = (float(Delta), round(float(nu), 9))
        if key not in res:
            t = time.time()
            r = A.action(float(nu), float(Delta), J, with_scalar=False)
            s = Cr.laplace_pref(float(nu), float(Delta), J)
            M = r['M']
            rec = dict(M0=M[:, :, 0], M2=M[:, :, 1] - M[:, :, 0], M4=M[:, :, 2] - M[:, :, 1], logC=r['logC'], logK=r['logK'],
                       slog=s['logY0'] - 3.5*np.log(2*np.pi), sS=np.array(s['S'][:3]), m=mphys(Delta, J), shadow_resid=r['shadow_resid'])
            res = load(J); res[key] = rec                                   # re-load: another process may have written
            pickle.dump(res, open(_fn(J) + '.tmp%d' % os.getpid(), 'wb')); os.replace(_fn(J) + '.tmp%d' % os.getpid(), _fn(J))
            if verbose:
                print(f"J={J} Delta={Delta} nu={nu:.4f}: {time.time() - t:.1f}s", flush=True)
        out.append(res[key])
    return out


if __name__ == '__main__':
    J = int(sys.argv[1]); Dl = float(sys.argv[2])
    get(J, Dl, [float(v) for v in sys.argv[3:]])
