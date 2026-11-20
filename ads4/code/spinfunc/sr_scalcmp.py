"""
sr_scalcmp (SUBLEADING_REGGE.md): exact scalar-skeleton heavy action (Mellin, exactskel/skelexact.ExactState, Df = 3) at the nu-nodes
of an EE K2 row, in the same c2 units as the row's Gamma:  Gs(x) = m^2 [ex_2/E_2 - P_J(1-2x^2)],  E_2 = (2-x^2)/(1-x^2)^2,
and the leading-Regge curvature part m^2 2(P_l(cosh b) - J0((2J+1)x))/E_2.  Read-only imports.  Caches cache/sr_scal_<D>_<J>.pkl.
Usage: python3 sr_scalcmp.py Delta J
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'exactskel'))
import numpy as np


def exact_scalar(Delta, J, nus, Df=3.0):
    fn = os.path.join(HERE, 'cache', 'sr_scal_%g_%d.pkl' % (Delta, J))
    db = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
    need = [nu for nu in nus if round(float(nu), 6) not in db]
    if need:
        from skelexact import ExactState
        es = ExactState(Delta, J, Df)
        v = es.values(np.array(need), ks=(2, 3))
        for i, nu in enumerate(need):
            db[round(float(nu), 6)] = (v[2][i], v[3][i])
        pickle.dump(db, open(fn + '.tmp', 'wb')); os.replace(fn + '.tmp', fn)
    return np.array([db[round(float(nu), 6)] for nu in nus])


if __name__ == '__main__':
    D, J = float(sys.argv[1]), int(sys.argv[2])
    al = D - J - 2; m = np.sqrt(al*(al + 2*J + 1))
    xs = np.array([float(a) for a in sys.argv[3].split(',')]) if len(sys.argv) > 3 else np.arange(0.025, 0.81, 0.025)
    t0 = time.time()
    v = exact_scalar(D, J, xs*m)
    print('time', time.time() - t0)
