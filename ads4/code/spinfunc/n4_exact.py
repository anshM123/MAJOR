"""
n4_exact (NEARGAP100_SMALLNU.md): EXACT scalar comparison records (Mellin, exactskel.skelexact.ExactState; read-only import) for
ARBITRARY blocks (Delta, J) -- the Delta = 100 row, the extra m-rows and the Delta_ref = 2e4 references of the graviton chain --
in e1's format (J, m, eta, Delta, nus = Lam(100) p on e1's 64 Gauss nodes, x, ex[2], rho[2], logq), so that n3_exact.X(nu, rec)
gives the exact shifted object X(nu) = g C(nu+i) + g~ C(nu-i) at the production nodes.
Usage: python3 n4_exact.py  (fills cache/n4_exact.pkl for the production states)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
EX = os.path.join(HERE, '..', 'exactskel'); OMU = os.path.join(HERE, '..', 'omegau')
import numpy as np
from scipy.special import j0, eval_legendre

PK = os.path.join(HERE, 'cache', 'n4_exact.pkl')
DF = 3.0


def _gl01(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(x + 1), 0.5*w


def record(Delta, J, nus=None):
    os.environ.setdefault('OMEGAU_DF', '3')
    for p in (EX, OMU):
        if p not in sys.path:
            sys.path.insert(0, p)
    from skelexact import ExactState
    import mpmath as mp
    Lam = np.sqrt(98.0*99.0)
    if nus is None:
        pn, _ = _gl01(64)
        nus = Lam*pn
        nus = nus[nus < 40.0]
    al = Delta - J - 2.0; ab = al + 2*J + 1
    m = np.sqrt(al*ab); eta = (al**2 + ab**2)/(2*al*ab)
    es = ExactState(Delta, J, DF)
    ex = es.values(nus, ks=(2,))
    ex2 = np.asarray(ex[2] if isinstance(ex, dict) else ex, complex).real
    x = nus/m
    E2 = (2 - x*x)/(1 - x*x)**2
    P = np.array([float(mp.legendre(-0.5 + 1j*v, eta).real) for v in nus])
    hyb = E2*eval_legendre(J, 1 - 2*x*x) + 2*(P - j0((2*J + 1)*x))
    return dict(J=J, m=m, eta=eta, Delta=float(Delta), nus=nus, x=x, ex={2: ex2}, rho={2: (ex2 - hyb)/E2}, logq=es.logq)


def get(Delta, J):
    db = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
    key = (round(float(Delta), 9), int(J))
    if key not in db:
        t = time.time()
        db = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
        db[key] = record(Delta, J)
        tmp = PK + '.tmp%d' % os.getpid()
        pickle.dump(db, open(tmp, 'wb')); os.replace(tmp, PK)
        print(f'   n4_exact: Delta={Delta:.4f} J={J}: {time.time() - t:.0f}s', flush=True)
    return db[key]


def X(nu, rec):
    import n3_exact as NX
    return NX.X(nu, rec)[0]


if __name__ == '__main__':
    import n3_run as NR
    Lam = np.sqrt(98.0*99.0)
    sts = []
    for J in (4, 5, 6, 7, 8, 9, 10):
        sts += [(100.0, J), (20000.0, J)]
    for J in (5, 7, 9):
        for mfl in (1.05, 1.1, 1.2, 1.3):
            m = mfl*Lam; Bq = 2*J + 1.0; al = (-Bq + np.sqrt(Bq*Bq + 4*m*m))/2
            sts.append((al + J + 2, J))
    for (D_, J) in sts:
        get(D_, J)
