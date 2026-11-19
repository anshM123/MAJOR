"""
p8_k2rows (PARITY_ODD_ROUTEA.md sect. 11): the Laplace (p3_odd_action, DEG 6 = K2; L1_DEG=8 -> K3) odd-J points on the extra near-gap
rows m(Delta, J) = mfl Lam(100) (p8_odd_pw.delta_of), state + reference (Delta_ref = 2e4, same nu), written to SEPARATE caches
cache/p8_k2_<kind>_J<J>.pkl (cache/p8_k3_... at L1_DEG=8) so that the p3 caches used by p4/p6 (sect. 6) are untouched.
These supply x < 0.2 (and the a_{+-} fit points) of the extra rows; the x >= 0.2 values are compared with the direct ones.
Usage: OMP_NUM_THREADS=1 P8_MFL=1.05,1.1,1.2,1.3 python3 p8_k2rows.py kind J [J ...]   (P8K_XS: x list, default 0.05 ... 0.8)
"""
import os, sys
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import p3_odd_action as P

LAM100 = 98.49873095629202
MFL = [float(v) for v in os.environ.get('P8_MFL', '1.05,1.1,1.2,1.3').split(',') if v]
XS = [float(v) for v in os.environ.get('P8K_XS', '0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8').split(',')]
DREF = 20000.0
TAG = 'k3' if P.T.DEG >= 8 else 'k2'


def delta_of(mfl, J):
    import numpy as np
    M = mfl*LAM100; b = 2*J + 1
    return round((-b + np.sqrt(b*b + 4*M*M))/2 + J + 2, 4)


P._fn = lambda kind, J: os.path.join(HERE, 'cache', 'p8_%s_%s_J%d.pkl' % (TAG, kind, J))

if __name__ == '__main__':
    kind = sys.argv[1]
    for J in [int(a) for a in sys.argv[2:]]:
        for r in MFL:
            D = delta_of(r, J)
            m = P.mphys(D, J)
            for x in XS:
                P.get(kind, J, D, [x*m])
                P.get(kind, J, DREF, [x*m])
        print(f"done {kind} J={J}", flush=True)
