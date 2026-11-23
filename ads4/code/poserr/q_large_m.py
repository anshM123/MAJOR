"""
poserr/q_large_m.py (POS_ERRORS.md §2.4): per-state budget of Result Q's upper functional at large mass, flat parts only (exact Wigner
functions by recurrence; Q3 weights from prl/data/q3_full_delta0.2.json, read-only; the (1-p)^3 rounding defect of CERT_CLOSE §2.1 is
irrelevant at m >= 1.3 Lam, u <= 0.6).  Per (m/Lam, J, h), Stage-0 normal form (u = Lam^2/m^2, y = u p^2):
   FF = B2 + B3 + FW - 2 c u^2,  B2 = int psi2 (2-y) D_h, B3 = s u int psi3 D_h, FW = a (2J(J+1) - c4 + s) u^3,
   N  = |B2| + |B3| + |FW| + 2 c u^2
and the sign-free costs (in FF/N units):
   Wq  = sum (|psi2| env2 + u |psi3| env3)/(m^2 N)       per uniform |Gamma| = 1 c2 on the smeared k = 2, 3 parts (split nu < 12 / >= 12)
   Wd  = 2 a u^3 / N                                     per unit s_d (derivative-piece correction, enters FF as a 2 u^3 s_d)
Output cache/q_large_m_<Dg>.json.  Usage: python3 q_large_m.py Dg [mfl,...] [bmax]
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
from large_m_map import jgrid, jmax
CA = os.path.join(HERE, 'cache')
QJ = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'q3_full_delta0.2.json')))
PSI2 = np.array([float(v) for v in QJ['psi2']]); PSI3 = np.array([float(v) for v in QJ['psi3']])
A = float(QJ['a']); CQ = float(QJ['c'])
SIG = {'++': 1.0, '+-': -1.0}; C4 = {'++': 1.0, '+-': 41.0}
NUSPLIT = 12.0


def lam_of(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def jac(n, a, b, x):
    p0 = np.ones_like(x)
    if n == 0:
        return p0
    p1 = 0.5*(a - b + (a + b + 2)*x)
    for k in range(2, n + 1):
        c1 = 2*k*(k + a + b)*(2*k + a + b - 2)
        c2 = (2*k + a + b - 1)*((2*k + a + b)*(2*k + a + b - 2)*x + a*a - b*b)
        c3 = 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
        p0, p1 = p1, (c2*p1 - c3*p0)/c1
    return p1


def gl_panels(n):
    k = max(1, int(np.ceil(n/64)))
    x, w = np.polynomial.legendre.leggauss(64)
    e = np.linspace(0, 1, k + 1)
    P = np.concatenate([(e[i] + e[i + 1])/2 + (e[i + 1] - e[i])/2*x for i in range(k)])
    W = np.concatenate([(e[i + 1] - e[i])/2*w for i in range(k)])
    return P, W


def state(Dg, J, mfl):
    Lm = lam_of(Dg); m = mfl*Lm; u = 1/mfl**2
    ph = (2*J + 1)*np.arcsin(min(1.0, 1/mfl))
    p, w = gl_panels(int(1.6*ph + 256))
    y = u*p*p; x = np.sqrt(y); nu = p*Lm
    ps2 = np.polynomial.polynomial.polyval(p, PSI2); ps3 = np.polynomial.polynomial.polyval(p, PSI3)
    out = dict(J=J, mfl=mfl, bhat=(2*J + 1)/mfl)
    for h in (['++'] if J < 4 else ['++', '+-']):
        if h == '++' and J % 2:
            continue
        s = SIG[h]
        D = jac(J, 0, 0, 1 - 2*y) if h == '++' else jac(J - 4, 0, 8, 1 - 2*y)
        B2 = (w*ps2*(2 - y)*D).sum(); B3 = s*u*(w*ps3*D).sum()
        FW = A*(2*J*(J + 1) - C4[h] + s)*u**3; V = -2*CQ*u*u
        N = abs(B2) + abs(B3) + abs(FW) + 2*CQ*u*u
        e2 = (2 - y) if h == '++' else (2 - y)/(1 - y)**4
        e3 = np.ones_like(y) if h == '++' else 1/(1 - y)**4
        wv = w*(np.abs(ps2)*e2 + u*np.abs(ps3)*e3)/m**2/N
        out[h] = dict(FF=float((B2 + B3 + FW + V)/N), N=float(N), B2=float(B2), B3=float(B3), FW=float(FW),
                      Wq_lo=float(wv[nu < NUSPLIT].sum()), Wq_hi=float(wv[nu >= NUSPLIT].sum()), Wd=float(2*A*u**3/N))
    return out


if __name__ == '__main__':
    Dg = float(sys.argv[1])
    mfls = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.3, 1.5, 2, 3, 5, 10, 20]
    bmax = float(sys.argv[3]) if len(sys.argv) > 3 else 400.0
    res = dict(Dg=Dg, states=[]); t0 = time.time()
    for mfl in mfls:
        for J in jgrid(Dg, mfl, bmax):
            res['states'].append(state(Dg, J, mfl))
        print('Dg %g mfl %g: %.0fs' % (Dg, mfl, time.time() - t0), flush=True)
    json.dump(res, open(os.path.join(CA, 'q_large_m_%g.json' % Dg), 'w'))
