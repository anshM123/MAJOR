"""
fullpos/fp_common.py (FULLPOS_CERT.md): exact data for the angular finite-mass action of Omega^EE.

Weight (read-only): prl/data/k2/k2_continuum_D_<Dg>_rep.json (endpoint-repaired; EPSW_CERT §3), all doubles taken exactly.
  chi_re(nu) = Q(nu) sum_k x_k [G_k(nu-c_k) + G_k(nu+c_k)] + sum_j y_j nu^{2j} (nu^2+1) log(nu^2+1) e^{-nu^2/2},  Q = nu^2/(nu^2+1/4)^2
  theta(nu)  = (1/2)[erf((nu+nu_c)/a) - erf((nu-nu_c)/a)],  nu_c = Lam + 100, a = 10   (UNIFIED_THEOREM_v2 §2.1, D's parameters)
  g(nu)      = -9216 (4nu^2+1)/((4nu^2+9)(4nu^2+49)) chi_re(nu) theta(nu) / Lam
CDH polynomials (verify_cdh.py): x = nu^2, s = x + 1/4
  d = 0:  f_0 = 0, f_1 = 1, (m+1)^2 f_{m+1} = -(s - 2m^2) f_m - (m-1)^2 f_{m-1}           (f_m = Q^{(0)}_{m-1})
  d = 1:  q_0 = 1, q_{-1} = 0, (n+3)^2 q_{n+1} = -(s + 4 - 2(n+2)^2) q_n - n^2 q_{n-1}       (q_n = Q^{(1)}_n)
  z_m = 12 f_m - (m-1) f_{m-1} + (m+1) f_{m+1}
  y_m = ((12-m^2)/2) f_m + ((13-7m)/2) f_{m-1} + ((13+7m)/2) f_{m+1} + ((m-2)^2/4) f_{m-2} + ((m+2)^2/4) f_{m+2}
Moments: Y_k = int g y_{k+2}, Z_k = int g z_{k+2}, Q_k = int_0^1 psi3(p) q_k(Lam^2 p^2) dp (exact rational).
"""
import os, sys, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
CACHE = os.path.join(HERE, 'cache')
K0 = Fr('27.48509817')


def wpath(Dg):
    return os.path.join(DATA, 'k2', f'k2_continuum_D_{Dg}_rep.json')


def load_weight(Dg):
    d = json.load(open(wpath(Dg)))
    assert d['Delta_g'] == Dg and 'repair' in d
    cs = d['centres_sigma']
    xs = [Fr(float(s)) for s in d['x_re']]
    ys = [Fr(float(s)) for s in d['y_re']]
    return dict(Dg=Dg, lam=Fr((Dg - 2)*(Dg - 1)), cs=[(Fr(c), Fr(s)) for c, s in cs], x=xs, y=ys)


def nf02():
    nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    p2 = [Fr(c)/K0 for c in nf['psi2_exact']]
    p3 = [Fr(c)/K0 for c in nf['psi3_exact']]
    return p2, p3


# ------------------------------------------------------------------ exact CDH polynomials in x = nu^2 (Fraction coefficient lists)
def padd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]


def pscale(a, c):
    return [c*t for t in a]


def pmulx(a):
    return [Fr(0)] + list(a)


def cdh0(M):
    """f_0..f_M as exact coefficient lists in x."""
    f = [[Fr(0)], [Fr(1)]]
    for m in range(1, M):
        # (m+1)^2 f_{m+1} = -(x + 1/4 - 2m^2) f_m - (m-1)^2 f_{m-1}
        t = padd(pmulx(f[m]), pscale(f[m], Fr(1, 4) - 2*m*m))
        t = padd(pscale(t, -1), pscale(f[m - 1], -(m - 1)**2))
        f.append(pscale(t, Fr(1, (m + 1)**2)))
    return f


def cdh1(N):
    """q_0..q_N (d = 1)."""
    q = [[Fr(1)]]
    prev = [Fr(0)]
    for n in range(0, N):
        t = padd(pmulx(q[n]), pscale(q[n], Fr(1, 4) + 4 - 2*(n + 2)**2))
        t = padd(pscale(t, -1), pscale(q[n - 1] if n >= 1 else prev, -(n*n)))
        q.append(pscale(t, Fr(1, (n + 3)**2)))
    return q


def zy_polys(K):
    """exact z_{k+2}, y_{k+2} coefficient lists, k = 0..K."""
    f = cdh0(K + 5)
    Z, Y = [], []
    for k in range(K + 1):
        m = k + 2
        Z.append(padd(padd(pscale(f[m], 12), pscale(f[m - 1], -(m - 1))), pscale(f[m + 1], m + 1)))
        y = pscale(f[m], Fr(12 - m*m, 2))
        y = padd(y, pscale(f[m - 1], Fr(13 - 7*m, 2)))
        y = padd(y, pscale(f[m + 1], Fr(13 + 7*m, 2)))
        y = padd(y, pscale(f[m - 2], Fr((m - 2)**2, 4)))
        y = padd(y, pscale(f[m + 2], Fr((m + 2)**2, 4)))
        Y.append(y)
    return Z, Y


def Q_moments(K, lam):
    """Q_k = int_0^1 psi3(p) q_k(lam p^2) dp, exact Fractions, k = 0..K."""
    _, p3 = nf02()
    mom = {}                                   # int_0^1 p^{2i} psi3 dp
    def M(i):
        if i not in mom:
            mom[i] = sum(c/(2*i + n + 1) for n, c in enumerate(p3))
        return mom[i]
    q = cdh1(K)
    return [sum(c*lam**i*M(i) for i, c in enumerate(q[k])) for k in range(K + 1)]


# ------------------------------------------------------------------ ladder data
def ladder(J, b, Delta=None, T=None):
    """T, C, u_k, D_k, V_{b,k}, sigma_b for channel b in {0,4}."""
    from math import comb
    C = J*(J + 1)
    if T is None:
        T = (Delta - 1)*(Delta - 2)
    n = J - b
    u = [(b + k)*(b + k + 1) for k in range(n + 1)]
    D = [comb(k + 3, 3) for k in range(n + 1)]
    if b == 0:
        V = [Fr(4*k, k + 3) for k in range(n + 1)]; sig = 1
    else:
        V = [4 + Fr(12, k + 1) for k in range(n + 1)]; sig = -1
    return dict(T=T, C=C, n=n, u=u, D=D, V=V, sig=sig)
