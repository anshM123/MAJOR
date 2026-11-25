"""
fullpos3/fp3_common.py: shared data and kernels for FULLPOS3.
  psi2, psi3: NF02 exact polynomials / k0 (prl/data/nf02_psi.json, read-only), as Fraction coefficient lists in p.
  F_{p,nu}(z) = (1-z)^{-p} P_{-1/2+i nu}((1+z)/(1-z)) satisfies (FULLPOS2 §1.1)
      z(1-z)^2 F'' + (1 - (2+2p)z + (1+2p)z^2) F' + (s - p + p^2 z) F = 0,   s = nu^2 + 1/4.
  taylorF(p, nu, c, K): Taylor coefficients t_j = F^{(j)}(c)/j!, j <= K, from F(c), F'(c) (Arb Legendre, type 3) and the ODE.
"""
import os, sys, json, math
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, acb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
CACHE = os.path.join(HERE, 'cache')
K0 = Fr('27.48509817')
DG = 1000
LAM2 = (DG - 2)*(DG - 1)


def nf02():
    nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    p2 = [Fr(c)/K0 for c in nf['psi2_exact']]
    p3 = [Fr(c)/K0 for c in nf['psi3_exact']]
    return p2, p3


def polyval_fr(coef, x):
    v = 0.0
    for a in reversed(coef):
        v = v*x + float(a)
    return v


def taylorF(p, nu, c, K):
    """t_0..t_K, F^{(j)}(c)/j! for F_{p,nu}; nu real or complex (acb/arb/float); c arb in (0,1)."""
    nu = acb(nu); c = arb(c)
    q = 1 - c
    X = acb((1 + c)/q)
    mu = acb(-0.5) + acb(0, 1)*nu
    P0 = X.legendre_p(mu, 0, 3)
    Pm = X.legendre_p(mu - 1, 0, 3)
    dP = mu*(X*P0 - Pm)/(X*X - 1)
    qp = q**(-p)
    t = [qp*P0, p*qp/q*P0 + qp*2/(q*q)*dP]
    s = nu*nu + arb(1)/4
    A = [c*q*q, q*q - 2*c*q, c - 2*q, arb(1)]
    B = [1 - (2 + 2*p)*c + (1 + 2*p)*c*c, -(2 + 2*p) + 2*(1 + 2*p)*c, arb(1 + 2*p)]
    Cc = [s - p + p*p*c, arb(p*p)]
    for k in range(0, K - 1):
        acc = acb(0)
        for i in range(1, 4):
            if k - i + 2 >= 0:
                acc += A[i]*(k - i + 2)*(k - i + 1)*t[k - i + 2]
        for i in range(0, 3):
            if k - i + 1 >= 0:
                acc += B[i]*(k - i + 1)*t[k - i + 1]
        for i in range(0, 2):
            if k - i >= 0:
                acc += Cc[i]*t[k - i]
        t.append(-acc/(A[0]*(k + 2)*(k + 1)))
    return t[:K + 1]


def stirling2(m, j, _c={}):
    if (m, j) in _c:
        return _c[(m, j)]
    if m == j:
        v = 1
    elif j == 0 or j > m:
        v = 0
    else:
        v = j*stirling2(m - 1, j) + stirling2(m - 1, j - 1)
    _c[(m, j)] = v
    return v


def mono_to_ff(coefs):
    """coefficients in k^m -> coefficients in falling factorials (k)_j."""
    out = [0]*len(coefs)
    for m, a in enumerate(coefs):
        for j in range(m + 1):
            out[j] = out[j] + a*stirling2(m, j)
    return out


def nfrom(a, c, s):
    """real n with n(n+s) = c (n+a)(n+a+s)."""
    A = 1 - c; B = s - c*(2*a + s); Cq = -c*a*(a + s)
    return (-B + math.sqrt(B*B - 4*A*Cq))/(2*A)
