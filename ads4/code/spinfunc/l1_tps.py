"""
l1_tps: truncated multivariate power series (TPS) in NV variables to total degree DEG, complex coefficients, plus Gaussian
expectations of products of TPS (for the Laplace expansion of M6 step 1/2).

A TPS is a complex numpy array of length N over the monomials MONS (graded by degree).  Exact arithmetic up to DEG
(no finite differences), so Taylor coefficients of compositions of +,*,/,log,pow are exact to round-off.
"""
import itertools, math
import numpy as np

NV, DEG = 7, 6
MONS = []
for d_ in range(DEG + 1):
    for c in itertools.combinations_with_replacement(range(NV), d_):
        e = [0]*NV
        for i in c:
            e[i] += 1
        MONS.append(tuple(e))
N = len(MONS)
EXP = np.array(MONS, dtype=np.int64)
DEGS = EXP.sum(1)
_B = DEG + 1
_CODE = (EXP*(_B**np.arange(NV))).sum(1)
_ORD = np.argsort(_CODE)
_SC = _CODE[_ORD]


def _lookup(codes):
    return _ORD[np.searchsorted(_SC, codes)]


_I, _J = [], []
for a in range(N):
    bs = np.where(DEGS <= DEG - DEGS[a])[0]
    _I.append(np.full(len(bs), a)); _J.append(bs)
_I = np.concatenate(_I); _J = np.concatenate(_J)
_K = _lookup(_CODE[_I] + _CODE[_J])


def const(c):
    a = np.zeros(N, complex); a[0] = c
    return a


def var(i, c0=0.0):
    a = const(c0); a[1 + i] = 1.0          # degree-1 monomials are MONS[1..NV] in order
    return a


assert all(MONS[1 + i][i] == 1 for i in range(NV))


def mul(a, b):
    w = a[_I]*b[_J]
    return np.bincount(_K, w.real, N) + 1j*np.bincount(_K, w.imag, N)


def _series(a, coefs):
    """sum_n coefs[n] * delta^n with delta = a/a0 - 1 (no constant term)."""
    d = a/a[0]; d[0] = 0.0
    out = const(coefs[0]); p = const(1.0)
    for n in range(1, DEG + 1):
        p = mul(p, d)
        out = out + coefs[n]*p
    return out


def recip(a):
    return _series(a, [(-1.0)**n for n in range(DEG + 1)])/a[0]


def div(a, b):
    return mul(a, recip(b))


def log(a, branch0=None):
    """log TPS; constant term = principal log(a0) unless branch0 given."""
    s = _series(a, [0.0] + [(-1.0)**(n + 1)/n for n in range(1, DEG + 1)])
    s[0] = np.log(a[0]) if branch0 is None else branch0
    return s


def powr(a, p):
    c = [1.0]
    for n in range(1, DEG + 1):
        c.append(c[-1]*(p - n + 1)/n)
    return _series(a, c)*np.exp(p*np.log(a[0]))


def homog(a, d_):
    b = np.zeros(N, complex); m = DEGS == d_; b[m] = a[m]
    return b


def grad_hess(a):
    g = np.array([a[1 + i] for i in range(NV)])
    H = np.zeros((NV, NV), complex)
    for k in np.where(DEGS == 2)[0]:
        e = EXP[k]; nz = np.nonzero(e)[0]
        if len(nz) == 1:
            H[nz[0], nz[0]] = 2*a[k]
        else:
            H[nz[0], nz[1]] = H[nz[1], nz[0]] = a[k]
    return g, H


def linsub(z0, L):
    """TPS variables z_i = z0_i + sum_j L_ij xi_j."""
    out = []
    for i in range(NV):
        a = const(z0[i])
        a[1:1 + NV] = L[i]
        out.append(a)
    return out


# ---------------------------------------------------------------- Gaussian expectations over xi ~ N(0, 1) (formal, complex)
_DF = np.zeros(4*DEG + 2)
for n in range(0, 4*DEG + 1, 2):
    _DF[n] = float(np.prod(np.arange(n - 1, 0, -2))) if n > 0 else 1.0


def expect(a):
    """E[a(xi)] for a single TPS."""
    return np.sum(a*np.prod(_DF[EXP], axis=1))


def expect2(a, b, chunk=200):
    """E[a(xi) b(xi)] including the degree > DEG part of the product (exact for polynomials of degree <= DEG each)."""
    ia = np.nonzero(a)[0]; ib = np.nonzero(b)[0]
    Eb = EXP[ib]; bb = b[ib]
    tot = 0.0 + 0.0j
    for s in range(0, len(ia), chunk):
        sl = ia[s:s + chunk]
        S = EXP[sl][:, None, :] + Eb[None, :, :]
        M = np.prod(_DF[S], axis=2)
        tot += np.sum(a[sl][:, None]*bb[None, :]*M)
    return tot


def csym_cholesky(S):
    """complex-symmetric S = L L^T (no conjugation), L lower triangular."""
    n = S.shape[0]; L = np.zeros_like(S, dtype=complex)
    for j in range(n):
        s = S[j, j] - np.sum(L[j, :j]**2)
        L[j, j] = np.sqrt(s)
        for i in range(j + 1, n):
            L[i, j] = (S[i, j] - np.sum(L[i, :j]*L[j, :j]))/L[j, j]
    return L
