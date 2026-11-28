"""Exact one-loop sigma bubble on S^1_beta x S^2_r (beta = 1) at the N = infinity saddle.

Pi_{n,J} = (1/(4 pi r^2)) sum_{l1,l2} (2l1+1)(2l2+1) (l1 l2 J;000)^2 K_n(E_l1, E_l2),
E_l^2 = (l+1/2)^2/r^2 + m(r)^2  (conformal coupling), m(r) from the S^2 gap equation <phi^2>_ren = 0.
K_n = T sum_{n1} 1/((w1^2+E1^2)((w-w1)^2+E2^2)) in closed form.
The l1 sum is truncated at L and Richardson-extrapolated (tail ~ 1/L).
"""
import numpy as np
from scipy.special import gammaln
from mpmath import mp, mpf, sqrt as msqrt, exp as mexp, quad as mquad, nsum, inf, findroot, log as mlog, pi as mpi

mp.dps = 25
M0 = float(2 * mlog((1 + msqrt(5)) / 2))


def G_s2(m, r):
    """renormalised <phi^2>/N on S^1 x S^2 (scout's Abel-Plana formula)."""
    M = m * r
    vac = (-M + 2 * mquad(lambda t: t / (msqrt(M ** 2 - t ** 2) * (mexp(2 * mpi * t) + 1)), [0, M])) / (4 * mpi * r)
    th = nsum(lambda l: (2 * l + 1) / (mexp(msqrt((l + 0.5) ** 2 / r ** 2 + m ** 2)) - 1)
              / msqrt((l + 0.5) ** 2 / r ** 2 + m ** 2), [0, inf]) / (4 * mpi * r ** 2)
    return vac + th


_mcache = {}


def m_of_r(r):
    if r not in _mcache:
        _mcache[r] = float(findroot(lambda mm: G_s2(mm, mpf(r)), mpf(M0)))
    return _mcache[r]


def log3j2(l1, l2, J):
    """log of (l1 l2 J; 0 0 0)^2 for arrays with l1+l2+J even and triangle satisfied."""
    g = (l1 + l2 + J) // 2
    a = gammaln(2 * g - 2 * l1 + 1) + gammaln(2 * g - 2 * l2 + 1) + gammaln(2 * g - 2 * J + 1) - gammaln(2 * g + 2)
    b = gammaln(g + 1) - gammaln(g - l1 + 1) - gammaln(g - l2 + 1) - gammaln(g - J + 1)
    return a + 2 * b


def Kn(w, E1, E2):
    with np.errstate(over='ignore'):
        n1 = 1 / np.expm1(E1)
        n2 = 1 / np.expm1(E2)
    d = E1 - E2
    s = E1 + E2
    with np.errstate(invalid='ignore', divide='ignore'):
        t2 = np.where(np.abs(d) > 1e-9 if w == 0 else True, (n2 - n1) * 2 * d / (w * w + d * d),
                      2 * np.exp(-E1) / np.expm1(-E1) ** 2)
    return ((1 + n1 + n2) * 2 * s / (w * w + s * s) + t2) / (4 * E1 * E2)


def bubble_partial(n, J, r, L, m=None):
    """partial sums S(L') for L' in L (list), l1 <= L'."""
    m = m_of_r(r) if m is None else m
    w = 2 * np.pi * n
    Lmax = max(L)
    l1 = np.arange(0, Lmax + 1)
    tot = np.zeros(Lmax + 1)
    # loop over offset l2 = l1 + d, d in [-J, J], parity of l1+l2+J even -> d = J mod 2
    for d in range(-J, J + 1, 2):
        l2 = l1 + d
        ok = (l2 >= 0) & (l2 >= np.abs(l1 - J) ) & (l2 <= l1 + J)
        L1 = l1[ok]; L2 = l2[ok]
        E1 = np.sqrt((L1 + 0.5) ** 2 / r ** 2 + m * m)
        E2 = np.sqrt((L2 + 0.5) ** 2 / r ** 2 + m * m)
        val = (2 * L1 + 1) * (2 * L2 + 1) * np.exp(log3j2(L1, L2, J)) * Kn(w, E1, E2)
        tot[ok] += val
    cs = np.cumsum(tot) / (4 * np.pi * r * r)
    return np.array([cs[x] for x in L])


def bubble(n, J, r, L0=4000, m=None, k=5):
    """polynomial extrapolation in 1/L of the truncated l1-sum (tail = a/L + b/L^2 + ...)."""
    Ls = [int(L0 * 2 ** (i / 2)) for i in range(k + 1)]
    S = bubble_partial(n, J, r, Ls, m)
    x = 1.0 / np.array(Ls, dtype=float)
    est = []
    for kk in (k, k + 1):
        A = np.vander(x[-kk:], kk, increasing=True)
        est.append(np.linalg.solve(A, S[-kk:])[0])
    return est[1], abs(est[1] - est[0])
