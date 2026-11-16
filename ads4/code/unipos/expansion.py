r"""Exact curvature expansion of F_M(beta) (Lemma 6 of UNIFORM_POSITIVITY.md).

With a = (beta+t)/2, c = (beta-t)/2:   cosh(beta) - cosh(t) = ((beta^2 - t^2)/2) * E,  E = (sinh a / a)(sinh c / c),
    (sinh a / a)^{-1/2} = sum_k d_k a^{2k},
    E^{-1/2}|_{t = beta s} = sum_n beta^{2n} q_n(s),   q_n(s) = sum_{k+l=n} d_k d_l ((1+s)/2)^{2k} ((1-s)/2)^{2l},
    F_M(beta) = sum_{n>=0} beta^{2n} G_n(M beta),
    G_n(b) = (2/pi) int_0^{pi/2} ftil(b sin th) q_n(sin th) dth = sum_{k>=0} g_{n,k} b^{2k},
    g_{n,k} = (-1)^k m_{2k}/(2k)! * sum_j q_{n,j} W_{k+j},   W_m = binom(2m,m)/4^m,   m_{2k} = int_0^1 f(p) p^{2k} dp.
Everything here is exact rational arithmetic.
"""
from fractions import Fraction as Fr
from math import comb, factorial
from common import moment


def d_coeffs(nmax):
    """d_k, k = 0..nmax: Taylor coefficients (in y = a^2) of (sinh a / a)^{-1/2}  (J.C.P. Miller recurrence)."""
    u = [Fr(1, factorial(2 * k + 1)) for k in range(nmax + 1)]
    al = Fr(-1, 2)
    d = [Fr(1)]
    for n in range(1, nmax + 1):
        d.append(sum((al * k - (n - k)) * u[k] * d[n - k] for k in range(1, n + 1)) / n)
    return d


def mu_coeffs(nmax):
    """mu_n = [y^{2n}] (y / sin y)  (>= 0);  computed as the Cauchy square of |d_k| = [y^{2k}] (y/sin y)^{1/2}."""
    d = d_coeffs(nmax)
    assert all((-1) ** k * d[k] > 0 for k in range(nmax + 1))
    return [sum(abs(d[k]) * abs(d[n - k]) for k in range(n + 1)) for n in range(nmax + 1)]


def polymul(a, b):
    r = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def q_polys(nmax):
    """q_n as lists of coefficients of s^{2j}, j = 0..n."""
    d = d_coeffs(nmax)
    plus = [[Fr(1)]]      # ((1+s)/2)^{2k}
    minus = [[Fr(1)]]     # ((1-s)/2)^{2k}
    p2 = [Fr(1, 4), Fr(1, 2), Fr(1, 4)]
    m2 = [Fr(1, 4), Fr(-1, 2), Fr(1, 4)]
    for k in range(nmax):
        plus.append(polymul(plus[-1], p2))
        minus.append(polymul(minus[-1], m2))
    qs = []
    for n in range(nmax + 1):
        tot = [Fr(0)] * (2 * n + 1)
        for k in range(n + 1):
            pr = polymul(plus[k], minus[n - k])
            c = d[k] * d[n - k]
            for i, x in enumerate(pr):
                tot[i] += c * x
        assert all(tot[i] == 0 for i in range(1, 2 * n + 1, 2))      # q_n is even in s
        qs.append(tot[0::2])
    return qs


def g_coeffs(nmax, kmax):
    """g[n][k], exact."""
    qs = q_polys(nmax)
    W = [Fr(comb(2 * m, m), 4 ** m) for m in range(nmax + kmax + 1)]
    g = []
    for n in range(nmax + 1):
        row = []
        for k in range(kmax + 1):
            row.append((-1) ** k * moment(2 * k) / factorial(2 * k) * sum(qs[n][j] * W[k + j] for j in range(n + 1)))
        g.append(row)
    return g
