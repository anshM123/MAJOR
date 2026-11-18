"""
Mack polynomials for identical external scalars in d = 3, in the conventions of CHMRS (arXiv:2106.10274) eq. (Bkt action)
and Caron-Huot--Mazac--Rastelli--Simmons-Duffin arXiv:2008.04931 eq. (qCal):

   M(S,T) ~ Q^n_{Delta,J}(T') / (S - tau - 2n),   tau = Delta - J,
   Q^n(T') = -sin^2[pi(tau-2Df)/2] K_n  Qpoly^n(T'),
   K_n = (Delta-1)_J Gamma(Delta+J) Gamma((tau-2Df+2n+2)/2)^2 / (2^{J-1} pi^2 n! (Delta-1/2)_n Gamma((Delta+J)/2)^4).

Qpoly^n is obtained from the explicit Mack polynomial of Costa--Goncalves--Penedones (arXiv:1209.4355, app. A, eq. (MackPol),
h = 3/2, Delta_12 = Delta_34 = 0, h + i nu = Delta) as
   Qpoly^n(T') = (-1)^J P_{Delta,J}(s, t),   s = 2Df - tau - 2n - T',  t = tau + 2n
(Costa's t is CHMRS's S, Costa's s = U - 2Df; the factor (-1)^J and 2^J come from the block normalisations
 G_Costa = (-1)^J 2^{-J} g_CHMRS; checked: J = 1 gives T' - 2Df + (Delta+2n-1)/2 as in 2008.04931).
In the B-hat sum rule T' = T + 2Df - tau - 2n, hence s = -T: Qpoly^n = (-1)^J P_{Delta,J}(-T, tau + 2n), independent of Df.

The polynomial P(s,t) is built exactly (Fractions; Delta is a binary float, converted exactly).
For numerics: coefficients of Qpoly^n as a polynomial in y = T - tau/2 - n (definite parity (-1)^J in y), computed per n in mpmath.
"""
from fractions import Fraction as Fr
from math import factorial
import functools
import mpmath as mp
import numpy as np


def _poch_poly(a, k):
    """(a)_k for a = linear poly given as dict {(i,j): coeff} in (s,t): returns poly dict."""
    out = {(0, 0): Fr(1)}
    for q in range(k):
        lin = dict(a)
        lin[(0, 0)] = lin.get((0, 0), Fr(0)) + q
        out = _mul(out, lin)
    return out


def _mul(p, q):
    o = {}
    for (i, j), c in p.items():
        for (k, l), d in q.items():
            key = (i + k, j + l)
            o[key] = o.get(key, Fr(0)) + c*d
    return {k: v for k, v in o.items() if v != 0}


def _add(p, q, f=Fr(1)):
    o = dict(p)
    for k, v in q.items():
        o[k] = o.get(k, Fr(0)) + f*v
    return {k: v for k, v in o.items() if v != 0}


def _rf(x, k):
    """rising factorial of a Fraction."""
    r = Fr(1)
    for q in range(k):
        r *= (x + q)
    return r


@functools.lru_cache(None)
def mack_P(Delta, J):
    """exact P_{Delta,J}(s,t) (Costa et al. normalisation, monic in s), d = 3, identical scalars. dict {(i,j): Fraction}: s^i t^j."""
    D = Fr(Delta)
    h = Fr(3, 2)
    hp = D            # h + i nu
    hm = 2*h - D      # h - i nu
    # delta_ij as linear polys in (s,t)
    d13 = {(1, 0): Fr(-1, 2)}
    d24 = {(1, 0): Fr(-1, 2)}
    d23 = {(1, 0): Fr(1, 2), (0, 1): Fr(1, 2)}
    d14 = {(1, 0): Fr(1, 2), (0, 1): Fr(1, 2)}
    al = [1 - (hp + J)/2, 1 - (hp + J)/2, 1 - (hm + J)/2, 1 - (hm + J)/2]
    den = _rf(hp - 1, J)*_rf(hm - 1, J)
    P = {}
    for r in range(J//2 + 1):
        # a_{J,r} = (-1)^r J! (h+J-1)_{-r} / (2^{2r} r! (J-2r)!),  (x)_{-r} = 1/((x-r)_r)
        aJr = Fr((-1)**r*factorial(J), 2**(2*r)*factorial(r)*factorial(J - 2*r))/_rf(h + J - 1 - r, r)
        pre = aJr*Fr(2)**(J + 2*r)*factorial(J - 2*r)/den
        # ((h + i nu - J - t)/2)_r ((h - i nu - J - t)/2)_r
        l1 = {(0, 0): (hp - J)/2, (0, 1): Fr(-1, 2)}
        l2 = {(0, 0): (hm - J)/2, (0, 1): Fr(-1, 2)}
        pr = _mul(_poch_poly(l1, r), _poch_poly(l2, r))
        inner = {}
        Kt = J - 2*r
        for k13 in range(Kt + 1):
            for k14 in range(Kt + 1 - k13):
                for k23 in range(Kt + 1 - k13 - k14):
                    k24 = Kt - k13 - k14 - k23
                    sgn = (-1)**(k13 + k24)
                    c = Fr(sgn, factorial(k13)*factorial(k14)*factorial(k23)*factorial(k24))
                    c *= _rf(al[0], J - r - (k13 + k14))
                    c *= _rf(al[1], J - r - (k23 + k24))
                    c *= _rf(al[2], J - r - (k13 + k23))
                    c *= _rf(al[3], J - r - (k14 + k24))
                    if c == 0:
                        continue
                    term = _mul(_mul(_poch_poly(d13, k13), _poch_poly(d14, k14)), _mul(_poch_poly(d23, k23), _poch_poly(d24, k24)))
                    inner = _add(inner, term, c)
        P = _add(P, _mul(pr, inner), pre)
    return P


def P_eval(P, s, t):
    return sum(c*s**i*t**j for (i, j), c in P.items())


def Qpoly_ycoeffs(Delta, J, ns, dps=40):
    """array [len(ns), J+1] (float64) of coefficients g_i(n) with Qpoly^n(T') = sum_i g_i(n) y^i,  y = T - tau/2 - n,
    T the B-hat variable (T' = T + 2Df - tau - 2n).  Qpoly^n = (-1)^J P(-T, tau + 2n)."""
    P = mack_P(float(Delta), int(J))
    tau = mp.mpf(Delta) - J
    out = np.zeros((len(ns), J + 1))
    with mp.workdps(dps):
        Pm = {k: mp.mpf(v.numerator)/v.denominator for k, v in P.items()}
        for a, n in enumerate(ns):
            n = mp.mpf(n)
            t = tau + 2*n
            c0 = tau/2 + n            # T = y + c0,  s = -T = -y - c0
            # P(s,t) = sum_{i,j} p_ij s^i t^j ; s^i = (-1)^i (y + c0)^i
            g = [mp.mpf(0)]*(J + 1)
            for (i, j), pij in Pm.items():
                base = pij*(-1)**i*t**j
                for l in range(i + 1):
                    g[l] += base*mp.binomial(i, l)*c0**(i - l)
            for l in range(J + 1):
                out[a, l] = float((-1)**J*g[l])
    return out


if __name__ == '__main__':
    # checks: monic in s; J = 1 against 2008.04931; parity in y
    for J in range(0, 5):
        P = mack_P(37.25, J)
        lead = P.get((J, 0), 0)
        print("J", J, "leading s^J coeff", lead, " #terms", len(P))
    D = 37.25; J = 1
    for n in (0, 3):
        g = Qpoly_ycoeffs(D, J, [n])[0]
        # 2008.04931: Q^n_{D,1}(T') = T' - 2Df + (D+2n-1)/2, with T' = T + 2Df - tau - 2n:  = T - tau/2 - n + [(D+2n-1)/2 - tau/2 - n] = y + 0
        print("J=1 n=%d y-coeffs" % n, g, "(expected [0, 1])")
    for J in (2, 3, 4, 6):
        g = Qpoly_ycoeffs(D, J, [0, 5, 50.5])
        print("J", J, "parity check (odd coeffs for even J should vanish):", np.abs(g[:, (J + 1) % 2::2]).max(), " leading", g[:, J])
