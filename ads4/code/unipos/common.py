r"""Common exact data and interval helpers for the uniform-positivity certificates.

  f(p)   = (1-p)^3 p (65 + 155 p + 47 p^2) = 65p - 40p^2 - 223p^3 + 259p^4 - 14p^5 - 47p^6
  phi(p) = f(p)/p
  S(tau) = int_0^1 f(p) sin(p tau) dp
         = 80/tau^3 + 6216/tau^5 + 33840/tau^7
           + (1602/tau^4 - 35520/tau^6) sin(tau) + (12384/tau^5 - 33840/tau^7) cos(tau)      (exact)

All rigorous computations use mpmath.iv (outward-rounded interval arithmetic) and exact Fractions.
"""
from fractions import Fraction as Fr
from math import comb, factorial
import mpmath
from mpmath import iv, mp

PREC = 120                      # bits, interval arithmetic working precision
iv.prec = PREC
mp.prec = PREC

# ---- exact data ---------------------------------------------------------------------------------
F_COEF = [0, 65, -40, -223, 259, -14, -47]          # f(p) = sum_i F_COEF[i] p^i
NORM_F = Fr(2591, 420)                              # int_0^1 f      (f >= 0 on [0,1])
NORM_PF = Fr(283, 120)                              # int_0^1 p f
D0 = {1: 65, 2: -80, 3: -1338, 4: 6216, 5: -1680, 6: -33840}     # f^{(j)}(0)
D1 = {3: -1602, 4: -12384, 5: -35520, 6: -33840}                 # f^{(j)}(1)  (f(1)=f'(1)=f''(1)=0)

TAU1 = Fr(96, 10)        # S > 0 on (0, TAU1]
TAU2 = Fr(1085, 100)     # S > 0 on [TAU2, oo)
TAU_INF = 25             # for tau >= 25 positivity of S is analytic
M1 = Fr(115, 100)        # small-M threshold (monotone-kernel argument)
MK = 6                   # kernel-positivity threshold: K_{f,M} > 0 certified for M <= MK


def moment(k):
    """int_0^1 f(p) p^k dp, exact."""
    return sum(Fr(c, i + k + 1) for i, c in enumerate(F_COEF))


assert moment(0) == NORM_F and moment(1) == NORM_PF


def fderiv(j, at):
    """exact f^{(j)}(at) for at in {0,1}."""
    c = F_COEF[:]
    for _ in range(j):
        c = [i * c[i] for i in range(1, len(c))]
    return sum(c) if at == 1 else (c[0] if c else 0)


for j in range(7):
    assert fderiv(j, 0) == D0.get(j, 0) and fderiv(j, 1) == D1.get(j, 0)


# ---- interval helpers ---------------------------------------------------------------------------
def frac_iv(q):
    """Enclosure of an exact rational (Fraction or int) by an mpmath interval."""
    q = Fr(q)
    return iv.mpf(q.numerator) / iv.mpf(q.denominator)


def lo(x):
    return mp.mpf(x.a)


def hi(x):
    return mp.mpf(x.b)


def ipow(x, n):
    """x**n for an interval x > 0 and integer n >= 0 (monotone, endpoints rounded outward)."""
    r = iv.mpf(1)
    for _ in range(n):
        r = r * x
    return r


def S_iv(tau):
    """Interval enclosure of S(tau) for an interval tau > 0 (closed form)."""
    t2 = tau * tau
    t3 = t2 * tau
    t4 = t2 * t2
    t5 = t4 * tau
    t6 = t4 * t2
    t7 = t6 * tau
    return (80 / t3 + 6216 / t5 + 33840 / t7
            + (1602 / t4 - 35520 / t6) * iv.sin(tau)
            + (12384 / t5 - 33840 / t7) * iv.cos(tau))


def K_low_iv(M, tau, mmax):
    r"""Interval enclosure of
        K_low^{(m)}(M,tau) = int_0^1 f(p) [1 + 2 sum_{k=1}^{m} e^{-2 pi k M p}] sin(p tau) dp
                           = S(tau) + 2 sum_{k=1}^m Im sum_j [f^{(j)}(0) - e^{-z_k} f^{(j)}(1)] z_k^{-j-1},
        z_k = 2 pi k M - i tau,
    for intervals M > 0, tau > 0.  (K_{f,M}(tau/M) >= K_low^{(m)}(M,tau) for every m >= 0, Lemma 4.)"""
    tot = S_iv(tau)
    c, s = iv.cos(tau), iv.sin(tau)
    for k in range(1, mmax + 1):
        a = 2 * iv.pi * k * M
        den = a * a + tau * tau
        wr, wi = a / den, tau / den                   # w = 1/z = (a + i tau)/(a^2+tau^2)
        ea = iv.exp(-a)
        er, ei = ea * c, ea * s                       # e^{-z} = e^{-a} (cos tau + i sin tau)
        pr, pi_ = wr, wi                              # w^1
        ar = iv.mpf(0); ai = iv.mpf(0)                # sum_j f^{(j)}(0) w^{j+1}
        br = iv.mpf(0); bi = iv.mpf(0)                # sum_j f^{(j)}(1) w^{j+1}
        for j in range(1, 7):
            pr, pi_ = pr * wr - pi_ * wi, pr * wi + pi_ * wr      # w^{j+1}
            ar = ar + D0[j] * pr; ai = ai + D0[j] * pi_
            if j in D1:
                br = br + D1[j] * pr; bi = bi + D1[j] * pi_
        im = ai - (er * bi + ei * br)
        tot = tot + 2 * im
    return tot
