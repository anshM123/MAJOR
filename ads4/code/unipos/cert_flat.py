r"""Certificates C5a / C5b (one-dimensional, interval arithmetic).

C5a (flat-space lemma, bounded part):   fhat(b) = int_0^1 f(p) J_0(p b) dp > 0   for 0 <= b <= TAU2 = 10.85.
C5b (independent second proof of the curved theorem for M >= M0C = 45, region M beta <= TAU2):
        fhat(b) - (2591/1680) sinh^2( b / (2*45) ) > 0     for 0 <= b <= 10.85,
     which with |F_M(beta) - fhat(M beta)| <= (int f) * sinh^2(beta/2)/4  (Lemma 9) gives F_M(beta) > 0 for
     M >= 45, M beta <= 10.85.

fhat is evaluated from its Taylor series  fhat(b) = sum_k (-1)^k m_{2k} (b/2)^{2k} / (k!)^2   (exact rational
coefficients, m_{2k} = int_0^1 f p^{2k} <= 2591/420), truncated at k = K with the rigorous tail bound
(2591/420) * (b^{2K+2}/(2K+2)!)/(1 - b^2/((2K+3)(2K+4)))   (using binom(2k,k)/4^k <= 1), and
|fhat'(b)| <= (2/pi) * int p f = (2/pi)(283/120) =: L   (|J_0'(x)| = |(2/pi) int_0^{pi/2} sin(x sin th) sin th dth| <= 2/pi).
On a cell [l, u] with centre c:  fhat >= fhat(c) - L*max(c-l, u-c), and the subtracted term is increasing in b.
"""
import sys
from fractions import Fraction as Fr
from math import factorial, comb
from common import *

K = 36
M0C = 45
COEF = [frac_iv((-1) ** k * moment(2 * k) * Fr(comb(2 * k, k), 4 ** k) / factorial(2 * k)) for k in range(K + 1)]
L_IV = 2 * frac_iv(NORM_PF) / iv.mpf(iv.pi.a)          # upper bound for (2/pi)(283/120)


def mpf_to_fr(x):
    s, man, exp, bc = x._mpf_
    v = Fr(int(man)) * (Fr(2) ** int(exp))
    return -v if s else v


def fhat_lower(c):
    y = frac_iv(c) ** 2
    v = iv.mpf(0)
    for a in reversed(COEF):
        v = v * y + a
    b2 = y
    tail = frac_iv(NORM_F) * ipow(b2, K + 1) / factorial(2 * K + 2) / (1 - b2 / ((2 * K + 3) * (2 * K + 4)))
    return mpf_to_fr(lo(v - tail))


def exp_padded(x):
    """Rigorous enclosure of exp on an iv interval.  mpmath.iv.exp is NOT a guaranteed enclosure (an independent
    referee found misses of up to 4e-40 relative at 120 bits: directed rounding of a non-directed approximation).
    We therefore widen its output by a relative 2^-60, which exceeds the observed defect by twenty orders of magnitude."""
    e = iv.exp(x)
    return e * iv.mpf([1 - iv.mpf(2) ** (-60), 1 + iv.mpf(2) ** (-60)])


def sub_upper(u, M0):
    if M0 is None:
        return Fr(0)
    x = frac_iv(u) / (2 * M0)
    s = (exp_padded(x) - exp_padded(-x)) / 2
    return mpf_to_fr(hi(frac_iv(NORM_F) / 4 * s * s))


def march(M0, tag, out):
    Lq = mpf_to_fr(hi(L_IV))
    l = Fr(0)
    r = Fr(1, 32)
    n = 0
    vmin = None
    while l < TAU2:
        ok = False
        for _ in range(60):
            c, u = l + r, l + 2 * r
            margin = fhat_lower(c) - sub_upper(u, M0) - Lq * r
            if margin > 0:
                ok = True
                break
            r = r / 2
        if not ok:
            out.write(f"FAIL {tag} at b = {float(l)}\n")
            return False, n, vmin
        n += 1
        vmin = margin if vmin is None else min(vmin, margin)
        out.write(f"PASS {tag} cell [{float(l):.9f}, {float(u):.9f}]  lower bound on the cell = {float(margin):.6e}\n")
        l = u
        r = min(r * 2, Fr(1, 16))
    return True, n, vmin


def main():
    out = sys.stdout
    ok1, n1, v1 = march(None, "fhat>0", out)
    out.write(f"SUMMARY C5a: {n1} cells; min certified cell lower bound {float(v1) if v1 else None:.4e}\n")
    out.write(f"VERDICT C5a (fhat > 0 on [0, {float(TAU2)}]): {'PASS' if ok1 else 'FAIL'}\n")
    ok2, n2, v2 = march(M0C, f"fhat-(2591/1680)sinh^2(b/{2*M0C})>0", out)
    out.write(f"SUMMARY C5b: {n2} cells; min certified cell lower bound {float(v2) if v2 else None:.4e}\n")
    out.write(f"VERDICT C5b (fhat(b) > (2591/1680) sinh^2(b/{2*M0C}) on [0, {float(TAU2)}]): {'PASS' if ok2 else 'FAIL'}\n")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main())
