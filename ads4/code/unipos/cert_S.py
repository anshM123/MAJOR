r"""Certificate C1:  S(tau) = int_0^1 f(p) sin(p tau) dp > 0  on  [pi, TAU1] u [TAU2, 25],
and exact check of the analytic tail inequality used for tau >= 25.

Method: |S'(tau)| = |int_0^1 p f(p) cos(p tau) dp| <= int_0^1 p f = 283/120 =: L.  If c is a point with an interval
lower bound s <= S(c), then S > 0 on (c - s/L, c + s/L).  The cells [l, u] below are chosen with
max(c-l, u-c) * L < s (checked in exact rational arithmetic) and cover the target interval without gaps.
"""
import sys
from fractions import Fraction as Fr
from common import *

L = NORM_PF
DY = 2 ** 40


def dyadic_floor(x):
    return Fr(int(x * DY), DY)          # x > 0 : floor to a dyadic rational


def mpf_to_fr(x):
    """exact conversion of an mpf to a Fraction."""
    s, man, exp, bc = x._mpf_
    v = Fr(int(man)) * (Fr(2) ** int(exp))
    return -v if s else v


def march(a, b, tag, out):
    """certify S > 0 on [a,b]; returns (ok, ncells, min lower bound)."""
    l = Fr(a)
    r = Fr(1, 64)
    n = 0
    smin = None
    while l < b:
        ok = False
        for _ in range(60):
            c = l + r
            sl = mpf_to_fr(lo(S_iv(frac_iv(c))))
            if sl > L * r:
                ok = True
                break
            r = r / 2
        if not ok:
            out.write(f"FAIL {tag} at tau = {float(l)}\n")
            return False, n, smin
        r2 = dyadic_floor(Fr(95, 100) * sl / L)
        assert r2 > 0 and L * r2 < sl and L * (c - l) < sl
        u = c + r2
        n += 1
        smin = sl if smin is None else min(smin, sl)
        out.write(f"PASS {tag} cell [{float(l):.9f}, {float(u):.9f}]  S(c) >= {float(sl):.6e}  L*halfwidth = {float(L*max(r2, c-l)):.6e}\n")
        l = u
        r = r2
    return True, n, smin


def main():
    out = sys.stdout
    ok_all = True
    # left end: an exact rational slightly below pi  (S > 0 on (0, pi] trivially, so any a <= pi works)
    a = Fr(314, 100)
    assert frac_iv(a).b < iv.pi.a
    ok1, n1, s1 = march(a, TAU1, "S>0 [3.14, 9.6]", out)
    ok2, n2, s2 = march(TAU2, Fr(TAU_INF), "S>0 [10.85, 25]", out)
    # analytic tail: for tau >= 25,  tau^6 S(tau) >= 80 tau^3 - 1602 tau^2 - 6168 tau - 35520 =: c(tau) > 0
    t = TAU_INF
    c_val = 80 * t ** 3 - 1602 * t ** 2 - 6168 * t - 35520
    c_der = 240 * t ** 2 - 3204 * t - 6168          # c'(25);  c'' = 480 tau - 3204 > 0 for tau >= 25
    ok3 = c_val > 0 and c_der > 0 and 480 * t - 3204 > 0
    out.write(f"{'PASS' if ok3 else 'FAIL'} tail polynomial: c(25) = {c_val} > 0, c'(25) = {c_der} > 0, c''(25) = {480*t-3204} > 0\n")
    ok_all = ok1 and ok2 and ok3
    out.write(f"SUMMARY C1: cells {n1} + {n2}; min certified lower bound of S on [3.14,9.6]: {float(s1):.4e}; on [10.85,25]: {float(s2):.4e}\n")
    out.write(f"VERDICT C1 (S > 0 on (0,9.6] u [10.85,oo)): {'PASS' if ok_all else 'FAIL'}\n")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
