r"""Certificate C3:  K_low^{(4)}(M, tau) > 0  for all (M, tau) in [M1, MK] x [TAU1, TAU2] = [1.15, 6] x [9.6, 10.85].

K_low^{(m)}(M,tau) = int_0^1 f(p) [1 + 2 sum_{k=1}^m e^{-2 pi k M p}] sin(p tau) dp   (explicit elementary closed form,
see common.K_low_iv) is a lower bound for the Mehler kernel K_{f,M}(t) = int_0^1 f(p) coth(pi M p) sin(M p t) dp at
t = tau/M (Lemma 5 of UNIFORM_POSITIVITY.md).  Adaptive bisection; every box is evaluated with outward-rounded interval
arithmetic (natural interval extension of the closed form) and is reported as PASS when the lower endpoint is > 0.
"""
import sys
from fractions import Fraction as Fr
from common import *

MMAX_TERMS = 4


def main():
    out = sys.stdout
    nM0, nT0 = 97, 25                     # initial grid: dM = 0.05, dtau = 0.05
    boxes = []
    for i in range(nM0):
        for j in range(nT0):
            boxes.append((M1 + Fr(MK - M1) * Fr(i, nM0), M1 + Fr(MK - M1) * Fr(i + 1, nM0),
                          TAU1 + (TAU2 - TAU1) * Fr(j, nT0), TAU1 + (TAU2 - TAU1) * Fr(j + 1, nT0), 0))
    npass = 0
    kmin = None
    ok = True
    while boxes:
        Ma, Mb, ta, tb, depth = boxes.pop()
        Mi = iv.mpf([frac_iv(Ma).a, frac_iv(Mb).b])
        ti = iv.mpf([frac_iv(ta).a, frac_iv(tb).b])
        v = K_low_iv(Mi, ti, MMAX_TERMS)
        if v.a > 0:
            npass += 1
            kmin = lo(v) if kmin is None else min(kmin, lo(v))
            out.write(f"PASS box M in [{float(Ma):.6f}, {float(Mb):.6f}] tau in [{float(ta):.6f}, {float(tb):.6f}]  K_low >= {mp.nstr(lo(v), 6)}\n")
        elif depth >= 30:
            ok = False
            out.write(f"FAIL box M in [{float(Ma):.6f}, {float(Mb):.6f}] tau in [{float(ta):.6f}, {float(tb):.6f}]  enclosure {v}\n")
        else:
            Mm, tm = (Ma + Mb) / 2, (ta + tb) / 2
            for (x0, x1) in ((Ma, Mm), (Mm, Mb)):
                for (y0, y1) in ((ta, tm), (tm, tb)):
                    boxes.append((x0, x1, y0, y1, depth + 1))
    out.write(f"SUMMARY C3: {npass} boxes PASS; min certified lower bound of K_low^({MMAX_TERMS}) = {mp.nstr(kmin, 6)}\n")
    out.write(f"VERDICT C3 (K_low > 0 on [{float(M1)}, {MK}] x [{float(TAU1)}, {float(TAU2)}]): {'PASS' if ok else 'FAIL'}\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
