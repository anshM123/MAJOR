"""[certified, Arb] The naive flat limit F_U of Omega_g3 (C_k -> B'_k) is NOT positive on heavy states:
   m^4 F_U(++, J=2, m) = int_0^1 [psi2(2-uq) + u psi3](1-uq)^-2 P_2(1-2uq) dp + (10*6 - 20) u^3 < 0  at m = 213/200, 6/5,
with u = 1/m^2, q = p^2.  Rigorous integration: flint.acb.integral (Arb, ball arithmetic; the integrand is analytic on a
neighbourhood of [0,1] since the pole p = m lies outside).
"""
import flint
from fractions import Fraction as Fr
flint.ctx.prec = 128
acb, arb = flint.acb, flint.arb

def FU_pp_J2(m):
    u = acb(flint.fmpq(m.denominator ** 2, m.numerator ** 2))
    def f(p, analytic):
        q = p * p
        psi2 = (1 - p) ** 3 * p * (65 + 155 * p + 47 * q)
        psi3 = 2 * (1 - p) ** 3 * (81 - 1927 * p + 4601 * q)
        x = 1 - 2 * u * q
        P2 = (3 * x * x - 1) / 2
        return (psi2 * (2 - u * q) + u * psi3) / (1 - u * q) ** 2 * P2
    I = acb.integral(f, 0, 1)
    return I + (10 * 6 - 20) * u ** 3

ok = True
for m in [Fr(213, 200), Fr(6, 5)]:
    v = FU_pp_J2(m)
    neg = bool(v.real < 0) and bool(abs(v.imag) < arb(1e-20))
    ok &= neg
    print('m^4 F_U(++, J=2, m=%s) = %s  -> %s' % (m, v.real.str(12), 'NEGATIVE (certified)' if neg else 'not certified'))
print('VERDICT cert_FU_negative: %s' % ('PASS' if ok else 'FAIL'))
