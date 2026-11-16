"""Certificate C (crude region): 1 <= m <= 2 (u in [1/4,1]) and J >= JP.

m^10 F = Phi + m^6 * int Psi D dp, m^6 <= 64.
 ++ : |P_J| <= 1, |Psi_+| <= 2 psi2 + |psi3|            => need 10JJ-5   > 64 * (2 int psi2 + int|psi3|)
 +- : |dtilde44| <= min( binom(J+4,8), (1-p^2)^{-4} )   (|d^J_44|<=1, d44=(1-y)^4 dtilde44, Szego 7.32.2)
      |Psi_-| <= 2 psi2 + |psi3| =: (1-p)^3 g(p)
      U(J) := int_0^1 (1-p)^3 g min(B,(1-p^2)^-4) <= int_0^{1/2} g(1-p)^3(1-p^2)^-4
              + gamma*(1/4) log B + g(1)/(4(1+p*)^4),  gamma = sup_[1/2,1] g/(1+p)^4,  p* = sqrt(1-B^-1/4) >= 1/2
      => need 10JJ-210 > 64 U(J); the difference is increasing in J, so J = JP suffices.
All bounds by interval arithmetic (arb) on subintervals; no floating point decisions.
"""
import sys
import flint
from fractions import Fraction as Fr
from common import PSI2, PSI3, peval, mom, set_prec

set_prec(200)
A = flint.arb

def ev(poly, x):
    r = A(0)
    for c in reversed(poly):
        r = r * x + A(flint.fmpq(c.numerator, c.denominator))
    return r

def g_hat(x):  # 2 p (65+155p+47p^2) + 2|81-1927p+4601p^2|  (interval x)
    q = 81 - 1927 * x + 4601 * x * x
    aq = q.abs_upper()  # upper bound of |q| over the ball
    return 2 * x * (65 + 155 * x + 47 * x * x) + 2 * A(aq)

def sup_on(f, a, b, ncell):
    best = None
    for i in range(ncell):
        lo = Fr(a) + (Fr(b) - Fr(a)) * i / ncell; hi = Fr(a) + (Fr(b) - Fr(a)) * (i + 1) / ncell
        x = A(flint.fmpq(lo.numerator, lo.denominator)).union(A(flint.fmpq(hi.numerator, hi.denominator)))
        v = f(x).upper()
        best = v if best is None else best.max(v)
    return best

def integral_upper(f, a, b, ncell):
    tot = A(0)
    for i in range(ncell):
        lo = Fr(a) + (Fr(b) - Fr(a)) * i / ncell; hi = Fr(a) + (Fr(b) - Fr(a)) * (i + 1) / ncell
        x = A(flint.fmpq(lo.numerator, lo.denominator)).union(A(flint.fmpq(hi.numerator, hi.denominator)))
        tot += A(f(x).upper()) * A(flint.fmpq((hi - lo).numerator, (hi - lo).denominator))
    return tot.upper()

def main():
    JP = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    ok_all = True
    # ++
    int_psi2 = mom(PSI2, 0)
    abs_psi3 = integral_upper(lambda x: A(ev(PSI3, x).abs_upper()), 0, 1, 4000)
    bound_pp = 64 * (2 * A(flint.fmpq(int_psi2.numerator, int_psi2.denominator)) + A(abs_psi3))
    lhs = A(10 * JP * (JP + 1) - 5)
    ok = bool(lhs > bound_pp)
    print('C ++ J>=%d: 10JJ-5=%s  > 64*(2int psi2+int|psi3|) <= %s : %s' % (JP, lhs.str(5), bound_pp.str(5), 'PASS' if ok else 'FAIL'))
    ok_all &= ok
    # +-
    c_a = integral_upper(lambda x: g_hat(x) / ((1 - x) * (1 + x) ** 4), 0, Fr(1, 2), 2000)
    gamma = sup_on(lambda x: g_hat(x) / (1 + x) ** 4, Fr(1, 2), 1, 4000)
    g1 = sup_on(g_hat, Fr(1, 2), 1, 4000)      # sup of g on [1/2,1] (contains [p*,1])
    J = JP
    from math import comb
    Bn = A(comb(J + 4, 8))
    assert bool(Bn > A(4) ** 4 / A(3) ** 4)      # => p* >= 1/2
    # second piece: B int_{p*}^1 (1-p)^3 g <= sup g * B (1-p*)^4/4 = sup g/(4(1+p*)^4) <= sup g/(4 (3/2)^4)
    U = A(c_a) + A(gamma) * Bn.log() / 4 + A(g1) / (4 * (A(3) / 2) ** 4)
    lhs = A(10 * J * (J + 1) - 210)
    ok = bool(lhs > 64 * U)
    print('C +- J>=%d: c_a<=%s gamma<=%s g(1)<=%s  U(J)<=%s ; 10JJ-210=%s > 64U=%s : %s' % (
        JP, A(c_a).str(6), A(gamma).str(6), A(g1).str(6), U.str(6), lhs.str(6), (64 * U).str(6), 'PASS' if ok else 'FAIL'))
    # monotonicity: d/dJ [10J(J+1)] = 20J+10 ; d/dJ [64*gamma/4*log binom(J+4,8)] <= 16*gamma*8/(J-3)
    mono = bool(A(20 * JP + 10) > 16 * A(gamma) * 8 / A(JP - 3))
    print('C +- monotonicity for J>=%d: %s' % (JP, 'PASS' if mono else 'FAIL'))
    ok_all &= ok and mono
    print('VERDICT cert_C: %s' % ('PASS' if ok_all else 'FAIL'))

if __name__ == '__main__':
    main()
