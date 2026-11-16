"""Ingredients for Lemma IR (I(lambda), the phi bounds).  NB (referee round 2): the Delta_IR table printed here, for
m_IR <= M/100 with the scaling-limit form of E, is only an ESTIMATE; the proved version (m_IR <= M/20000, E on the
certified region Omega_eps, finite-mass corrections bounded) is cert_eps.py.

Under E: Delta_IR <= (1+eta) (16/k(eps)) int_0^{1/mu} [fhat_eps(b)]_- b log(1/(mu b)) db,  mu = lambda*eps,
fhat_eps(b) = int_eps^1 psi2(p) J0(bp) dp.  Using fhat >= 0 (UNIFORM_POSITIVITY.md, Theorem 2),
psi2 = p phi(p), 0 <= 65 - phi(p) <= 40p + 223p^2 on [0,1/100] (checked here exactly) and |J0| <= 1:
   [fhat_eps(b)]_- <= [65 eps J1(eps b)/b]_+ + 14 eps^3,
   int_0^{1/mu} (...) b log(1/(mu b)) db <= 65 I(lambda) + 3.5 eps/lambda^2,
   I(lambda) = int_0^{1/lambda} [J1(x)]_+ log(1/(lambda x)) dx
             = sum_k (-1)^k X^{2k+2} / (2^{2k+1} k! (k+1)! (2k+2)^2),  X = 1/lambda < j_{1,1}  (J1 > 0 on (0,X]).
Prints the resulting Delta_IR bound and the corrected constant B_E = B + Delta for eps <= 1/100, eta = 0.
"""
import sympy as sp
import flint
from fractions import Fraction as Fr
from math import factorial
from common import set_prec

set_prec(200)
A = flint.arb
def fq(x): x = Fr(x); return A(flint.fmpq(x.numerator, x.denominator))

ok = True
# 0 <= 65 - phi <= 40p + 223p^2 on [0, 1/100]
p = sp.Symbol('p')
phi = sp.expand((1 - p)**3*(65 + 155*p + 47*p**2))
# divide out the zeros at p = 0 so that the root counts cover the closed interval [0, 1/100]
d1 = sp.Poly(sp.cancel(sp.expand(65 - phi)/p), p)                         # 40 + 223p - 259p^2 + ...
d2 = sp.Poly(sp.cancel(sp.expand(40*p + 223*p**2 - (65 - phi))/p**3), p)   # 259 - 14p - 47p^2
c1 = sp.count_roots(d1, 0, sp.Rational(1, 100)) == 0 and d1.eval(0) > 0
c2 = sp.count_roots(d2, 0, sp.Rational(1, 100)) == 0 and d2.eval(0) > 0
MSG = ['0 <= 65-phi <= 40p+223p^2 on [0,1/100] (zeros at p=0 divided out): %s' % (c1 and c2)]
ok &= c1 and c2
eps0 = Fr(1, 100)
MSG.append('int_0^eps (40p+223p^2) p dp <= 14 eps^3 for eps<=1/100: %s' % (Fr(40, 3) + Fr(223, 4) * eps0 <= 14))
ok &= Fr(40, 3) + Fr(223, 4) * eps0 <= 14

def J1_positive(X, n=400):
    """J1 > 0 on (0, X]: near 0 J1(x)/x >= 1/2 - x^2/16 > 0 for x <= 1; Arb on cells of [1, X]."""
    if X <= 1: return True
    for i in range(n):
        lo = 1 + (Fr(X) - 1) * i / n; hi = 1 + (Fr(X) - 1) * (i + 1) / n
        x = fq(lo).union(fq(hi))
        if not bool(x.bessel_j(1) > 0): return False
    return True

def I_lambda(lam, K=60):
    X = 1 / Fr(lam)
    tot = Fr(0); terms = []
    for k in range(K):
        t = Fr((-1) ** k) * X ** (2 * k + 2) / (2 ** (2 * k + 1) * factorial(k) * factorial(k + 1) * (2 * k + 2) ** 2)
        tot += t; terms.append(t)
    # alternating tail with ratio < 1/2 after k > X^2: |tail| <= 2|term_K|
    nxt = X ** (2 * K + 2) / (2 ** (2 * K + 1) * factorial(K) * factorial(K + 1) * (2 * K + 2) ** 2)
    assert X * X / (4 * (K + 1) * (K + 2)) < Fr(1, 2)
    return fq(tot) + (-fq(2 * nxt)).union(fq(2 * nxt))

k0 = Fr(330151, 12012)
# k_lo for eps <= 1/100 (as in cert_low): k0 - 2 int_0^{1/100} psi3 p^4
psi3 = sp.expand(2*(1 - p)**3*(81 - 1927*p + 4601*p**2))
klo = k0 - 2 * Fr(str(sp.integrate(psi3*p**4, (p, 0, sp.Rational(1, 100)))))
psi2 = sp.expand((1 - p)**3*p*(65 + 155*p + 47*p**2))
khi = k0 + 2 * Fr(str(sp.integrate(psi2*p**6, (p, 0, sp.Rational(1, 100)))))
e = eps0
ell = 40*e + Fr(223, 2)*e**2 - Fr(259, 3)*e**3 + Fr(7, 2)*e**4 + Fr(47, 5)*e**5
Bexact = 16 * (ell - Fr(1171, 15)) / khi
Bconst = fq(Bexact)
if __name__ == '__main__':
    for line in MSG: print(line)
    print('B (eps<=1/100) =', Bconst.str(10))
    print('lambda   X=1/lambda   I(lambda)            Delta_IR bound (eta=0)    B_E = B + Delta')
    for lam in [Fr(2, 7), Fr(1, 2), Fr(1), Fr(2), Fr(4)]:
        X = 1 / lam
        pos = J1_positive(X)
        ok &= pos and X < Fr(383, 100)
        I = I_lambda(lam)
        D = 16 / fq(klo) * (65 * I + fq(Fr(7, 2)) * fq(eps0) / fq(lam) ** 2)
        print('%-7s  %-11s  %-20s %-25s %s' % (lam, float(X), I.str(8), D.str(6), (Bconst + D).str(6)))
    print('VERDICT cert_ir: %s' % ('PASS' if ok else 'FAIL'))
