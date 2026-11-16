"""t3: the graviton term of -Omega_light in closed form, its large-Lambda expansion, and the conditional AdS4 constants
(ADS4_LIGHT.md sections 2.5 and 4).

-Omega_light|_Einstein = Lam^4 int_0^1 psi2(p) 16 pi G/(Lam^2 p^2 + 9/4) dp = 16 pi G Lam^2 L_AdS(Lam),
L_AdS(Lam) = int_0^1 psi2(p)/(p^2 + a^2) dp,  a = 3/(2 Lam).
Exact: psi2(p) = (p^2 + a^2) Q_a(p) + alpha(a) p + beta(a)  (polynomial division), so
L_AdS = int_0^1 Q_a + (alpha/2) log(1 + 1/a^2) + (beta/a) arctan(1/a).   [D, sympy; checked against quadrature]
"""
import sympy as sp
import mpmath as mp

q, a = sp.symbols('q a', positive=True)
psi2 = sp.expand((1 - q) ** 3 * q * (65 + 155 * q + 47 * q ** 2))
psi3 = sp.expand(2 * (1 - q) ** 3 * (81 - 1927 * q + 4601 * q ** 2))


def L_exact():
    Qa, rem = sp.div(sp.Poly(psi2, q), sp.Poly(q ** 2 + a ** 2, q))
    rem = rem.as_expr()
    alpha = sp.expand(rem.coeff(q, 1))
    beta = sp.expand(rem.coeff(q, 0))
    L = sp.integrate(Qa.as_expr(), (q, 0, 1)) + alpha / 2 * sp.log(1 + 1 / a ** 2) + beta / a * sp.atan(1 / a)
    return sp.simplify(L), alpha, beta


def main():
    L, alpha, beta = L_exact()
    print('alpha(a) =', sp.factor(alpha), ';  beta(a) =', sp.factor(beta))
    ser = sp.series(L, a, 0, 4).removeO()
    ser = sp.simplify(sp.expand(ser))
    print('L_AdS small-a expansion:', sp.expand(ser))
    k0 = sp.Rational(330151, 12012)
    Lam = sp.symbols('Lambda', positive=True)
    serL = sp.expand(ser.subs(a, sp.Rational(3, 2) / Lam))
    print('L_AdS(Lam) =', serL, '+ O(Lam^-4 log Lam)')
    # numerical check vs quadrature
    Lf = sp.lambdify(a, L, 'mpmath')
    for Lv in [3, 10, 100, 1000]:
        av = mp.mpf(1.5) / Lv
        quad = mp.quad(lambda x: sp.lambdify(q, psi2, 'mpmath')(x) / (x * x + av * av), [0, av, 10 * av, 1])
        print('  Lam=%5d  exact=%s  quad=%s   3-term asympt=%s' % (Lv, mp.nstr(Lf(av), 15), mp.nstr(quad, 15),
                                                                 mp.nstr(sp.N(serL.subs(Lam, Lv)), 15)))
    A = 16 * 65 / k0
    B = 16 * (65 * sp.log(sp.Rational(2, 3)) - sp.Rational(1171, 15)) / k0
    C1 = sp.simplify(16 * serL.coeff(Lam, -1) / k0)
    print('\nConditional bound (k_AdS = k(0) + O(Lam^-2), all other light terms as in ADS4_LIGHT.md section 2):')
    print('  |g3|^2 Lam^8 <= A log Lam + B + C1/Lam + O(log Lam/Lam^2):  A = 1040/k(0) = %s,  B = %s,  C1 = %s = %s'
          % (sp.N(A, 10), sp.N(B, 10), C1, sp.N(C1, 8)))
    print('  flat (3.4a)-type constant for comparison: 16(-1171/15)/k(0) = %s; AdS shift 1040 log(2/3)/k(0) = %s'
          % (sp.N(16 * sp.Rational(-1171, 15) / k0, 10), sp.N(A * sp.log(sp.Rational(2, 3)), 10)))
    print('  CFT3 form: t4^2 + alpha_T^2 = 57600 |g3|^2 (AdS units):  t4^2 + alpha_T^2 <= 57600 (A log Lam + B + C1/Lam)/Lam^8')
    for Dg in [10, 20, 50, 100, 1000]:
        Lv = mp.sqrt((Dg - 2) * (Dg - 1))
        rhs = float(A) * mp.log(Lv) + float(B) + float(C1) / Lv
        print('    Delta_g=%5d  Lam=%9.3f  A log Lam + B + C1/Lam = %9.4f  =>  sqrt(t4^2+alpha_T^2) <= %.3e'
              % (Dg, Lv, rhs, mp.sqrt(57600 * max(rhs, 0)) / Lv ** 4 if rhs > 0 else float('nan')))
    print('  (for Delta_g <~ 10 the right side is not positive/meaningful: the O(1/Lam^2) terms are uncontrolled there)')
    # sensitivity to the unknown AdS coefficients of the |g3|^2 polynomials (t2 (ii)):
    print('\nSensitivity: k_AdS = k(0) + kappa2/Lam^2 + ..., kappa2 = -9677 b2/13860 + 225 f1/14 (b2, f1 not computed).')
    print('  relative change of the bound = -kappa2/(k(0) Lam^2):  e.g. |kappa2| = 100 at Lam = 100 -> %.2e' % (100 / float(k0) / 1e4))


if __name__ == '__main__':
    main()
