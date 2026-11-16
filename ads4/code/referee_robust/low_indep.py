"""Referee: independent exact low side for the robust functionals (coefficients typed from ROBUST_FUNCTIONAL.md section 3,
not read from the author's cache).  -B_k|low computed by Laurent series (not sympy.residue), then the functional applied."""
import sympy as sp
from fractions import Fraction as Fr

s, p, G, g3, g4, g5, g6, g6p, ml, a0, a2, e = sp.symbols('s p G g3 g4 g5 g6 g6p m_l a0 a2 e')

def flow(S, U):
    T = -S - U
    return 8*sp.pi*G/(S*T*U) + 2*sp.pi*G*S*U/T*g3 + g4 + g5*T + g6*T**2 - g6p*S*U

def fmat(S, U):
    T = -S - U
    return a0/(ml**2 - T) + a2*(1 - 6*S*U/ml**4)/(ml**2 - T)

def res(expr, z):
    # residue via Laurent coefficient: substitute s = z + e, series in e, coefficient of 1/e
    ser = sp.series(sp.together(expr.subs(s, z + e)), e, 0, 1).removeO()
    return sp.simplify(ser.coeff(e, -1))

def mBlow(k, f, poles):
    ev = f(s, -p**2) + f(p**2 - s, -p**2); od = -f(s, -p**2) + f(p**2 - s, -p**2)
    if k % 2 == 0:
        integrand = (2*s - p**2) / (s*(s - p**2))**sp.Rational(k - 2, 2) * ev
    else:
        integrand = od / (s*(s - p**2))**sp.Rational(k - 3, 2)
    return sp.expand(sp.simplify(sum(res(integrand, z) for z in poles) / 2))

def omp3(c):  # (1-p)^3 * sum c_i p^i
    return sp.expand((1 - p)**3 * sum(sp.Rational(str(ci)) * p**i for i, ci in enumerate(c)))

FUNS = {
  'Fg3': ([0, 65, 155, 47, 0, 0], [162, -3854, 9202, 0, 0], 5),   # psi2 = (1-p)^3 p (65+155p+47p^2), psi3 = 2(1-p)^3(81-1927p+4601p^2)
  'M02': ([0, Fr(711, 8), Fr(129649, 1000), Fr(561721, 1000), Fr(-1013557, 1000), Fr(997013, 1000)],
          [Fr(54971, 1000), Fr(-102583, 100), Fr(-31903, 25), Fr(946687, 1000), Fr(7534331, 500)], Fr(1129, 250)),
  'M01': ([0, Fr(64493, 1000), Fr(3779, 40), Fr(208699, 500), Fr(-790783, 1000), Fr(884501, 1000)],
          [Fr(1001, 50), Fr(-74971, 250), Fr(-3284623, 1000), Fr(-2054369, 1000), Fr(2261339, 100)], Fr(4259, 1000)),
}

if __name__ == '__main__':
    low = {k: mBlow(k, flow, [0, p**2]) for k in (2, 3, 4, 5)}
    mat = {k: mBlow(k, fmat, [0, p**2, ml**2, p**2 - ml**2]) for k in (2, 3, 4, 5)}
    for k in (2, 3, 4, 5): print('-B%d|low =' % k, low[k], ' | matter:', sp.factor(mat[k]))
    xx = sp.Symbol('xx', positive=True)
    for name, (x, y, a) in FUNS.items():
        a = sp.Rational(str(a))
        psi2 = omp3(x); psi3 = omp3(y)
        if name == 'Fg3': psi3 = omp3([162, -3854, 9202])
        def act(ex):
            e2 = sp.expand(ex[2]); pole = e2.coeff(p, -2); reg = sp.expand(e2 - pole/p**2)
            t = sp.integrate(sp.expand(psi2*reg), (p, 0, 1)) + sp.integrate(sp.expand(psi3*ex[3]), (p, 0, 1))
            t += -a*sp.diff(sp.expand(ex[4]).subs(p, sp.sqrt(xx)), xx).subs(xx, 0) + a*sp.expand(ex[5]).subs(p, 0)
            return sp.expand(t), pole
        tot, pole = act(low)
        k = sp.nsimplify(-tot.coeff(g3)/(sp.pi*G))
        x0 = psi2.coeff(p, 1)
        Lfin = sp.limit(sp.integrate(sp.expand(psi2/p**2), (p, e, 1)) - x0*sp.log(1/e), e, 0)
        A = 16*x0/k; B = 16*Lfin/k
        print('\n==', name, ': coeffs g4,g5,g6,g6p =', [tot.coeff(c) for c in (g4, g5, g6, g6p)], ' pole', pole)
        print('   k =', k, float(k), ' A =', float(A), ' B =', float(B), ' B_AdS =', float(B - A*sp.log(sp.Rational(3, 2))))
        tm, _ = act(mat)
        C0 = sp.factor(tm.coeff(a0)); C2 = sp.factor(tm.coeff(a2))
        print('   C0 =', C0); print('   C2 =', C2)
        for nm, cc in (('C0', C0), ('C2', C2)):
            vals = [float(cc.subs(ml, sp.Rational(i, 200))) for i in range(1, 201)]
            num = sp.Poly(sp.numer(sp.together(cc)), ml)
            print('   %s: max on grid (0,1] = %.4g ; real roots of numerator in [0,1]: %d' % (nm, max(vals), sp.count_roots(num, 0, 1)))
