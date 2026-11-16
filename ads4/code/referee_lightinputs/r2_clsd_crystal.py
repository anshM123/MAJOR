"""Referee r2: CLS-D (arXiv:2311.04271) subtraction factors, independent re-derivation of the crystal values.

From the f^+ lines alone one can solve for the crystal values of u', v':  u' = m^4/nu^4, v' = (m^2-nu^2)^2/nu^4.
Check: all four cases of (106)/(107) [arXiv v-numbering (3.75)/(3.76)] reproduce (108) [(3.77)] with these values,
and read off the Regge-spin shift n (f ~ m^{-2n}) for each parity.  Also evaluate f^-_{3;2}.
"""
import sympy as sp
m, nu = sp.symbols('m nu', positive=True)
up, vp = m**4/nu**4, (m**2-nu**2)**2/nu**4
def fp(k, J):
    if (k-2*J-2) % 4 == 0: return (vp-up)/(up*vp)**sp.Rational(k-2*J+2, 4)
    return ((vp-up)**3-(vp-up))/(4*(up*vp)**sp.Rational(k-2*J+4, 4))
def fm(k, J):
    if (k-2*J-3) % 4 == 0: return -1/(up*vp)**sp.Rational(k-2*J+1, 4)
    return -((up-vp)**2-1)/(4*(up*vp)**sp.Rational(k-2*J+3, 4))
crys_p = lambda k, J: -nu**(2*(k-2*J+1))*(2*m**2-nu**2)/(m**2*(m**2-nu**2))/(m**2*(m**2-nu**2))**sp.Rational(k-2*J, 2)
crys_m = lambda k, J: -nu**(2*(k-2*J))*nu**2/(m**2*(m**2-nu**2))/(m**2*(m**2-nu**2))**sp.Rational(k-2*J-1, 2)
def npow(e):
    t = sp.symbols('t', positive=True)
    return -sp.limit(sp.log(sp.Abs(sp.simplify(e.subs(m, 1/t))))/sp.log(t), t, 0)/(-2)  # f ~ m^{-2n}
print(" J  k  par  matches(108)  n=Regge shift  base 2J-1  => J0 <")
for J in [1, 2, 3]:
    for k in range(2, 10):
        f, c = (fp(k, J), crys_p(k, J)) if k % 2 == 0 else (fm(k, J), crys_m(k, J))
        same = sp.simplify(sp.factor(f - c)) == 0
        n = npow(c)
        print(f" {J}  {k}  {'+' if k%2==0 else '-'}    {same}          {n}             {2*J-1}        {2*J-1+n}")
print("f^-_{3;2} from (107):", sp.simplify(fm(3, 2)), "  (CLS-D (111)/(3.80) prints +1)")
