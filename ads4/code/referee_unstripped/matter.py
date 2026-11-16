"""Referee: matter low sides of B'_2, B'_3 from (B.1) via residues at s=0,q and the light poles."""
import sympy as sp
s, q, e = sp.symbols('s q e'); a0, a2, ml = sp.symbols('a0 a2 m_l', positive=True)
def fm(S, U):
    T = -S - U
    return a0/(ml**2 - T) + a2*(1 - 6*S*U/ml**4)/(ml**2 - T)
def K(k): return (2*s - q)*(s*(s - q))**sp.Rational(-(k + 2), 2) if k % 2 == 0 else (s*(s - q))**sp.Rational(-(k + 1), 2)
def low(k):
    A = (q - s)**4*fm(s, -q); B = s**4*fm(q - s, -q)
    ex = sp.together(K(k)*((A + B) if k % 2 == 0 else (-A + B)))
    tot = 0
    for a in [0, q, q - ml**2, ml**2]:
        tot += sp.series(ex.subs(s, a + e), e, 0, 1).removeO().coeff(e, -1)
    return sp.factor(sp.simplify(tot/2))
c2 = (-2*a0*ml**6 - 3*a0*ml**4*q - 2*a2*ml**6 + 9*a2*ml**4*q + 6*a2*ml**2*q**2 + 6*a2*q**3)/ml**4
c3 = (-(a0 + a2)*ml**4 + 6*a2*ml**2*q + 6*a2*q**2)/ml**4
L2, L3 = low(2), low(3)
print(L2); print(L3)
print('k=2', 'PASS' if sp.simplify(L2 - c2) == 0 else 'FAIL', ' k=3', 'PASS' if sp.simplify(L3 - c3) == 0 else 'FAIL')
