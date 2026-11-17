"""s0: mode weights of the project kernels, checked against REGGE_MOMENT_TTTT.md Lemma L table (l = 2):
D1^2 ~ 4[(2,-2)+(-2,2)] + (2/9)(0,0) (SPINNING_E2E.md 3.1)  ->  M0 = N1/(64 dd52 dd92), Mflip = 1/(4 dd52 dd92);
A = D5 D1 ~ (2,-2) - (-2,2) (up to a constant)  ->  M0 ∝ g1 (a-1)(nu^2 + 2(a-3)^2 + 1/4)/dd72."""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e')); os.chdir(os.path.join(HERE, '..', 'e2e'))
import sympy as sp
from laplace import M_table, nu, a, g1
T = M_table(2)
def resp(w, h1, h3):
    return sp.simplify(sum(ww*(T[(q[0], q[1], h1, h3)][0] + g1*T[(q[0], q[1], h1, h3)][1]) for q, ww in w.items()))
W = {(2, -2): 4, (-2, 2): 4, (0, 0): sp.Rational(2, 9)}
dd = lambda b: nu**2 + (a - b)**2
N1 = 64*a**4 - 768*a**3 + 64*a**2*nu**2 + 3280*a**2 - 320*a*nu**2 - 5840*a + 16*nu**4 + 488*nu**2 + 3705
pe = sp.simplify((resp(W, 1, -1) + resp(W, -1, 1))/2); fe = sp.simplify((resp(W, 1, 1) + resp(W, -1, -1))/2)
print("D1^2 pres (sum over the two parity images /2) / [N1/(64 dd dd)] =", sp.simplify(pe/(N1/(64*dd(sp.Rational(5,2))*dd(sp.Rational(9,2))))))
print("D1^2 flip / [1/(4 dd dd)] =", sp.simplify(fe/(1/(4*dd(sp.Rational(5,2))*dd(sp.Rational(9,2))))))
WA = {(2, -2): 1, (-2, 2): -1}
po = sp.simplify((resp(WA, 1, -1) - resp(WA, -1, 1))/2)
print("A pres odd / [g1 (a-1)(nu^2+2(a-3)^2+1/4)/dd72] =", sp.simplify(po/(g1*(a-1)*(nu**2 + 2*(a-3)**2 + sp.Rational(1,4))/dd(sp.Rational(7,2)))))
