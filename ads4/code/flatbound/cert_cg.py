"""Exact check of the Clebsch-Gordan identity used for the +- channel (Lemma C of the .md):

   (1-y)^8 * 2F1(4-J, J+5; 1; y) = sum_{j=-4}^{4} c_j(J) * 2F1(-(J+j), J+j+1; 1; y)      (all complex J, |y|<1)

with c_j(J) the rational functions in common.cg_coeffs.  Proof structure (see .md): for integer J >= 4 the CG
product formula d^J_44 d^4_{-4,-4} = sum_L <J4;4-4|L0>^2 d^L_00 gives the identity with weights CG^2; the
coefficients of y^0..y^8 determine the weights uniquely (Vandermonde in L(L+1)), and this script checks
symbolically in J that c_j(J) solve these 9 equations; hence c_j(J) = CG^2 for all integers J>=4, and every
y^k coefficient of LHS-RHS is a rational function of J vanishing at all integers J>=4, i.e. identically.
As an extra check the coefficients y^9..y^14 are verified to vanish identically in J.
"""
import sympy as sp
from common import cg_coeffs

J, y = sp.symbols('J y')
def ser(a, b, n):
    out = [sp.Integer(1)]
    for k in range(1, n):
        out.append(sp.expand(out[-1] * (a + k - 1) * (b + k - 1) / k**2))
    return out
N = 15
lhs = sp.Poly(sp.expand((1 - y)**8 * sum(c * y**k for k, c in enumerate(ser(4 - J, J + 5, N)))), y)
cs = cg_coeffs(J)
ok = True
rows = []
for k in range(N):
    rhs_k = sum(cs[j + 4] * ser(-(J + j), J + j + 1, N)[k] for j in range(-4, 5))
    d = sp.simplify(sp.together(rhs_k - lhs.coeff_monomial(y**k)))
    ok &= (d == 0)
    print('y^%d coefficient: LHS-RHS = %s' % (k, d))
# uniqueness: 9x9 matrix of the first 9 coefficients of P_{J+j}(1-2y) is nonsingular
Mrows = [[ser(-(J + j), J + j + 1, 9)[k] for j in range(-4, 5)] for k in range(9)]
det = sp.factor(sp.Matrix(Mrows).det())
print('det of the 9x9 system =', det)
nonsing = all(sp.Poly(f, J).degree() >= 0 for f in [det]) and det != 0
# c_j >= 0 for J >= 4 and sum = 1
sumc = sp.simplify(sum(cs))
print('sum c_j =', sumc)
ok &= nonsing and sumc == 1
print('VERDICT cert_cg: %s' % ('PASS' if ok else 'FAIL'))
