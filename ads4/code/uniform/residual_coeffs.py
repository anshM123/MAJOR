"""
Curvature-subtracted NLO/NNLO remainder coefficients (Df = 3, k = 2..5), from t1's R1_j, R2_j:
  exact fixed-J series  m^{2k}C = BP(x,J) + r1(x,J)/m^2 + r2(x,J)/m^4 + ...,  r1 = sum_j R1_j x^{2j-2},  r2 = sum_j R2_j x^{2j-4};
  Regge-curvature part (exact): 2 P_{-1/2+i x m}(1 + e/m^2) - 2 J_0((2J+1)x),  e = 2JJ + 1/2,  whose 1/m^2, 1/m^4 coefficients are
  2 (-e/2)^n/(n!)^2 e_1(n) x^{2n-2},  2 (-e/2)^n/(n!)^2 e_2(n) x^{2n-4},  e_r = elementary symmetric polynomials of {(i+1/2)^2, i < n}.
  r_ires = r_i - curvature part: the remainder relative to the hybrid (flat + leading-Regge curvature).
Output: cache/residual_coeffs_Df3.pkl (coefficients as polynomials in JJ = J(J+1)), printed table.
"""
import os, sys, pickle
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sympy as sp
from chmrs_series import k, Df, J2
d = pickle.load(open(os.path.join(HERE, 'cache', 'remainder_coeffs.pkl'), 'rb'))
R1 = [sp.sympify(eval(s, sp.__dict__.copy())) for s in d['R1']]
R2 = [sp.sympify(eval(s, sp.__dict__.copy())) for s in d['R2']]
e = 2*J2 + sp.Rational(1, 2)
esym = lambda vals, r: sum(sp.prod(c) for c in combinations(vals, r))
out = {}
for kk in (2, 3, 4, 5):
    r1, r2 = [], []
    for j in range(1, 5):
        vals = [sp.Rational(2*i + 1, 2)**2 for i in range(j)]
        r1.append(sp.expand((R1[j - 1] - 2*(-e/2)**j/sp.factorial(j)**2*esym(vals, 1)).subs({k: kk, Df: 3})))
    for j in range(2, 5):
        vals = [sp.Rational(2*i + 1, 2)**2 for i in range(j)]
        r2.append(sp.expand((R2[j - 2] - 2*(-e/2)**j/sp.factorial(j)**2*esym(vals, 2)).subs({k: kk, Df: 3})))
    out[kk] = {'r1': [[float(v) for v in sp.Poly(c, J2).all_coeffs()[::-1]] for c in r1],
               'r2': [[float(v) for v in sp.Poly(c, J2).all_coeffs()[::-1]] for c in r2]}
    print(kk, 'r1 (coefficients of JJ^0, JJ^1, ... for x^0, x^2, x^4, x^6):', [[round(v, 3) for v in c] for c in out[kk]['r1']])
    print('   r2 (x^0, x^2, x^4):', [[round(v, 1) for v in c] for c in out[kk]['r2']])
pickle.dump(out, open(os.path.join(HERE, 'cache', 'residual_coeffs_Df3.pkl'), 'wb'))
