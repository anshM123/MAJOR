"""u2: curvature-subtracted remainder coefficients of the exact scalar skeleton (CHMRS order-1/m^8 fixed-J series) for k = 2..7.

Same construction as uniform/residual_coeffs.py (which covers k = 2..5): from the symbolic NLO/NNLO coefficients R1_j, R2_j of
uniform/t1_structure.py (functions of k, Df, J2; read from uniform/cache/remainder_coeffs.pkl, not modified) subtract the exact
leading-Regge-curvature part  2 P_{-1/2+i x m}(1 + e/m^2) - 2 J_0((2J+1)x),  e = 2JJ + 1/2.  Odd k: CHMRS's kernel continued in k [I].
Output: cache/rcoef_Df<Df>.pkl  ({k: {'r1': [[coeffs in JJ] for x^0,x^2,x^4,x^6], 'r2': [... x^0,x^2,x^4]}}), and a table (cache/u2.log).
"""
import os, sys, pickle
from itertools import combinations
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'uniform'))
import sympy as sp
from chmrs_series import k, Df, J2          # symbols only
d = pickle.load(open(os.path.join(HERE, '..', 'uniform', 'cache', 'remainder_coeffs.pkl'), 'rb'))
R1 = [sp.sympify(eval(s, sp.__dict__.copy())) for s in d['R1']]
R2 = [sp.sympify(eval(s, sp.__dict__.copy())) for s in d['R2']]
e = 2*J2 + sp.Rational(1, 2)
esym = lambda vals, r: sum(sp.prod(c) for c in combinations(vals, r))
log = []
for DFV in (sp.Integer(3), sp.Integer(2), sp.Rational(5, 2)):
    out = {}
    for kk in range(2, 8):
        r1, r2 = [], []
        for j in range(1, 5):
            vals = [sp.Rational(2*i + 1, 2)**2 for i in range(j)]
            r1.append(sp.expand((R1[j - 1] - 2*(-e/2)**j/sp.factorial(j)**2*esym(vals, 1)).subs({k: kk, Df: DFV})))
        for j in range(2, 5):
            vals = [sp.Rational(2*i + 1, 2)**2 for i in range(j)]
            r2.append(sp.expand((R2[j - 2] - 2*(-e/2)**j/sp.factorial(j)**2*esym(vals, 2)).subs({k: kk, Df: DFV})))
        out[kk] = {'r1': [[float(v) for v in sp.Poly(c, J2).all_coeffs()[::-1]] for c in r1],
                   'r2': [[float(v) for v in sp.Poly(c, J2).all_coeffs()[::-1]] for c in r2]}
        s = (f"Df={float(DFV):g} k={kk}: r1res(x=0,J=0) = {float(r1[0].subs(J2, 0)):+9.3f}  r1res x^2..x^6 (J=0) = "
             f"{[round(float(c.subs(J2, 0)), 1) for c in r1[1:]]}   r2res x^0..x^4 (J=0) = {[round(float(c.subs(J2, 0)), 0) for c in r2]}")
        print(s, flush=True); log.append(s)
    pickle.dump(out, open(os.path.join(HERE, 'cache', 'rcoef_Df%g.pkl' % float(DFV)), 'wb'))
open(os.path.join(HERE, 'cache', 'u2.log'), 'w').write("\n".join(log) + "\n")
