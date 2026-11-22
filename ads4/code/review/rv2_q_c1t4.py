"""rv2_q_c1t4 (REVIEW_UNIFIED): leverage of a nu-independent light term 16 pi G c1 t4 in Result Q's delta=0.2 functional,
using Result D's t4 bound; and 2L_AdS re-check.  psi2 read as a polynomial in p (coefficient list, low order first)."""
import json, numpy as np
from scipy.integrate import quad
d = json.load(open('../resultq/functionals_q3.json'))['Q3-full-delta0.2']
c2 = np.array(d['psi2']); P = np.polynomial.Polynomial(c2)
I1 = quad(P, 0, 1, limit=200)[0]
for Dg, tb in ((300, 3.13e-6), (1000, 2.89e-8)):
    Lam = np.sqrt((Dg-2)*(Dg-1))
    L2 = 2*quad(lambda p: P(p)/(p*p + (1.5/Lam)**2), 0, 1, limit=400, points=[1e-3, 1e-2])[0]
    shift = 2*tb*Lam**2*I1
    print(f'Dg={Dg}: Lam={Lam:.2f} 2L_AdS={L2:.3f} (doc 72.04/92.81)  int psi2={I1:.4f}  |t4|max*Lam^2={tb*Lam**2:.4f}  '
          f'shift per unit c1 = {shift:.3f} = {100*shift/L2:.2f}% of 2L_AdS')
