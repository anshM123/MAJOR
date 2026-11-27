"""validate hev.py against the reviewer's exact-rational closed form (rv_closed.py, imported read-only)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'review_uniform_gap'))
from fractions import Fraction as Fr
from rv_closed import closed_form
from hev import Coeffs, S_state
lam2 = 997002
cf = Coeffs(lam2, 80)
for (J, b, a) in [(10, 0, 998), (11, 4, 998), (30, 0, 1500), (40, 4, 998), (64, 0, 998), (64, 4, 1203)]:
    Delta = Fr(a)+J+2; T = (Delta-1)*(Delta-2); m2 = T-J*(J+1)
    ex = closed_form(J, b, m2, lam2)
    S, w, info = S_state(cf, J, b, a)
    print(J, b, a, float(ex), S.str(25), float((S - (ex.numerator/__import__('flint').arb(ex.denominator))).mid())/float(ex))
