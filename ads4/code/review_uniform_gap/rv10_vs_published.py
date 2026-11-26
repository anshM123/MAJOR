"""rv10_vs_published.py -- context: new-family ideal action S = H/w^3 vs the PUBLISHED Delta_g=1000 functional's H/w^3
(certified moments code/fullpos/cache/moments_1000_K1000.json via xl_common.load_cert_moments, read-only), twist edge."""
import sys, os
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'xcheck_ladder')))
import xl_common as xl
from flint import arb, fmpq, ctx
from fractions import Fraction as Fr
from rv_closed import closed_form
ctx.prec = 2000
Z, Y = xl.load_cert_moments(1000, 1000, kmax=70)
Q = xl.Q_exact(70, 997002)
print("J b | S_new (ideal, uncut) | H_pub/w^3 | (S_new - H_pub/w^3)/S_new")
for J, b in [(4, 4), (5, 4), (6, 0), (12, 0), (12, 4), (20, 0), (20, 4), (32, 0), (32, 4), (64, 0), (64, 4)]:
    Dl = J+1000; T = Fr((Dl-1)*(Dl-2)); m2 = T-J*(J+1)
    Ta = arb(fmpq(T.numerator, T.denominator))
    H = xl.H_ang(J, b, Ta, Z, Y, Q, 997002); w3 = (1-arb(J*(J+1))/Ta)**3
    S = closed_form(J, b, m2, 997002); Sa = arb(fmpq(S.numerator, S.denominator))
    print(J, b, '|', Sa.str(10), '|', (H/w3).str(10), '|', ((Sa-H/w3)/Sa).str(3))
