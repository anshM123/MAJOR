r"""no_eft/nl1b_sanity.py [DV/NE]: consistency checks of the light-state evaluation of nl1_light_action.py.
 (i) heavy states: hhat/T*(1-r)^3 equals fp_poly.Action.H_ang (the certified FP evaluator) -- same formula, different bookkeeping;
 (ii) J = 0, 2 light states (the only light single-trace spins allowed by A1'): ladder of length <= 3, moderate values;
 (iii) continuity in twist at fixed J: hhat at twist 2, 10, 100, 500, 900, 990, 998, 1000, 1010.
Usage: python3 nl1b_sanity.py   (K = 300 moments; Delta_g = 1000)."""
import sys, os, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'fullpos'))
import fp_poly as FP
from flint import arb
import nl1_light_action as N
Z, Y = FP.load_moments(1000, 300); lam = 998*999
A = FP.Action(1000, Z, Y, FP.Q_exact(300, lam))
def lg(x):
    return float((abs(x)).log()/arb(10).log()) if not (0 in x) else float('nan')
print('(i) heavy states: H_ang (fp_poly) vs hhat*(1-r)^3/T')
for J, tw, b in ((10, 1000, 0), (64, 1000, 4), (200, 1500, 0), (290, 1000, 4)):
    h, T, r = N.hhat(A, J, b, J + tw)
    Ha = A.H_ang(J, b, J + tw and (J + tw - 1)*(J + tw - 2))
    print(f'  J={J} tw={tw} b={b}: H_ang={float(Ha.mid()):.12e}  hhat(1-r)^3/T={float((h*(1-r)**3/T).mid()):.12e}')
print('(ii) light J = 0, 2 (b = 0): log10|hhat| at twist 1(J=2 only: T itself), 2, 4, 10, 100, 500')
for J in (0, 2):
    for tw in ((1, 2, 4, 10, 100, 500) if J == 2 else (2, 4, 10, 100, 500)):
        h, T, r = N.hhat(A, J, 0, J + tw)
        print(f'  J={J} tw={tw}: hhat = {float(h.mid()): .6e}')
print('(iii) continuity in twist, J = 64, b = 0: log10|hhat|, sign')
for tw in (2, 10, 100, 300, 500, 700, 900, 990, 998, 1000, 1010, 1100):
    h, T, r = N.hhat(A, 64, 0, 64 + tw)
    print(f'  tw={tw:5d}: log10|hhat| = {lg(h):9.3f}  sign {"+" if h > 0 else "-"}')
