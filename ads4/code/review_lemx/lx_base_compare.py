"""lx_base_compare.py: could the reviewed base-multiplier lemma (code/norming/base_bounds.py, README §3) replace Lemma X on the
right edge of the [3.65, 3.75] contour?  Compare, for b = 0 and b = 4, log of rt_pp.PP0 / rt_pp.PP (Lemma X based, as used by
rt_cert.piece_right) with base_bounds.power_s1 / base_bound-style helicity on the same boxes (tau in [Lam-1, Lam_up],
d in [e^-3.75, e^-3.65], a >= 998), and with the exact maximum from lx_exact (b = 0: a = 998, n at L = 3.75)."""
import math, sys
from fractions import Fraction as F
sys.path.insert(0, '/home/user/MAJOR/ads4/code/norming'); sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic')
import base_bounds as BB, rt_pp
from lx_common import LAMf, S_hat_vec, n_at
from flint import arb
Lup = arb(998*999).sqrt().upper(); Lup_f = float(Lup)
tlo = F(LAMf) - 1; thi = F(Lup_f) + F(1, 10**9)
dlo = F(math.exp(-3.75)) * F(999999, 1000000); dhi = F(math.exp(-3.65)) * F(1000001, 1000000)
n = n_at(998.0, math.exp(-3.75))
print(' sigma box        | log Lemma X PP0(sig_hi) | log base_bounds.power_s1 | exact max log|S^| (a=998, L=3.75, tau=Lam)')
for (sl, sh) in [(F(0), F(1, 4)), (F(1, 4), F(1, 2)), (F(3, 4), F(1)), (F(9, 4), F(5, 2)), (F(13, 4), F(7, 2)), (F(17, 4), F(9, 2)), (F(21, 4), F(11, 2))]:
    lx = math.log(rt_pp.PP0(float(sh), LAMf - 1, Lup_f))
    bb = math.log(float(BB.power_s1(sl, sh, tlo, thi, dlo, dhi, F(998))))
    ex = max(S_hat_vec(n, 998.0, [complex(float(sl), LAMf), complex(float(sh), LAMf)]))
    print(f' [{float(sl):.2f}, {float(sh):.2f}]     | {lx:8.3f}               | {bb:8.3f}                 | {ex:.3f}')
