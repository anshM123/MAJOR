"""
Calibration for t5_helicity_l2: CHL's momentum-space T^{E,h} (chl_momentum) versus the exact Weyl-covariant helicity families of currents
(hhat_h of SPINNING_4PT.md sect. 2.3): N_h with T^{E,h} = N_h hhat_h at (Delta,J) and at (3-Delta,-1-J); kappa_h = N_h(Delta)/N_h(Delta~);
kappa_h A_h with A_opp = (-J)_2/(J+1)_2, A_same = (Delta-1)_2/(2-Delta)_2 (CHL's MFT helicity factors).  For currents the MFT normalisation
is known exactly (SPINNING_4PT.md sect. 4.6) and Weyl-even: this run fixes the reference value of kappa_h A_h.
"""
import os, sys
import sympy as sp
from chl_momentum import chl_components, our_components
from t1_mft_l1 import hopp, hsame
from infframe import Delta, J
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
R = sp.Rational
for (Dv, Jv) in ((R(97, 10), R(13, 3)), (R(53, 4), R(11, 2))):
    N = {}
    for (D_, J_) in ((Dv, Jv), (3 - Dv, -1 - Jv)):
        C = chl_components(1, D_, J_)
        for name, e, idx in (('same', hsame, 0), ('opp', hopp, 1)):
            O = our_components(1, {w: c.subs({Delta: D_, J: J_}) for w, c in e.items()}, D_, J_)
            r = [C[h][idx]/O[h] for h in O if abs(O[h]) > 1e-25*max(abs(x) for x in O.values())]
            wrong = max(abs(C[h][1 - idx]) for h in O if abs(O[h]) > 1e-25*max(abs(x) for x in O.values()))
            pr(f"  ({D_},{J_}) {name}: N_h on the two nonzero components {['%.12g' % z.real for z in r]}; CHL's other helicity there: {wrong:.1e}")
            N[(name, D_)] = r[0]
    for name in ('opp', 'same'):
        A = (sp.rf(-Jv, 2)/sp.rf(Jv + 1, 2)) if name == 'opp' else (sp.rf(Dv - 1, 2)/sp.rf(2 - Dv, 2))
        k = N[(name, Dv)]/N[(name, 3 - Dv)]
        pr(f"  ({Dv},{Jv}) {name}: kappa_h = {k.real:.12g},  kappa_h A_h = {k.real*float(A):.12g}")
open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t6_kappa_l1.txt'), 'w').write("\n".join(out) + "\n")
