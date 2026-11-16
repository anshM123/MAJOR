"""
Exact helicity basis of <J J O> in the differential basis, from Caron-Huot--Li's parity-even helicity structures (eq. VVstcor of
arXiv:2102.08160, without its overall factor 2^{(J+1)/2} n), identified at exact rational (Delta, J) by chl_map.identify:
     T^{E,h} = K (a_h v_1 + b_h v_2)   (K common to both h).
Closed forms (fitted at 6 exact points with a 6-parameter polynomial ansatz, then CHECKED exactly at further points, including
Weyl-reflected labels (3-Delta, -1-J)):
     b_opp/a_opp  = J(J+1)(Delta^2 - 3 Delta - J^2 - J + 4)        = (A_-^2 - 1/4)(m^2 + 2)
     b_same/a_same = -(Delta-1)(Delta-2)(Delta^2 - 3 Delta - J^2 - J) = -(A_+^2 - 1/4)(m^2 - 2)
     a_opp = -2/(J(J-1) alpha Y),   a_same = -2/(Delta(Delta-1) alpha Y),   Y = Delta^2 - 3Delta + J^2 + J,  alpha = Delta - J - 2.
"""
import os, pickle, sys
import sympy as sp
from chl_map import frame_Q, identify
from embed2 import CACHE
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def closed(Dv, Jv):
    Y = Dv**2 - 3*Dv + Jv**2 + Jv; al = Dv - Jv - 2
    return dict(opp=(-sp.Integer(2)/(Jv*(Jv - 1)*al*Y), Jv*(Jv + 1)*(Dv**2 - 3*Dv - Jv**2 - Jv + 4)),
                same=(-sp.Integer(2)/(Dv*(Dv - 1)*al*Y), -(Dv - 1)*(Dv - 2)*(Dv**2 - 3*Dv - Jv**2 - Jv)))

if __name__ == "__main__":
    ws, vecs = pickle.load(open(os.path.join(CACHE, 'consbasis_l1.pkl'), 'rb'))
    Qs, BV, syms = frame_Q()
    R = sp.Rational
    pr("# CHL helicity structures in the (v_1, v_2) basis: exact identification vs closed forms (fourpt/t1b_helicity_l1.py)")
    pts = [(R(97, 10), 4), (R(23, 2), 3), (R(37, 3), 9), (R(19, 2), 2), (3 - R(97, 10), -5), (3 - R(37, 3), -10), (R(61, 4), 11)]
    for (Dv, Jv) in pts:
        try:
            r = identify(Dv, Jv, Qs, BV, syms, ws, vecs)
        except AssertionError:
            pr(f"  ({Dv},{Jv}): degenerate random kinematics, skipped"); continue
        c = closed(Dv, Jv)
        ok = all(sp.simplify(r[h][0] - c[h][0]) == 0 and sp.simplify(r[h][1]/r[h][0] - c[h][1]) == 0 for h in ('opp', 'same'))
        pr(f"  (Delta,J) = ({Dv},{Jv}): a, b/a  opp = {r['opp'][0]}, {r['opp'][1]/r['opp'][0]};  same = {r['same'][0]}, {r['same'][1]/r['same'][0]};  closed forms exact: {ok}")
    open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'out_t1b_helicity_l1.txt'), 'w').write("\n".join(out) + "\n")
