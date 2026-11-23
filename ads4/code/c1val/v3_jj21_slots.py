#!/usr/bin/env python3
"""
c1val / v3: consistency of the ASD boundary relation with the published 3d spinor-helicity <TT>, <TTT>
(Jain, John, Mehta, Nizami, Suresh, arXiv:2106.00016, eqs. JsJstwopointfn, evenTTTso, tttoddsphresult, TTTshal and the
display after eq. tttsph1; LaTeX read from source).

Prediction of the ASD-closure analysis (C1_VALUE.md sec. 3.3): with two ASD legs of label L, the Witten (Einstein) 3-pt
function equals kappa * (gravitational Chern-Simons contact) EXACTLY in the all-L slot, and differs from it in the mixed
slot (L, L, -L), the difference being the self-dual content Phi_34 that survives in AdS.
In JJ21 the CS contact is the c_T' (parity-odd, non-homogeneous) structure.  The prediction then reads:
  (a) <T-T-> and <T-T-T->_nh depend on (c_T, c_T') only through the SAME combination c_T + i c_T'
      (so the shift c_T' -> c_T' + i c_T, i.e. subtracting kappa*CS, removes both);
  (b) in the mixed slot <T-T-T+> the even and odd form factors differ (X != Y), so the same shift leaves c_T (X - Y) != 0.
Writes cache/v3_jj21_slots.txt.
"""
import os, sys, random
import sympy as sp

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cache', 'v3_jj21_slots.txt')
lines = []
def log(s=''):
    print(s); lines.append(str(s))

k1, k2, k3 = sp.symbols('k1 k2 k3', positive=True)
cT, cTp = sp.symbols('c_T c_Tp')
E = k1 + k2 + k3; b = k1*k2 + k2*k3 + k3*k1; c = k1*k2*k3; a12 = k1 + k2
F = (E**3 - E*b - c)/c**2                                     # nh form factor of <T-T-T-> (both parities)
twopt_mm = cT + sp.I*cTp                                       # <T-T-> coefficient (eq. JsJstwopointfn)
ttt_mmm_nh = cT*F + sp.I*cTp*F                                 # eq. evenTTTso (c_T part) + tttoddsphresult (x c123)
X = (E - 2*k3)**2*(E**3 - E*b - c)/(E**2*c**2)                 # even nh <T-T-T+> (eq. evenTTTso)
Y = ((E - 2*k3)**3 - (E - 2*k3)*(b - 2*k3*a12) + c)/c**2       # odd nh <T-T-T+> (display after eq. tttsph1; the
                                                                # TTTshal line has a dimensionally inconsistent typo)
ttt_mmp_nh = cT*X + sp.I*cTp*Y

shift = {cTp: cTp + sp.I*cT}   # subtract kappa*CS with the ASD normalisation: c_T' -> c_T' + i c_T
log('v3: JJ21 <TT>, <TTT> non-homogeneous helicity components under subtraction of the ASD-normalised CS contact')
log(f'(a) <T-T-> coeff after shift: {sp.simplify(twopt_mm.subs(shift) - sp.I*cTp)}  (residual beyond the odd piece; 0 = removed)')
log(f'    <T-T-T->_nh after shift : {sp.simplify(ttt_mmm_nh.subs(shift) - sp.I*cTp*F)}')
mixed = sp.simplify(ttt_mmp_nh.subs(shift) - sp.I*cTp*Y)
log(f'(b) <T-T-T+>_nh after shift: c_T*(X - Y) with X - Y = {sp.factor(sp.simplify(X - Y))}')
random.seed(1)
vals = []
for _ in range(5):
    kv = sorted([random.uniform(0.2, 2.0) for _ in range(3)])
    # triangle inequality for |k_i| of a closed momentum triangle
    if kv[2] >= kv[0] + kv[1]:
        kv[2] = 0.9*(kv[0] + kv[1])
    sub = {k1: kv[0], k2: kv[1], k3: kv[2]}
    vals.append((kv, float((X - Y).subs(sub)), float(X.subs(sub)), float(Y.subs(sub))))
for kv, d, xv, yv in vals:
    log(f'    k = ({kv[0]:.3f},{kv[1]:.3f},{kv[2]:.3f}):  X = {xv:+.5e}  Y = {yv:+.5e}  X-Y = {d:+.5e}')
# total-energy pole: X has a double pole at E = 0 (flat-space amplitude), Y has none
Es = sp.symbols('Es')
Xe = sp.simplify(X.subs(k3, Es - k1 - k2)); Ye = sp.simplify(Y.subs(k3, Es - k1 - k2))
lx = sp.limit(Xe*Es**2, Es, 0); ly = sp.limit(Ye*Es**2, Es, 0)
log(f'    E^2 X at E->0: {sp.simplify(lx)}   E^2 Y at E->0: {sp.simplify(ly)}  (X carries the flat-space E^-2 pole; Y does not)')
ok = (sp.simplify(twopt_mm.subs(shift) - sp.I*cTp) == 0 and sp.simplify(ttt_mmm_nh.subs(shift) - sp.I*cTp*F) == 0
      and all(abs(v[1]) > 1e-8 for v in vals) and sp.simplify(ly) == 0 and sp.simplify(lx) != 0)
log(f'ALL CHECKS PASS: {ok}')
with open(OUT, 'w') as fh:
    fh.write('\n'.join(lines) + '\n')
sys.exit(0 if ok else 1)
