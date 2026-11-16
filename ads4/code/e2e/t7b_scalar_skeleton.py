"""
First correction of the scalar "skeleton" <phi_a phi phi_a phi> with dim(phi_a) = l + 1 (the dimension of the conserved operator),
dim(phi) = Dphi, at the same (Delta, J) as the spinning runs:  R = exact/leading,  c2 = lim m^2 (R - 1).
Used as a reference for the spin dependence of the O(1/m^2) terms (t7_corrections.py).   Usage: python3 t7b_scalar_skeleton.py
"""
import sys, pickle, numpy as np, mpmath as mp
from moments import ScalarMoments
from t3_scalar_mixed import prediction
from t4_spinning_e2e import main_cases, DPHI
from cases import info
mp.mp.dps = 25
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
summary = {}
for l in (1, 2):
    D1, D2 = float(l + 1), DPHI
    data = {}
    for (Delta, J) in main_cases():
        m, b_, eta = info(Delta, J)
        mom = ScalarMoments(Delta, J, D1, D2, int(min(14000, max(800, 60*Delta))))
        for k in (2, 4, 6):
            for nu in (0.0, 1.0, 3.0):
                sg, lg = mom.moment(k, nu)
                R = sg*mp.exp(mp.mpf(lg))/prediction(Delta, J, D1, D2, k, nu).real
                data.setdefault((round(float(b_), 1), k, nu), []).append((m, R))
    pr(f"# scalar skeleton, dimensions ({D1}, {D2}):  m^2 (R-1) for the available m, and c2 from the two largest m (ansatz c2/m^2 + c4/m^4)")
    pr("#  beta   k   nu    m^2(R-1) ...                                                  c2")
    for key in sorted(data):
        seq = sorted(data[key], key=lambda t: t[0])
        (m1, r1), (m2, r2) = seq[-2], seq[-1]
        a11, a12, a21, a22 = 1/m1**2, 1/m1**4, 1/m2**2, 1/m2**4
        c2 = ((r1-1)*a22 - (r2-1)*a12)/(a11*a22 - a12*a21)
        summary[(l,) + key] = float(c2)
        pr(f"  {key[0]:4.1f}   {key[1]}  {key[2]:4.1f}   " + " ".join(f"{float(m*m*(R-1)):9.3f}" for m, R in seq) + f"    {float(c2):9.3f}")
open('out_t7b_scalar_skeleton.txt', 'w').write("\n".join(out) + "\n")
pickle.dump(summary, open('data_t7b.pkl', 'wb'))
