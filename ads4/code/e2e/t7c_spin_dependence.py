"""
Spin-dependent part of the O(1/m^2) correction (task 4): c2(spinning, helicity preserved / flipped) - c2(scalar skeleton with the same
dimensions (l+1, Dphi), the same k, nu, beta), compared with
  (i)  the exact normalisation polynomials of the conserved structures,   pol = -kappa_l [1/sinh^2(beta/2) + 1/cosh^2(beta/2)],
       kappa_1 = 1/4, kappa_2 = 5/2    [(J+1/2)^{2l} -> J(J+1) resp. (J-1)J(J+1)(J+2), and the same for Delta - 3/2];
  (ii) the helicity-dependent term of flat-space partial waves in its natural AdS form,   hel = -2 l^2 P'_{-1/2+i nu}(eta)/P_{-1/2+i nu}(eta)
       (-> 2 l^2 (nu/beta) J_1(nu beta)/J_0(nu beta) for beta -> 0;  REGGE_MOMENT_TTTT.md sect. 6.2 with lambda = l).
Usage: python3 t7c_spin_dependence.py     (needs data_t7_l*.pkl, data_t7b.pkl)
"""
import pickle, os, numpy as np, mpmath as mp
from formulas import Pnu
mp.mp.dps = 20
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
sb = pickle.load(open('data_t7b.pkl', 'rb'))
for l in (1, 2):
    if not os.path.exists(f'data_t7_l{l}.pkl'): continue
    s7 = pickle.load(open(f'data_t7_l{l}.pkl', 'rb'))
    kap = 0.25 if l == 1 else 2.5
    pr(f"# l = {l}:   D_pres = c2_pres - c2_scalar,  D_flip = c2_flip - c2_scalar")
    pr("#  beta   k   nu    c2_pres    c2_flip   c2_scalar |   D_pres   pol+hel  (pol, hel) |   D_flip")
    for (beta, k, nu), v in sorted(s7.items()):
        if k == 2: continue
        cs = sb[(l, beta, k, nu)]
        eta = mp.cosh(beta)
        P = Pnu(nu, eta); dP = mp.diff(lambda x: Pnu(nu, x), eta)
        pol = -kap*(1/np.sinh(beta/2)**2 + 1/np.cosh(beta/2)**2)
        hel = float(mp.re(-2*l*l*dP/P))
        pr(f"  {beta:4.1f}   {k}  {nu:4.1f}  {v[0]:9.2f}  {v[1]:9.2f}  {cs:9.2f}  | {v[0]-cs:8.2f}  {pol+hel:8.2f}  ({pol:7.2f}, {hel:7.2f}) | {v[1]-cs:8.2f}")
open('out_t7c_spin_dependence.txt', 'w').write("\n".join(out) + "\n")
