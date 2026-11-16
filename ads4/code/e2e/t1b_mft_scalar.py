"""
Generalised-free (mean-field) decomposition for scalars with unequal dimensions, in the channel (12)(34) of
<phi_1 phi_2 phi_2 phi_1>:   sum_{n,J} P_{n,J} g^{a,a}_{D1+D2+2n+J, J}(z,zbar) = u^{(D1+D2)/2} v^{-D2},   a = b = (D2-D1)/2
(the correlator 1/(x_14^{2 D1} x_23^{2 D2}) divided by the standard prefactor K).  This fixes the conventions (a, b) and the
normalisation of the seed blocks  g ~ z^{(Delta-J)/2} zbar^{(Delta+J)/2}.
Mean-field coefficients (Fitzpatrick-Kaplan arXiv:1111.6972, eq. (2.32) form; h = d/2), in this block normalisation:
   P_{n,J} = (D1-h+1)_n (D2-h+1)_n (D1)_{J+n} (D2)_{J+n} / [ n! J! (J+h)_n (D1+D2+n-2h+1)_n (D1+D2+2n+J-1)_J (D1+D2+n+J-h)_n ]
(i) at D1 = D2 the sum over even J of 2 P_{n,J} is CHMRS (2.33) exactly (checked below), which fixes the block normalisation;
(ii) the decomposition is then tested at D1 != D2, where odd J contribute with the sign (+1) for the Wick contraction (14)(23)
     (and (-1)^J for (13)(24), i.e. for u^{(D1+D2)/2} in the ordering <phi_1 phi_2 phi_1 phi_2>).
"""
import sys, mpmath as mp, numpy as np
sys.path.insert(0, '../reggemoment')
from blocks import CBlock
mp.mp.dps = 30
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
h = mp.mpf(3)/2
def P_FK(n, J, D1, D2):
    rf = mp.rf
    return (rf(D1-h+1, n)*rf(D2-h+1, n)*rf(D1, J+n)*rf(D2, J+n)
            /(mp.factorial(n)*mp.factorial(J)*rf(J+h, n)*rf(D1+D2+n-2*h+1, n)*rf(D1+D2+2*n+J-1, J)*rf(D1+D2+n+J-h, n)))
def log_pMFT_CHMRS(Delta, J, Dphi, d=3):
    g = mp.loggamma; hh = mp.mpf(d)/2
    r = mp.log(2) + g(Delta-1) + g(J+hh) - 2*g(Dphi) - 2*g(Dphi-hh+1) - g(Delta-hh) - g(J+1)
    r += 2*g((J+Delta)/2) + 2*g((Delta-J-d+2)/2) - g(J+Delta-1) - g(Delta-J-d+1)
    r += g((Delta-J-2*d+2*Dphi+2)/2) + g((Delta+J-d+2*Dphi)/2) - g((Delta-J-2*Dphi+2)/2) - g((Delta+J+d-2*Dphi)/2)
    return r
pr("# (i) convention factor: CHMRS (2.33) [identical scalars, sum over even J only, factor 2] / P_FK at D1 = D2, n = 0, 1, 2")
Dp = mp.mpf('1.75')
for J in (0, 2, 4, 6):
    row = []
    for n in (0, 1, 2):
        pc = mp.exp(log_pMFT_CHMRS(2*Dp + 2*n + J, J, Dp))
        row.append(mp.nstr(pc/P_FK(n, J, Dp, Dp), 12))
    pr(f"   J = {J}: {row}")
kap = lambda J: 1
pr("# (ii) test at unequal dimensions:  sum_{n,J} P_{n,J} g^{a,a}  vs  u^{(D1+D2)/2} v^{-D2}")
for (D1, D2) in ((1.75, 1.75), (3.0, 1.75), (2.0, 2.75), (1.25, 2.5)):
    for (z, zb) in ((0.05, 0.08), (0.11, 0.03)):
        a = (D2 - D1)/2
        tot = mp.mpf(0)
        rho = z/(1+np.sqrt(1-z))**2; rhob = zb/(1+np.sqrt(1-zb))**2
        for n in range(0, 9):
            for J in range(0, 40):
                Delta = D1 + D2 + 2*n + J
                if (4*np.sqrt(rho*rhob))**Delta < 1e-22: continue
                g = CBlock(Delta, J, a, a, 40, prec='q').g(rho, rhob)[0]
                tot += kap(J)*P_FK(n, J, mp.mpf(D1), mp.mpf(D2))*g
        u = mp.mpf(z)*zb; v = (1-mp.mpf(z))*(1-mp.mpf(zb))
        rhs = u**((D1+D2)/mp.mpf(2))*v**(-mp.mpf(D2))
        pr(f"   D1={D1} D2={D2} z={z} zbar={zb}:  sum = {mp.nstr(tot, 15)}   exact = {mp.nstr(rhs, 15)}   rel.diff = {mp.nstr(abs(tot/rhs-1), 3)}")
open('out_t1b_mft_scalar.txt', 'w').write("\n".join(out) + "\n")
