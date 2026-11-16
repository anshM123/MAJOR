"""
<TTT> in d = 3: O = T (J = 2, Delta = 3, s = 5/2).  Conserved structures: 2 parity-even + 1 parity-odd.
Computes (ahat^T_{2,2}, ahat^T_{-2,2}) for the free-boson, free-fermion and parity-odd structures.
Result:  ahat_{2,2} = (n_B + n_F)/(n_B + n_F) = 1,   c = ahat_{-2,2} = (n_B - n_F + i kappa p_odd)/(n_B + n_F).
"""
from fractions import Fraction as Fr
import sympy as sp, mpmath as mp
import structures as st, wick_free as wf, pipeline as pl, kernels_lib as kl, free_checks as fc
mp.mp.dps = 30
J, Delta, sig = 2, 3, Fr(5, 2)
NT = kl.norm_T_plus()
print("parity-even structures at s = 5/2:")
mon_e, bas_e = kl.conserved(2, 'even', sval=sig)
print("   number:", len(bas_e))
_, fS = fc.fit('scalar', 2, 3, wf.c_scalar_T())
_, fF = fc.fit('fermion', 2, 3, wf.c_fermion_T(), 'z')
lamS, lamF = fS['even'][1], fF['even'][1]
print("   free scalar  <TTT> (units (4 pi)^-3):", dict((m, c) for m, c in zip(mon_e, lamS) if c != 0))
print("   free fermion <TTT> (units c_F^3 kappa^3):", dict((m, c) for m, c in zip(mon_e, lamF) if c != 0))
print("   independent:", sp.Matrix([lamS, lamF]).rank() == 2)
for nm, lam, cal in (('boson', lamS, (4*mp.pi)**(-3)/((mp.mpf(3)/2)*(4*mp.pi)**(-2))),
                     ('fermion', lamF, -1/(8*mp.pi)/12)):
    amp = fc.amplitudes(2, sig, 'even', mon_e, lam, 'plus')
    # ahat = 2 pi M / (N_T^+ C_T) for the canonically normalised correlator; here unit N_{T,m} = N_T^+/16 * 16...
    a2 = 2*mp.pi*amp[2]*cal/mp.sqrt(kl.norm_formula(3, 2, 2)*NT)
    am2 = 2*mp.pi*amp[-2]*cal/mp.sqrt(kl.norm_formula(3, 2, -2)*NT)
    print(f"   {nm:8s}: ahat_(2,2) = {mp.nstr(mp.chop(a2), 12)},  ahat_(-2,2) = c = {mp.nstr(mp.chop(am2), 12)}")
print("parity-odd structure at s = 5/2:")
mon_o, bas_o = kl.conserved(2, 'odd', sval=sig)
print("   number:", len(bas_o))
for b in bas_o:
    print("   coefficients:", dict((m, c) for m, c in zip(mon_o, b) if c != 0))
    amp = fc.amplitudes(2, sig, 'odd', mon_o, b, 'plus')
    for m in (2, 1, 0, -1, -2):
        print(f"   m = {m:+d}: M = {mp.nstr(mp.chop(amp.get(m, mp.mpc(0))), 15)}")
    print("   => contributes nothing to ahat_(2,2) (Ward identity) and an imaginary amount to c:  2 pi M_{-2}/(N_{2,2} N_T^+)^(1/2) =",
          mp.nstr(mp.chop(2*mp.pi*amp[-2]/mp.sqrt(kl.norm_formula(3, 2, 2)*NT)), 15), " per unit coefficient / C_T")
