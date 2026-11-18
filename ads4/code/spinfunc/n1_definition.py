"""
n1_definition (doc 24.1): the CHMRS magic-block coefficient g(nu) (l2_scalar_validate.gfac, a = 7, d = 3) and what it implies for the
definition of Omega^{A'}:
  * Re g(nu) (= (g + g~)/2 for real nu) has a zero nu_0: the heavy-limit (fixed nu, m -> oo) normalisation of the literal CLS-D k = 2
    functional, N(nu)(g + g~) C_heavy, vanishes there, so Psi^{CLS}/N^(nu) (doc R.3) has a pole at p_0 = nu_0/Lambda;
  * the leading Regge kernel of Psi^{CLS} implied by the identity, [g P_{l-1} + g~ P_{l+1}]/(g + g~), vs the eta P_l assumed by the doc-R.5
    eta-inversion: sigma(nu) := Im g/(g + g~) vs 1/(4 nu).
Usage: python3 n1_definition.py -> cache/n1_definition.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np, mpmath as mp
import l2_scalar_validate as V

def reg(nu):
    gp, gm = V.gfac(nu); return (gp + gm).real

lines = []
def pr(s=""):
    print(s, flush=True); lines.append(s)
pr("# g(nu) = ((2-d)/2 + i nu)/(i nu) (gamma_7(nu+i)/gamma_7(nu))^2 (CHMRS; l2_scalar_validate.gfac)")
pr("  nu      g+g~ (real)   Im g/(g+g~)   1/(4nu)   nu*Im g/(g+g~)")
for nu in (0.5, 1, 2, 3, 4, 5, 5.5, 6, 6.5, 7, 8, 10, 15, 20, 30, 50, 100, 300):
    gp, gm = V.gfac(nu)
    s = (gp/(gp + gm)).imag
    pr(f"  {nu:6.1f}  {(gp+gm).real:+10.5f}   {s:+12.5f}   {1/(4*nu):8.5f}   {nu*s:+9.4f}")
nu0 = float(mp.findroot(lambda v: reg(float(v)), 6.0))
pr(f"zero of g + g~:  nu_0 = {nu0:.6f}")
for Dg in (100, 300, 1000):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0))
    pr(f"  Delta_g = {Dg}: p_0 = nu_0/Lambda = {nu0/Lam:.5f}")
# asymptotic coefficient
gp, gm = V.gfac(2000.0)
pr(f"large nu: nu Im g/(g+g~) at nu = 2000: {2000*(gp/(gp+gm)).imag:.4f} (eta P_l would give 0.25)")
open(os.path.join(HERE, 'cache', 'n1_definition.txt'), 'w').write("\n".join(lines) + "\n")
