"""Result Q, Stage 7 (RESULT_Q_STAGE7.md §4.2): validity checks of the extended-grid C4 weights (functionals_q7_xprec.json) in double precision:
shifted light-side moments M_n = Re sum_j int chi R_j (nu - ij)^{2n} (targets: value (1,0,0,0), derivative (0,1,0,0)), R^4 normalisation
int chi Nhat, Stage-6 covariant test (flat scalar skeleton g at five states: image must be g(0) / g'(0)), amplification A[1], A[|nu-ij|^2]
and the Stage-6 stress coefficients (s_a, s_b), compared with the stored H0 weights.  Output: cache/stage7_xpcheck.txt"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import stage4_hj as HJ4
import stage6_gauge as G6

F4 = HJ4.numfuncs('EE3')
lines = []


def pr(s):
    print(s, flush=True); lines.append(s)


def check(lab, chi, dnu, deriv):
    nu = (np.arange(len(chi)) + 0.5)*dnu
    R = {j: F4[j](nu.astype(complex)) for j in F4}; n = sum(np.abs(R[j]) for j in R); Nh = sum(R.values()).real
    Mn = [dnu*np.sum(chi*sum(R[j]*(nu - 1j*j)**(2*k) for j in R)).real for k in range(4)]
    A0 = dnu*np.sum(np.abs(chi)*n); A2 = dnu*np.sum(np.abs(chi)*sum(np.abs(R[j])*np.abs(nu - 1j*j)**2 for j in R))
    mu1 = dnu*np.sum(chi*Nh*nu**2); sa = mu1 - Mn[1]; sb = dnu*np.sum((chi*sum((j*j/9.)*R[j] for j in R)).real*nu**2)
    cov = []
    for (J, m) in ((5, 1.2*298.5), (4, 298.5*1.013), (8, 1.1*298.5), (24, 1.3*298.5), (32, 3*298.5)):   # Dg = 300 (m > nu_max)
        g = lambda mu2: G6.gfun(J, m, mu2)
        img = dnu*np.sum(chi*sum(R[j]*g((nu - 1j*j)**2) for j in R)).real
        eps = 1e-6*m*m
        tgt = (g(eps) - g(-eps)).real/(2*eps) if deriv else g(0.0).real
        cov.append(abs(img - tgt))
    pr(f"  {lab:28s}: int chi Nhat = {dnu*np.sum(chi*Nh):+.12f}; M0..3 = " + " ".join(f"{x:+.2e}" for x in Mn) +
       f"; A[1] = {A0:.4f}, A[|nu-ij|^2] = {A2:.3f}; s_a = {sa:+.4f}, s_b = {sb:+.4f}; covariant test max err {max(cov):.1e}")


st = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']
xp = json.load(open(os.path.join(HERE, 'functionals_q7_xprec.json')))
pr("# stage7_xpcheck: validity of the extended-grid (dnu = 1/25, N = 3072) C4 weights vs the stored (dnu = 1/50, N = 2048) ones")
check('stored value (H0)', np.array(st['chi_value']), 0.02, False)
check('XP value', np.array([float(x) for x in xp['C4_value']['chi']]), 0.04, False)
check('stored derivative (H0)', np.array(st['chi_deriv']), 0.02, True)
check('XP derivative', np.array([float(x) for x in xp['C4_deriv']['chi']]), 0.04, True)
open(os.path.join(HERE, 'cache', 'stage7_xpcheck.txt'), 'w').write("\n".join(lines) + "\n")
