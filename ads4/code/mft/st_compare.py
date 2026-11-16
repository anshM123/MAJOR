"""Novelty/consistency check: leading-twist (tau_0) MFT coefficients vs Sleight-Taronna arXiv:1804.09334 eq. (opejphi),
a_(l,0) = 2^{l-1} (tau_1+2)_{l-1} (tau_2)_{l-1} / ((l-1)! (l+tau_1+tau_2)_{l-1}), at tau_1 = 1 (conserved J in d = 3), tau_2 = Dphi.
Prints p[L](N=J-1, J) / a_(J,0) and its product with L_1^2/4 (L = L_1 S_1 on tau_0, so p[S_1] = L_1^2 p[L])."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from mft_exact import closed_even
for Dphi in (sp.Rational(7, 4), sp.Rational(5, 2), sp.Rational(13, 10)):
    r1, r2 = [], []
    for J in range(1, 10):
        D = Dphi + 1 + J
        ours = closed_even(1, D, J, Dphi)
        st = 2**(J-1)*sp.rf(3, J-1)*sp.rf(Dphi, J-1)/(sp.factorial(J-1)*sp.rf(J+1+Dphi, J-1))
        L1 = 2*J*(D-1)*(D-Dphi)
        r1.append(sp.nsimplify(sp.simplify(ours/st))); r2.append(sp.nsimplify(sp.simplify(ours*L1**2/st)))
    print("Dphi = %s:  p[L]/a_ST = %s\n             p[S_1]/a_ST = %s" % (Dphi, r1, r2))
