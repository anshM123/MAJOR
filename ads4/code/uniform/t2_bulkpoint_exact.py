"""
t2: exact-block test of the bulk-point (flat-anchored) description at small J and nu/m up to ~0.5 (identical scalars, Df = 3,
the skeleton of <TTTT>; also Df = 1.75).  Exact Pihat moments from the e2e machinery (ScalarMoments, exact d = 3 blocks).

Normalised density  pi_ex = a_Df m^{2k}/gamma_{2Df+k-1}(nu)^2 * Pihat_{k,nu}[Delta,J]  (Pihat[Delta,J] = (-1)^J Pihat[G]/(2 sin^2 q)).
Compared with
  Regge leading     2 P_{-1/2+i nu}(eta_AdS)
  bulk-point (BP)   2 P_J(1-2x^2) / [(1-x^2) ((1+sqrt(1-x^2))/2)^{2k-2}],   x = nu/m      (CHMRS (3.38), normalised)
  CHMRS series      seriesPi at eta_AdS, truncated at 1/m^8
Output columns: x, pi_ex/Regge, pi_ex/BP, m^2 (pi_ex/BP - 1), pi_ex/series.  Precision: the moment is suppressed by
gamma^2 ~ e^{-pi nu} relative to the integrand; with double precision nu <= ~8 is reliable (checked by the series column).
"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
import skel
from moments import ScalarMoments

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def BP_pi(k, x, J):
    return 2*skel.legP(J, 1 - 2*x*x)/((1 - x*x)*((1 + np.sqrt(1 - x*x))/2)**(2*k - 2))

pr("# t2: exact Pihat (normalised) vs Regge leading, bulk-point formula and CHMRS series; small J (bulk-point regime)")
pr("#  Df   J   Delta     m      k  nu    x      ex/Regge     ex/BP        m^2(ex/BP-1)   ex/series8")
for Df in (3.0, 1.75):
    for J in (0, 2, 4):
        for mt in (12, 16, 24, 32):
            Delta = float(skel.Delta_of(mt*mt, J))
            m2, eta = skel.m2_eta(Delta, J); m = np.sqrt(m2)
            N = int(min(12000, max(800, 60*Delta)))
            t0 = time.time()
            mom = ScalarMoments(Delta, J, Df, Df, N)
            lq, la = skel.log_q(Delta, J, Df), mp.log(skel.a_Df(Df))
            for k in (4, 6):
                for nu in (0.0, 2.0, 4.0, 6.0, 8.0):
                    if nu > 0.6*m:
                        continue
                    sg, lg = mom.moment(k, nu)
                    ex = float(sg)*mp.exp(mp.mpf(lg))
                    pi_ex = float(ex*mp.exp(la - lq + k*mp.log(m2))/skel.gamma_a(2*Df + k - 1, nu)**2)     # (-1)^J already in the e2e moment
                    Pd = skel.Pders(nu, eta)
                    rg = 2*Pd[0]
                    x = nu/m
                    bp = BP_pi(k, x, J)
                    ser = skel.series_value('seriesPi', k, Df, nu, Delta, J, order=4, Pd=Pd)
                    pr(f"  {Df:4.2f} {J:2d} {Delta:8.3f} {m:6.2f}  {k}  {nu:3.0f}  {x:.3f}  {pi_ex/rg:11.6f}  {pi_ex/bp:11.6f}  {m2*(pi_ex/bp-1):11.3f}   {pi_ex/ser:11.7f}")
            pr(f"#   (N = {N}, {time.time()-t0:.1f} s)")
open(os.path.join(HERE, 'cache', 'out_t2.txt'), 'w').write("\n".join(out) + "\n")
