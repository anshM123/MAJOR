"""
v1: validation of the CHMRS order-1/m^8 Regge series (ancillary file, parsed by chmrs_series.py) against exact-block Regge
moments computed with the e2e machinery (ads4/code/e2e/moments.py ScalarMoments; exact d = 3 scalar blocks by radial Casimir
recursion; nothing approximated on this side), for IDENTICAL external scalars (D1 = D2 = Df).

measured  M = (-1)^J Pihat_{k,nu}[G]/(2 sin^2)            (e2e normalisation = CHMRS (3.17))
series    M_j = q_{Delta,J} gamma_{2Df+k-1}(nu)^2/(a_Df m^{2k}) * seriesPi truncated at 1/m^{2j}
Output: R_j = M/M_j and m^{2(j+1)}(R_j - 1) (should converge to the next coefficient, i.e. stay bounded).
Also checks the identification al = 2 Df + k - 1 (al = 2Df+k-3/2 is tested as a mutation).
Usage: python3 v1_series_vs_exact.py [quick]
"""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np, mpmath as mp
import skel
from moments import ScalarMoments
from cases import cases_fixed_beta, info

out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
mlist = [10, 20, 40] if quick else [10, 20, 30, 45, 65]
pr("# v1: CHMRS order-8 Regge series vs exact-block moments, identical scalars.  R_j = exact/series_j")
for Df in ((3.0,) if quick else (3.0, 1.75)):
    for beta in ((1.0,) if quick else (0.4, 1.0)):
        pr(f"# Df = {Df}, beta ~ {beta}")
        pr("#  Delta      J     m     k  nu   R_0          R_1          R_2          R_4          m^2(R_0-1)  m^4(R_1-1)  m^6(R_2-1)  mut:R_1(al-1/2)  N=leadratio")
        for (Delta, J) in cases_fixed_beta(beta, mlist):
            t0 = time.time()
            m, b_, eta = info(Delta, J); m = float(m); eta = float(eta)
            N = int(min(12000, max(600, 45*Delta)))
            mom = ScalarMoments(Delta, J, Df, Df, N)
            lq, la = skel.log_q(Delta, J, Df), mp.log(skel.a_Df(Df))
            for k in (4, 6):
                for nu in (0.0, 1.0, 3.0):
                    sg, lg = mom.moment(k, nu)
                    ex = float(sg)*mp.exp(mp.mpf(lg))
                    g2 = skel.gamma_a(2*Df + k - 1, nu)**2
                    pref = mp.exp(lq - la - k*mp.log(m*m))*g2
                    Pd = skel.Pders(nu, eta)
                    t = skel.series_terms('seriesPi', k, Df, nu, eta, Pd=Pd)
                    R = [float(ex/(pref*sum(t[i]/(m*m)**i for i in range(j + 1)))) for j in range(5)]
                    # mutation: al -> al - 1/2 (i.e. the CHMRS 'a' = 2Df+k-d/2) at order 1
                    from chmrs_series import k as ks
                    al_m = 2*Df + k - 1.5
                    tm = [skel._lamb('seriesPi', j)(k, Df, al_m, nu, eta, 0, *Pd) for j in range(2)]
                    Rm = float(ex/(pref*(tm[0] + tm[1]/(m*m))))
                    lead_e2e = mp.exp(mp.log(mp.pi/4) + 2*Df*mp.log(4) - skel.log_b(Delta, J) - (2*Df + k - 1)*mp.log(m*m))*g2*Pd[0]
                    Nrat = float(pref*t[0]/lead_e2e)
                    pr(f"  {Delta:9.4f} {J:4d} {m:6.2f}  {k}  {nu:3.1f}  {R[0]:.9f}  {R[1]:.9f}  {R[2]:.9f}  {R[4]:.9f}  {(R[0]-1)*m**2:10.3f}  {(R[1]-1)*m**4:10.2f}  {(R[2]-1)*m**6:11.1f}  {Rm:.6f}  {Nrat:.8f}")
            pr(f"#    (N = {N}, {time.time()-t0:.0f} s)")
open(os.path.join(HERE, 'cache', 'out_v1%s.txt' % ('_quick' if quick else '')), 'w').write("\n".join(out) + "\n")
