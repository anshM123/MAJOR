"""
Stage S: end-to-end test of the heavy-block Regge-moment formula for scalars with pairwise different dimensions,
<phi_1 phi_2 phi_1 phi_2> (dimensions D1, D2; s-channel 12 -> 34), which is the scalar skeleton of the spinning case.
Exact:      (-1)^J Pihat_{k,nu}[ u'^{-(D1+D2)/2} g^{a,b}_{Delta,J} ],   a = -(D1-D2)/2 = -b   (no approximation).
Prediction (derived in SPINNING_E2E.md sect. 4.1; no free constant):
            (pi/4) 4^{D1+D2} / b^{(D34)}_{Delta,J} * gamma_{2 D2+k-1}(nu) gamma_{2 D1+k-1}(nu) * P_{-1/2+i nu}(eta_AdS) / m^{2(D1+D2+k-1)},
with the shadow coefficient b^{(D34)} of formulas.log_b (D34 = D1 - D2).  For comparison the same with b^{(0)} (which
differs at relative order 1/m).
Usage: python3 t3_scalar_mixed.py [quick]
"""
import sys, time, numpy as np, mpmath as mp
from moments import ScalarMoments
from formulas import gamma_a, log_b, Pnu, m2_eta
from cases import cases_fixed_beta, info
mp.mp.dps = 25
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def prediction(Delta, J, D1, D2, k, nu, D34=None):
    Delta = mp.mpf(Delta)
    m2, eta = m2_eta(Delta, J)
    ax, ay = 2*D2 + k - 1, 2*D1 + k - 1
    D34 = (D1 - D2) if D34 is None else D34
    return mp.exp(mp.log(mp.pi/4) + (D1 + D2)*mp.log(4) - log_b(Delta, J, D34) - (ax + ay)/2*mp.log(m2))*gamma_a(ax, nu)*gamma_a(ay, nu)*Pnu(nu, eta)

if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    mlist = [10, 20, 40] if quick else [10, 20, 40, 80, 120]
    for (D1, D2) in ((3.0, 1.75), (2.0, 1.75), (1.75, 3.0)):
        pr(f"# D1 = {D1}, D2 = {D2};  R = exact/prediction,  R0 = the same with b^(D34=0)")
        pr("#  Delta        J      m      beta   k   nu       exact (-1)^J Pihat       R           m^2 (R-1)     R0        m (R0-1)")
        for beta in ((1.0,) if quick else (0.4, 1.0, 2.0)):
            for (Delta, J) in cases_fixed_beta(beta, mlist):
                t0 = time.time()
                m, b_, eta = info(Delta, J)
                N = int(min(12000, max(600, 45*Delta)))
                mom = ScalarMoments(Delta, J, D1, D2, N)
                for k in (2,):
                    for nu in (0.0, 1.0, 3.0):
                        sg, lg = mom.moment(k, nu)
                        ex = sg*mp.exp(mp.mpf(lg))
                        R = ex/prediction(Delta, J, D1, D2, k, nu).real
                        R0 = ex/prediction(Delta, J, D1, D2, k, nu, D34=0).real
                        pr(f"  {Delta:10.5f} {J:4d} {float(m):7.2f} {float(b_):6.3f}  {k}  {nu:4.1f}  {mp.nstr(ex, 12):>22s}  {float(R):11.8f}  {float((R-1)*m*m):10.4f}  {float(R0):10.6f}  {float((R0-1)*m):9.4f}")
                pr(f"#    (N = {N}, {time.time()-t0:.0f} s)")
    if not quick:
        open('out_t3_scalar_mixed.txt', 'w').write("\n".join(out) + "\n")
