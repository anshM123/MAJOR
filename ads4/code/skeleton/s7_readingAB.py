"""
s7 (follow-up 3): Reading A vs Reading B for the opposite-helicity (J_z = 2) component of <JJJJ>, at FIXED impact parameter
(fixed beta, i.e. b-hat = (2J+1)/m = 2 sinh(beta/2) fixed, J grows with m) and fixed x = nu/m, as m grows.

Prediction for the exact MFT-normalised opp moment (odd kernel; its leading value is purely n00 P):
    Kpred_D = S0 * Nrat * R_skel(Df) * Kp * D(J, x)/P_J(1 - 2x^2)        (flat density D in place of the scalar P_J)
  A : D = P_J                                    (scalar flat density; = s4 c = 0)
  B : D = d^J_{22}(theta) = cos^4(theta/2) P^{(0,4)}_{J-2}(cos theta),  cos theta = 1 - 2x^2   (photon +- Wigner function, J_z = 2)
  Bt: D = d~^J_{22} = P^{(0,4)}_{J-2}            (the reduced function with the scalar prefactor, i.e. F_flat's D_{+-} literally)
  C : curved helicity shift eta -> eta - 8/m^2 (s4 c = 8), for comparison
  d_D = (Q R Kp - Kpred_D)/(A E_k/2) as in s4;  only points where P_J(1-2x^2) is not near a zero (|P_J|/amp >= 0.3) are used for A, B, Bt.
Prints m^2 d and m d at x = 0.05, 0.1, 0.15, 0.2 per (beta, m), k = 4, 6, Df = 2 and 2.5.  Writes cache/out_s7.txt.
"""
import os, pickle
from common import CACHE, info
import numpy as np
from scipy.special import eval_jacobi, eval_legendre
import spin
from skeleton import Skeleton, Nrat
from s3_compare import envelope
from s4_summary import cases, Ek, TAILMAX, dPc

def D_of(name, J, x):
    z = 1 - 2*x*x
    if name == 'A':
        return eval_legendre(J, z)
    dt = eval_jacobi(J - 2, 0, 4, z)
    return dt if name == 'Bt' else (1 - x*x)**2*dt

def PJ_cond(J, x):
    z = 1 - 2*x*x
    P = eval_legendre(J, z)
    dP = J*(z*eval_legendre(J, z) - eval_legendre(J - 1, z))/(z*z - 1) if J > 0 else 0.0
    amp = np.sqrt(P*P + ((1 - z*z)*dP*dP)/(J + 0.5)**2)
    return abs(P)/amp

out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

XS = (0.05, 0.1, 0.15, 0.2)
pr("# s7: opp helicity, odd kernel. m^2 d_D (and m d_D in brackets) at x = " + ", ".join(map(str, XS)) + ";  '*' = P_J near a zero (skipped)")
for k in (6, 4):
    for Df in (2.0, 2.5):
        pr(f"\n## k={k} Df={Df}")
        pr("#  beta     m    J   | " + " | ".join(f"x={x}: A / B / Bt / C" for x in XS))
        for (D, J, m, beta) in cases():
            sp = pickle.load(open(os.path.join(CACHE, f"s1_D{D:.6f}_J{J}_std.pkl"), 'rb'))
            sk = Skeleton(D, J, 'std')
            NR = Nrat(D, J, Df); Q = spin.Q(D, J, 'opp')
            if (k, 0.0) not in sp:
                continue
            R0, _, Kp0 = sp[(k, 0.0)][('odd', 'opp')]
            S0 = Q*R0/(NR*sk.R(Df, k, 0.0))
            cells = []
            for x0 in XS:
                nus = [nu for (kk, nu) in sp if kk == k and abs(nu/m - x0) < 1e-3]
                if not nus or sp[(k, nus[0])][('odd', 'opp')][1] >= TAILMAX:
                    cells.append("        -        "); continue
                nu = nus[0]; x = nu/m
                R, tail, Kp = sp[(k, nu)][('odd', 'opp')]
                _, A, n00, P = envelope(D, J, k, nu, 'odd', 'opp')
                base = S0*NR*sk.R(Df, k, nu)
                vals = []
                ok = PJ_cond(J, x) >= 0.3
                for name in ('A', 'B', 'Bt'):
                    if not ok:
                        vals.append("   *"); continue
                    Kd = base*Kp*D_of(name, J, x)/eval_legendre(J, 1 - 2*x*x)
                    d = (Q*R*Kp - Kd)/(A*Ek(k, x)/2)
                    vals.append(f"{m*m*d:6.1f}")
                dC = (Q*R*Kp - base*(Kp + n00*dPc(8, nu, D, J)))/(A*Ek(k, x)/2)
                vals.append(f"{m*m*dC:6.1f}")
                cells.append("/".join(vals))
            sk.save()
            pr(f"   {beta:.2f} {m:6.1f} {J:4d} | " + " | ".join(cells))
open(os.path.join(CACHE, 'out_s7.txt'), 'w').write("\n".join(out) + "\n")
