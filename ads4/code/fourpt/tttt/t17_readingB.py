"""
H2 of RESULT_D.md for gravitons: Reading A / B / B~ for the opposite-helicity <TTTT> moment at FIXED impact parameter (fixed beta = 1,
b-hat = (2J+1)/m = 2 sinh(1/2) = 1.04, J grows with m), the analogue of SKELETON_TEST.md sect. 8 (currents, s7_readingAB.py).
Observable: the exact opposite-helicity wave W[X_o,X_o] with the odd kernel A x A (t16; its leading value is purely the preserved
component, ~ P_{-1/2+i nu}(cosh beta)), R(nu) = measured/leading prediction.  At fixed x = nu/m:
   r(x) = R(nu) / (R(0) * rho_skel(nu)),   rho_skel = R_skel(Df,k,nu)/R_skel(Df,k,0)  (exact scalar skeleton, skeleton/skeleton.py; it removes
   the scalar O(nu^2/m^2) structure common to all readings; Df = 3 = Delta_T, and Df = 2 as a check),
and the readings predict r = D(J,x)/P_J(1 - 2x^2) with
   A : D = P_J (scalar);  B : D = d^J_{44}(theta) = cos^8(theta/2) P^{(0,8)}_{J-4}(cos theta), cos theta = 1 - 2x^2 (full graviton Wigner function,
   J_z = 4);  B~ : D = d~^J_{44} = P^{(0,8)}_{J-4} (reduced function, F_flat's D_{+-} literally).
Printed: m^2 d_X with d_X = r - D_X/P_J  (multiplicative at fixed b-hat), and m d_A; points with P_J(1-2x^2) near a zero are skipped.
Usage: python3 t17_readingB.py   Output: out_t17_readingB.txt
"""
import os, sys, glob
import numpy as np, mpmath as mp
from scipy.special import eval_jacobi, eval_legendre
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', 'oddk'))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'skeleton')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from t13_analysis_l2 import Modes, design, Tpred, CLASSES
from t16_kernels_ext import KER
from predict4 import weyl
from skeleton import Skeleton
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
XS = (0.05, 0.1)
K_ = 6

def D_of(name, J, x):
    z = 1 - 2*x*x
    if name == 'A':
        return eval_legendre(J, z)
    dt = eval_jacobi(J - 4, 0, 8, z)
    return dt if name == 'Bt' else (1 - x*x)**4*dt

def PJ_cond(J, x):
    z = 1 - 2*x*x
    P = eval_legendre(J, z)
    dP = J*(z*eval_legendre(J, z) - eval_legendre(J - 1, z))/(z*z - 1)
    amp = np.sqrt(P*P + ((1 - z*z)*dP*dP)/(J + 0.5)**2)
    return abs(P)/amp

def R_opp(mo, k, nu, kern='O'):
    D, J = mo.Delta, mo.J
    A = design(D, J, k, nu); Tp = Tpred('opp', 'opp', D, J, 8.0)
    K = KER[kern]; terms = []; pred = 0j
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4)
            v, sh = mo.moment(k, nu, ('opp', 'opp'), q)
            terms.append((w13*w24, v, sh)); pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    shm = max(t[2] for t in terms)
    meas = sum(w*v*np.exp(sh - shm) for w, v, sh in terms)*(-1)**J
    return complex(meas*np.exp(shm - np.log(abs(pred)))/(pred/abs(pred))).real

if __name__ == "__main__":
    files = [f for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'modes_l2_*_d.npz')))]
    pr("# H2 for gravitons: reading A / B / B~ at fixed b-hat (fourpt/tttt/t17_readingB.py); opp helicity, odd kernel A x A, k = 6")
    pr("#  beta    m    J  Df | x: nu, R(nu), r, m^2 d_A, m^2 d_B, m^2 d_Bt, (m d_A), [P_J condition]")
    for fn in files:
        mo = Modes(fn)
        Ap, Am, m, beta = weyl(mo.Delta, mo.J); m = float(m); beta = float(beta)
        if abs(beta - 1.0) > 0.05: continue
        sk = Skeleton(mo.Delta, mo.J, 'std')
        R0 = R_opp(mo, K_, 0.0)
        for Df in (3.0, 2.0):
            s0 = sk.R(Df, K_, 0.0)
            cells = []
            for x in XS:
                nu = round(x*m, 6)
                R = R_opp(mo, K_, nu)
                r = R/(R0*sk.R(Df, K_, nu)/s0)
                c = PJ_cond(mo.J, x)
                P = D_of('A', mo.J, x)
                dA = r - 1.0; dB = r - D_of('B', mo.J, x)/P; dBt = r - D_of('Bt', mo.J, x)/P
                star = '' if c >= 0.3 else ' *near zero*'
                cells.append(f"x={x}: nu={nu:.3f} R={R:.4f} r={r:.4f} A {m*m*dA:+.1f} B {m*m*dB:+.1f} Bt {m*m*dBt:+.1f} (m dA {m*dA:+.3f}) [{c:.2f}]{star}")
            pr(f"  {beta:.2f} {m:6.1f} {mo.J:4d} {Df:.0f} | " + " | ".join(cells))
            sk.save()
    pr("# predicted flat helicity term at fixed b: d^J_44/P_J - 1 -> 2*16 (x^2/b) J1(b)/J0(b) (lambda = 4), i.e. O(1/m) at fixed b-hat for reading A")
    open(os.path.join(HERE, 'out_t17_readingB.txt'), 'w').write("\n".join(out) + "\n")
