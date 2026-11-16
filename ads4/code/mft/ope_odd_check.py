"""
Block-level test of the parity-odd normalisation and of its sign (referee issues 1, 2), l = 1, Dphi = 7/4, real section.

OPE expansion x_12 = delta v -> 0 of the MFT correlator F = <O_1 O_3><phi_2 phi_4> (|x_ij^2| powers, as in lightcone_check.py).
Primaries with N <= 1 (Delta = Delta_O + Dphi + N):  (N,J) = (0,1) even,  (1,1) odd,  (1,2) even.  Hence
   R(delta) := F - sum_{even, N<=1} p[L] (-1/2)^J W_diff   [exact e2e blocks, Fock coefficients]
             = (odd (1,1) leading OPE term) + O(delta^2 relative to F),
and the leading OPE term of an exchange with <O phi X> = lam S is  lam^2 (-1)^J <S, S>  with the spin-J contraction <,> of
SPINNING_E2E.md sect. 2.3 (t2_spinning_blocks.ope_check), the rule that the even exchanges obey (W_diff = 8^J <sum c W, sum c W>,
sum_w c_w W_w = 2^{-J} L).  For the odd structure (real section, Lorentzian eps)
   Ob = eps(Z1,P1,Z5,P2,P5) v5^{J-1} / (|s12|^al |s15|^be |s25|^ga),  exponents of S3(2, Dphi),
the Fock space gives |lam|^2 = p[Ob] = 2^{3J-3} P_FK(2,Dphi) H_1 (= 7/11 at (1,1)), and the test measures
   sigma := R / (p[Ob] (-1)^J <Ob, Ob>)       (prediction: |sigma| = 1; its sign is the real-section sign of the odd exchange).
A control is the same with the (1,2) even subtraction omitted (R must then not be of order delta).
"""
import os, sys
import numpy as np, mpmath as mp, sympy as sp
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here); sys.path.insert(0, os.path.join(here, '..', 'e2e'))
from spinblock import SpinBlock
from lightcone_check import wfull, F_mft
from mft_exact import closed_even, closed_odd_l1
mp.mp.dps = 40
DPHI = mp.mpf(7)/4
md = lambda u, v: -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def P(x): return mp.matrix([1, md(x, x)] + list(x))
def Zv(x, z): return mp.matrix([0, 2*md(x, z)] + list(z))
def edot(A, B): return -(A[0]*B[1] + A[1]*B[0])/2 + md(A[2:5], B[2:5])
def eps5(*vs):
    M = mp.matrix(5, 5)
    for j, v in enumerate(vs):
        for i in range(5): M[i, j] = v[i]
    return mp.det(M)

def odd_struct(Delta, J):
    al = (2 + DPHI - Delta + J)/2; be = (Delta + J + 2 - DPHI)/2; ga = (Delta + J - 2 + DPHI)/2
    def S(x1, x2, x5, z1, z5):
        x1, x2, x5 = [list(v) for v in (x1, x2, x5)]
        P1, P2, P5, Z1, Z5 = P(x1), P(x2), P(x5), Zv(x1, z1), Zv(x5, z5)
        v5 = edot(Z5, P1)*edot(P2, P5) - edot(Z5, P2)*edot(P1, P5)
        s = lambda A, B: abs(2*edot(A, B))
        return eps5(Z1, P1, Z5, P2, P5)*v5**(J - 1)/(s(P1, P2)**al*s(P1, P5)**be*s(P2, P5)**ga)
    return S

def refl(x, n):
    x2 = md(x, x); c = md(x, n)
    return [n[i] - 2*x[i]*c/x2 for i in range(3)]

def pairing(S, Delta, J, x1, x2, x3, x4, z1, z3):
    """<C_L, R>: leading OPE tensor of the left structure (x_5 -> infinity) contracted with the right structure (exchanged at x_2)."""
    x5 = [mp.mpf('1e16'), mp.mpf('3e15'), mp.mpf('-2e15')]
    M = 2*J + 1
    phis = [2*mp.pi*k/M for k in range(M)]
    x25 = [x2[i] - x5[i] for i in range(3)]; s25 = -md(x25, x25)
    CL, RR = [], []
    with mp.workdps(100):
        for ph in phis:
            n = [mp.mpf(1), mp.cos(ph), mp.sin(ph)]
            CL.append(s25**Delta*S(x1, x2, x5, z1, refl(x25, n)))
            RR.append(S(x3, x4, x2, z3, n))
    am = {m: sum(CL[k]*mp.exp(-1j*m*phis[k]) for k in range(M))/M for m in range(-J, J+1)}
    bm = {m: sum(RR[k]*mp.exp(-1j*m*phis[k]) for k in range(M))/M for m in range(-J, J+1)}
    return mp.re(sum(am[-m]*bm[m]*2**J*(-1)**(J-m)/mp.binomial(2*J, J-m) for m in range(-J, J+1)))

def main():
    v = [mp.mpf('1.0'), mp.mpf('0.3'), mp.mpf('-0.2')]
    x2 = [mp.mpf(0)]*3
    x3 = [mp.mpf('1.3'), mp.mpf('0.2'), mp.mpf('0.4')]; x4 = [mp.mpf('2.9'), mp.mpf('-0.3'), mp.mpf('0.5')]
    pols = [((1, 0.7), (1, 2.1)), ((1, 1.9), (1, 0.4)), ((1, 3.0), (1, 4.4))]
    Dp = 1.75
    E = [(0, 1), (1, 2)]
    pE = {(N, J): float(closed_even(1, sp.Rational(15, 4) + N, J, sp.Rational(7, 4))) for (N, J) in E}
    DeltaO, JO = mp.mpf(19)/4, 1
    pO = float(closed_odd_l1(sp.Rational(19, 4), 1, sp.Rational(7, 4)))*4**(JO - 1)     # p[Ob] = 4^{J-1} p[eps(z1,zeta,x)(zeta.x)^{J-1}...]
    S = odd_struct(DeltaO, JO)
    print("p[L](0,1) = %s, p[L](1,2) = %s, p[Ob](1,1) = %s  (Fock-space values, Dphi = 7/4)" % (pE[(0, 1)], pE[(1, 2)], pO))
    print("| polarisations (th1, th3) | delta | R/F | sigma = R/(p (-1)^J <Ob,Ob>) | control: R without (1,2) / F |")
    print("|---|---|---|---|---|")
    for (a, t1), (b, t3) in pols:
        z1 = [mp.mpf(1), mp.cos(t1), mp.sin(t1)]; z3 = [mp.mpf(1), mp.cos(t3), mp.sin(t3)]
        for dl in ('1e-2', '1e-3', '1e-4'):
            d = mp.mpf(dl)
            x1 = [x2[i] + d*v[i] for i in range(3)]
            fl = lambda w: [float(c) for c in w]
            F = F_mft(1, fl(x1), fl(x2), fl(x3), fl(x4), [fl(z1)], [fl(z3)])[0, 0]
            Ew = {}
            for (N, J) in E:
                W, _ = wfull(SpinBlock(1, 3.75 + N, J, Dp, 80), fl(x1), fl(x2), fl(x3), fl(x4), [fl(z1)], [fl(z3)])
                Ew[(N, J)] = pE[(N, J)]*(-0.5)**J*W[0, 0]
            R = F - Ew[(0, 1)] - Ew[(1, 2)]
            pred = pO*(-1)**JO*pairing(S, DeltaO, JO, x1, x2, x3, x4, z1, z3)
            print("| (%.1f, %.1f) | %s | %.3e | %.6f | %.3e |" % (t1, t3, dl, R/F, float(R/pred), (F - Ew[(0, 1)])/F), flush=True)

if __name__ == "__main__":
    import warnings; warnings.filterwarnings('ignore')
    main()
