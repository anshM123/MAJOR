"""
Audit of ANGULAR_ROUTE §21 items 3-4: numerical stress test of the combination step.

Saturating model: one parity-even operator per unit cell of the (h, J) lattice, h = n + 1/4 (n >= 5), even J >= 4,
with 4^D 2^-J |lambda_{++}|^2 = (s h)^{1/2} (so that every dyadic shell sum is ~ S^{3/2} U^{3/2}, the claimed maximum),
|a0|^2 = |lambda|^2 |K^{++}_0|^2 (C_T = 1), K from the KERNELS_D3 closed form (verified vs the pipeline, ca3).
Then
   mu((L, inf)) = sum_{R > L} C |a0|^2,   R = T/(T - C),  T = (D-1)(D-2), C = J(J+1),
   B_N = ((N+2)(N+3)/6) sum C |a0|^2 R_N,  R_N = (X_N/C) prod_{k<N} (C - k(k+1))/(T - k(k+1)),  X_N = C - N(N+1)  (b = 0).
Checks: L^3 mu((L,inf)) bounded; N B_N bounded; and the intermediate band (tau >= kappa sqrt J) gives L^4 mu bounded.
"""
import numpy as np
from scipy.special import gammaln

def logK2_over_base(D, J, s, h):
    # log of |K^{++}_0|^2 / (4^D 2^-J), h = n + 1/4  (sin^2(pi h) = 1/2)
    lp = (np.log(27) + 3.5*np.log(np.pi) + 4*np.log(2)
          + gammaln(D - 1) - gammaln(D - J - 1)                     # (D-J-1)_J
          + gammaln(D + J) + gammaln(D - 0.5)
          - (gammaln(2*J + 1) - 2*gammaln(J + 1))                   # C(2J, J)
          - 2*gammaln(s) - 2*gammaln(D + 3)
          + np.log(0.5) + 2*gammaln(h) - 2*np.log(np.pi))           # 1/Gamma(1-h)^2
    return lp

Jmax, nmax = 6000, 3000
J = np.arange(4, Jmax + 1, 2, dtype=float)
n = np.arange(5, nmax + 1, dtype=float)
JJ, NN = np.meshgrid(J, n, indexing='ij')
h = NN + 0.25; s = h + JJ; D = s + h
lam2_base = np.sqrt(s*h)                    # 4^D 2^-J |lambda|^2
a02 = lam2_base*np.exp(logK2_over_base(D, JJ, s, h))
T = (D - 1)*(D - 2); C = JJ*(JJ + 1); R = T/(T - C)
w = C*a02
print(f"total sum Delta^2 |a0|^2 = {np.sum(D**2*a02):.4e},  sum C|a0|^2 = {np.sum(w):.4e}")
# convergence of sum Delta^2 |a0|^2 in the cutoff
for cut in (500, 1000, 2000, 3000):
    m = (D < cut)
    print(f"   sum_{{D<{cut}}} D^2|a0|^2 = {np.sum((D**2*a02)[m]):.6e}")
# dyadic shell sums of |a0|^2 vs S^-5 U^2
print("dyadic shells: sum |a0|^2 / (S^-5 U^2)")
for S in (64, 128, 256, 512, 1024):
    row = []
    for U in (8, 32, 128, 512):
        if U > S:
            continue
        m = (s >= S) & (s < 2*S) & (h >= U) & (h < 2*U)
        row.append(f"U={U}:{np.sum(a02[m])/(S**-5.0*U**2):.3g}")
    print(f"   S={S}: " + " ".join(row))
print("tail of mu: L^3 mu((L,inf)) [all tau>=10]   and  L^4 mu_band((L,inf)) [tau >= sqrt(J)]")
band = (2*h >= np.sqrt(JJ))
for L in (2, 4, 8, 16, 32, 64):
    mu = np.sum(w[R > L]); mub = np.sum(w[(R > L) & band])
    print(f"   L={L:3d}: L^3 mu = {L**3*mu:.4g}   L^2 mu = {L**2*mu:.3g}   L^4 mu_band = {L**4*mub:.4g}")
print("B_N (b = 0 ladder), N B_N and N^2 B_N")
logRN = np.zeros_like(C);
for N in range(0, 41):
    if N > 0:
        k = N - 1
        num = C - k*(k + 1); den = T - k*(k + 1)
        ok = num > 0
        logRN = np.where(ok, logRN + np.log(np.where(ok, num, 1)/den), -np.inf)
    XN = C - N*(N + 1)
    RN = np.where(XN > 0, XN/C*np.exp(logRN), 0.0)
    BN = (N + 2)*(N + 3)/6*np.sum(w*RN)
    if N in (1, 2, 4, 8, 16, 24, 32, 40):
        print(f"   N={N:2d}: B_N = {BN:.4e}  N B_N = {N*BN:.4e}  N^2 B_N = {N*N*BN:.4e}")
