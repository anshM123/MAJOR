#!/usr/bin/env python3
"""BBP Ising-chain checks (1703.08460 sec 4.1) for U(N_c), N_f = 3N_c+2, and the complex-mass flows
T_M[N_f] -> T_M[N_f-1] that reach smaller N_f.  Checks relevance, unitarity, option-2 exclusion, F monotonicity."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from u2_ct import Grid, UN_logZ_and_d2, chiral_logZ_b1, PI

Fch = lambda R: -chiral_logZ_b1(R)


def FS(Nc, Nf, D):
    R0 = Nf * (1 - D) - (Nc - 1)
    L = min(60.0, 45 / (2 * PI * max(R0, 0.05)) + 2)
    g = Grid(D, L=L, h=0.02, need_p2=False)
    return -UN_logZ_and_d2(g, Nc, Nf)[0]


def TMF(Nc, Nf):
    D = 1 - (Nc + 1) / Nf
    if Nf == Nc + 2:
        return 9 * Fch(2 / 3) if Nc == 1 else Nf * Nf * Fch(0.5)
    F = FS(Nc, Nf, D)
    if 2 * D < 0.5:
        F += Nf * Nf * (Fch(0.5) - Fch(2 * D))
    return F


lines = []
for Nc in [1, 2, 3, 4, 5]:
    Nf = 3 * Nc + 2
    n2 = Nf * Nf
    hi = 1 - (Nc - 1 + 0.03) / Nf
    Rs = minimize_scalar(lambda D: -FS(Nc, Nf, D), bounds=(0.05, hi), method="bounded", options={"xatol": 1e-6}).x
    F0 = FS(Nc, Nf, Rs) + n2 * (Fch(2 / 3) + Fch(0.5))
    F1 = FS(Nc, Nf, 2 / 3) + n2 * (Fch(2 / 3) + Fch(0.5))
    RQ = 1 - (Nc + 1) / Nf
    F2 = FS(Nc, Nf, RQ) + n2 * (Fch(2 - 2 * RQ) + Fch(0.5))
    F3 = FS(Nc, Nf, RQ)
    r = minimize_scalar(lambda D: -(FS(Nc, Nf, D) + n2 * Fch(2 - 2 * D)), bounds=(0.40, hi), method="bounded",
                        options={"xatol": 1e-6})
    R2 = r.x
    s = (f"Nc={Nc} Nf={Nf}: R*(T0)={Rs:.4f}; PhiQQ R={2/3+2*Rs:.3f}; R[V](R_Q=2/3)={Nf/3-Nc+1:.3f}; "
         f"R_Phi(T_M+Phi)={2-2*RQ:.3f}; sigmaPhi R={0.5+2-2*RQ:.3f}\n"
         f"   F chain: {F0:.4f} -> {F1:.4f} -> {F2:.4f} -> {F3:.4f}  monotone={F0 > F1 > F2 > F3}\n"
         f"   option-2 T_(PhiQQ): R_Q={R2:.4f} R_Phi={2-2*R2:.4f} -> Phi^3 {'relevant: option 1' if 3*(2-2*R2) < 2 else 'IRRELEVANT: option 2'}"
         f"; R[V] there={Nf*(1-R2)-Nc+1:.3f}")
    Fm = [TMF(Nc, n) for n in range(Nc + 2, Nf + 2)]
    s += "\n   F(T_M) for N_f=Nc+2..3Nc+3: " + " ".join(f"{x:.3f}" for x in Fm) + \
         f"  increasing in N_f: {all(np.diff(Fm) > 0)}"
    print(s, flush=True)
    lines.append(s)
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bbp_chain_general_output.txt"), "w") as f:
    f.write("\n".join(lines) + "\n")
