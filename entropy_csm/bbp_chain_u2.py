#!/usr/bin/env python3
"""F-theorem / relevance checks along BBP's Ising-SCFT chain (1703.08460 sec 4.1) for U(2), N_f=8 = 3N_c+2:
T0 + 64 Ising(Phi^3) + 64 sigma -> T_{PhiQQt+Phi^3} -> T_{V+PhiQQt} -> T_M."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from u2_ct import Grid, UN_logZ_and_d2, chiral_logZ_b1

Fch = lambda R: -chiral_logZ_b1(R)
def FS(D):
    g = Grid(D, L=10.0, h=0.01, need_p2=False)
    return -UN_logZ_and_d2(g, 2, 8)[0]
Nf2 = 64
Rs = minimize_scalar(lambda D: -FS(D), bounds=(0.3, 0.6), method="bounded", options={"xatol": 1e-7}).x
F0 = FS(Rs) + Nf2 * Fch(2 / 3) + Nf2 * Fch(0.5)
F1 = FS(2 / 3) + Nf2 * Fch(2 / 3) + Nf2 * Fch(0.5)
F2 = FS(5 / 8) + Nf2 * Fch(3 / 4) + Nf2 * Fch(0.5)
F3 = FS(5 / 8)
print(f"R_Q*(T0)={Rs:.6f}")
print(f"F: T0+Ising+sigma {F0:.6f} -> T_(PhiQQ+Phi3)+sigma {F1:.6f} -> T_(V+PhiQQ)+sigma {F2:.6f} -> T_M {F3:.6f}")
print("monotone decreasing:", F0 > F1 > F2 > F3)
print("relevance: PhiQQt at T0+Ising: R =", 2 / 3 + 2 * Rs, "; V at R_Q=2/3: R =", 8 * (1 / 3) - 1,
      "; sigma*Phi at R_Phi=3/4: R =", 0.5 + 0.75)
# option 2 of BBP: T_{PhiQQt} with Phi^3 dropped; F-maximise SQCD + 64 Phi (R_Phi = 2 - 2R)
r = minimize_scalar(lambda D: -(FS(D) + Nf2 * Fch(2 - 2 * D)), bounds=(0.40, 0.85), method="bounded",
                    options={"xatol": 1e-7})
R2 = r.x
print(f"T_(PhiQQt) without Phi^3: R_Q*={R2:.6f}, R_Phi={2-2*R2:.6f} -> Phi^3 {'relevant (option 1 selected)' if 3*(2-2*R2) < 2 else 'irrelevant (option 2)'}; "
      f"R[V] there = {8*(1-R2)-1:.4f} -> {'relevant' if 8*(1-R2)-1 < 2 else 'irrelevant'}; F={-r.fun:.6f}")
print("unitarity at intermediate points: R[V](R_Q=2/3)=", 8/3-1, " R_Phi>=1/2 ok; mesons flipped")
