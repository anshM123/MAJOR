#!/usr/bin/env python3
"""Audit of the U(2), N_f=8, W=V_+ + V_- candidate (and its BBP dual U(4)+M). Uses u2_ct.py."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from u2_ct import *

F28 = 28 * np.log(2) / 2
out = []


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    out.append(s)


def U2_data(D, L=8.0, h=0.01, Nf=8, N=2):
    g = Grid(D, L=L, h=h)
    lz, d1, d2, sg = UN_logZ_and_d2(g, N, Nf)
    return g, -lz, CT_from_d2logZ(d2), d1, sg


# ---------------- calibrations (repeat key ones)
say("== calibrations ==")
say("free chiral C_T:", CT_from_d2logZ(chiral_d2(0.5)), "(4)")
say("XYZ C_T:", 3 * CT_from_d2logZ(chiral_d2(2 / 3)), "(8.7214)")
g = Grid(0.5)
for n in [2, 4]:
    lz, d1, d2, _ = UN_logZ_and_d2(g, 1, n)
    ex = 12 * (2 * n * n * polygamma(1, n / 2 + 1) + (PI ** 2 * n + 4) * n + 4) / (PI ** 2 * (n + 1)) / 1.5
    say(f"N=4 SQED_{n}: C_T={CT_from_d2logZ(d2):.10f}  1910.03600: {ex:.10f}")
s, h = g.s, g.h
w = 0.5 * 4 * np.sinh(2 * PI * s) ** 2 * np.exp(6 * g.lp)
obs = np.array([V2(2 * x) for x in s]) + 6 * g.p2
say(f"N=4 SU(2) w/ 3 flavours (vector-multiplet test; must equal SQED_4 by Sp(1)N_f=3 = U(1)N_f=4): "
    f"{CT_from_d2logZ(h * (w * obs).sum() / (h * w.sum())):.10f}")

# ---------------- 1+2: the candidate
say("\n== U(2), N_f=8, R_Q=5/8 ==")
res = {}
for (L, hh) in [(8.0, 0.01), (10.0, 0.005)]:
    g, F, CT, d1, sg = U2_data(5 / 8, L=L, h=hh)
    res[(L, hh)] = (F, CT)
    if hh >= 0.01:
        lzd, d2d = U2_direct(g, 8)
        say(f"   2D direct quadrature: C_T={CT_from_d2logZ(d2d):.10f}  F={-lzd:.12f}")
    say(f"L={L} h={hh}: F={F:.12f}  C_T(det)={CT:.10f}  dlogZ/db={d1:.1e}  sign det={sg}")
F_c, CT_c = res[(10.0, 0.005)]
say(f"claimed C_T = 92.01494393 ; difference = {CT_c - 92.01494393:.3e}")
say(f"F = {F_c:.10f}  vs  F(28 free chirals) = {F28:.10f}  -> F-theorem {'PASS' if F_c > F28 else 'FAIL'}")

# ---------------- plain U(2) N_f=8 SQCD: F-maximise
say("\n== plain U(2) N_f=8 (no W): F-maximisation over R_Q ==")
Fof = lambda D: U2_data(D, L=8.0, h=0.01)[1]
r = minimize_scalar(lambda D: -Fof(D), bounds=(0.15, 0.85), method="bounded", options={"xatol": 1e-7})
Dstar = r.x
g, Fs, CTs, _, _ = U2_data(Dstar)
say(f"R_Q* = {Dstar:.7f}  F* = {Fs:.8f}  C_T = {CTs:.6f}  R[V_pm] = 8(1-R_Q)-1 = {8*(1-Dstar)-1:.6f} "
    f"-> monopole W {'RELEVANT' if 8*(1-Dstar)-1 < 2 else 'IRRELEVANT'} at electric SQCD point")
say(f"   (F at R_Q=5/8 = {F_c:.8f}; F-max requires F* >= F(5/8): {Fs >= F_c - 1e-9})")

# ---------------- dual: U(4), N_f=8, + 64 M, W = M q qt (+ monopoles)
say("\n== dual U(4) N_f=8 + 64 singlets M ==")


def U4M(Dq, L=8.0, h=0.01):
    g = Grid(Dq, L=L, h=h)
    lz, d1, d2, sg = UN_logZ_and_d2(g, 4, 8)
    DM = 2 - 2 * Dq
    F = -lz - 64 * chiral_logZ_b1(DM)
    CT = CT_from_d2logZ(d2) + 64 * CT_from_d2logZ(chiral_d2(DM))
    return F, CT, d1, sg


Fd, CTd, d1, sg = U4M(3 / 8, L=10.0, h=0.005)
say(f"at R_q=3/8 (BBP point): F={Fd:.12f}  C_T={CTd:.10f}  (sign det {sg}, dlogZ/db {d1:.1e})")
say(f"BBP check: F_U2 - F_U4M = {F_c - Fd:.3e} ; C_T diff = {CT_c - CTd:.3e}")
r = minimize_scalar(lambda D: -U4M(D)[0], bounds=(0.15, 0.62), method="bounded", options={"xatol": 1e-7})
Dq = r.x
Fq, CTq, _, _ = U4M(Dq)
say(f"undeformed U(4)+M: R_q* = {Dq:.7f}  F* = {Fq:.8f}  C_T = {CTq:.6f}  R[M]={2-2*Dq:.5f}  "
    f"R[Vhat_pm] = 8(1-R_q)-3 = {8*(1-Dq)-3:.6f} -> " +
    ('BELOW UNITARITY (Vhat free; naive F-max invalid)' if 8*(1-Dq)-3 < 0.5 else ('RELEVANT' if 8*(1-Dq)-3 < 2 else 'IRRELEVANT')))
say(f"   F-theorem along U(4)+M --(Vhat)--> candidate: F_UV={Fq:.6f} > F_IR={F_c:.6f}: {Fq > F_c}")

# corrected F-maximisation of U(4)+M: decouple Vhat_pm (and flux (1,-1,0,0)) when below unitarity
say("  F(R_q) scan of undeformed U(4)+M (uncorrected | V-hat-decoupling corrected):")
def U4M_corr(D):
    F = U4M(D, L=8.0, h=0.02)[0]
    RV = 8 * (1 - D) - 3
    Fc = F
    if RV < 0.5:
        Fc += 2 * (-chiral_logZ_b1(0.5)) - 2 * (-chiral_logZ_b1(RV))
    return F, Fc, RV, 16 * (1 - D) - 6
for D in [0.375, 0.45, 0.5, 0.54, 0.56, 0.58, 0.60, 0.62]:
    F, Fc, RV, R11 = U4M_corr(D)
    say(f"   R_q={D:.3f}: F={F:.5f}  F_corr={Fc:.5f}  R[Vhat]={RV:.3f}  R[flux(1,-1,0,0)]={R11:.3f}")

# BBP at N_f = N_c + 2 (the IR of the generic Higgs-branch point on the dual side):
#   U(4), 6 flavours, W = Vhat_+ + Vhat_-  <->  36 free singlets (U(0)); R_q = 1/6, R_M = 1/3
g6 = Grid(1 / 6, L=8.0, h=0.01)
lz6, _, _, _ = UN_logZ_and_d2(g6, 4, 6)
say(f"\nBBP U(4),N_f=6 at R_q=1/6: log Z = {lz6:.12f} ; 36*l(1-1/3) = {36*chiral_logZ_b1(1/3):.12f}")

# ---------------- spectrum of low-lying chiral operators
say("\n== chiral operators (R-charges at the candidate point) ==")
RQ = 5 / 8
for name, R in [("meson Q Qt (64, adj+1 of SU(8)^2)", 2 * RQ),
                ("monopoles V_pm, flux (+-1,0) [in W]", 8 * (1 - RQ) - 1),
                ("monopole flux (1,-1)", 8 * (1 - RQ) * 2 - 2),
                ("monopole flux (1,1)", 8 * (1 - RQ) * 2),
                ("dual: M", 2 * RQ), ("dual: q qt (F-term zero)", 2 * (1 - RQ)),
                ("dual monopoles Vhat_pm [in W]", 8 * (1 - (1 - RQ)) - 3)]:
    say(f"  {name:45s} R = {R:.4f}   {'ok' if R >= 0.5 else 'VIOLATES unitarity'}")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "candidate_u2_output.txt"), "w") as f:
    f.write("\n".join(out) + "\n")
