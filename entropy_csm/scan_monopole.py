#!/usr/bin/env python3
"""
Scan of U(N_c), N_f SQCD with monopole superpotential (BBP 1703.08460) for an existence-safe violation of
c_S/C_T <= 16 pi zeta3 via the Higgs-branch bound  c_S >= 3.5 (2 N_c N_f - N_c^2)  [real-scalar units].

Columns: R_Q* of plain SQCD T_0 (F-maximised), R[V] there, relevance, R_Q at T_M, C_T(T_M) (with BBP's
accidental-symmetry corrections), threshold 3.5*dimH, ratio C_T/(4 dimH).  Violation iff ratio < 7/8.
Also W = V_+ only (U(1)_A x U(1)_T -> diagonal; one mixing parameter F-maximised).
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import minimize_scalar
from u2_ct import Grid, UN_logZ_and_d2, CT_from_d2logZ, chiral_d2, chiral_logZ_b1, PI

out = []


def say(s):
    print(s, flush=True)
    out.append(s)


def Lfor(decay):           # integration half-width so that exp(-2 pi decay L) < e^-45
    return min(60.0, 45 / (2 * PI * max(decay, 0.05)) + 2.0)


def FZ(Nc, Nf, D, eta=0.0, need_p2=False, h=0.02):
    R0 = Nf * (1 - D) - (Nc - 1)
    g = Grid(D, L=Lfor(R0 - abs(eta)), h=h, need_p2=need_p2)
    lz, d1, d2, sg = UN_logZ_and_d2(g, Nc, Nf, eta)
    return -lz, (CT_from_d2logZ(d2) if need_p2 else None), sg


CTch = lambda R: CT_from_d2logZ(chiral_d2(R))
Fch = lambda R: -chiral_logZ_b1(R)


def T0(Nc, Nf):
    """F-maximise plain SQCD; returns R*, F*, R[V]; handles decoupled free mesons/monopoles."""
    lo = max(0.02, 1 - (Nc - 1 + 0.02) / Nf) if False else 0.05
    hi = 1 - (Nc - 1 + 0.03) / Nf          # keep R[V] > 0.03 for convergence
    f = lambda D: -FZ(Nc, Nf, D)[0]
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
    D = r.x
    RV = Nf * (1 - D) - (Nc - 1)
    flags = []
    if 2 * D < 0.5:
        flags.append("meson<1/2")
    if RV < 0.5:
        flags.append("V<1/2")
    if flags:   # decoupling correction: F_corr = F - sum F_ch(R_O) + sum F_ch(1/2)
        def fc(D):
            F = FZ(Nc, Nf, D)[0]
            RV = Nf * (1 - D) - (Nc - 1)
            if 2 * D < 0.5:
                F += Nf * Nf * (Fch(0.5) - Fch(2 * D))
            if RV < 0.5:
                F += 2 * (Fch(0.5) - Fch(RV))
            return F
        r = minimize_scalar(lambda D: -fc(D), bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
        D = r.x
        RV = Nf * (1 - D) - (Nc - 1)
        flags.append("corrected")
    return D, -r.fun, RV, flags


def TM(Nc, Nf):
    """C_T and F at the W=V_+ + V_- point, with BBP's IR identifications."""
    D = 1 - (Nc + 1) / Nf
    dimH = 2 * Nc * Nf - Nc * Nc
    note = ""
    if Nf == Nc + 2:
        if Nc == 1:
            CT = 9 * CTch(2 / 3); F = 9 * Fch(2 / 3); note = "WZ W=det M (9 chirals, R=2/3)"
        else:
            CT = 4 * Nf * Nf; F = Nf * Nf * Fch(0.5); note = "free N_f^2 mesons (det M irrelevant)"
        return D, CT, F, dimH, note
    F, CT, sg = FZ(Nc, Nf, D, need_p2=True, h=0.01)
    if 2 * D < 0.5:   # N_c+3 <= N_f <= 4(N_c+1)/3 : M decouple as free fields
        CT = CT - Nf * Nf * CTch(2 * D) + 4 * Nf * Nf
        F = F - Nf * Nf * Fch(2 * D) + Nf * Nf * Fch(0.5)
        note = "M free (accidental)"
    return D, CT, F, dimH, note


def onesided(Nc, Nf):
    """W = V_+ : eta = R0 - 2, F-maximise over R_Q with R[V_-] = 2 R0 - 2 > 0."""
    def pars(D):
        R0 = Nf * (1 - D) - (Nc - 1)
        return R0 - 2.0, 2 * R0 - 2
    # need R[V_-] = 2R0-2 > 0.03  ->  R0 > 1.015
    hi = 1 - (Nc - 1 + 1.02) / Nf
    lo = 0.05
    if hi <= lo:
        return None
    f = lambda D: -FZ(Nc, Nf, D, eta=pars(D)[0])[0]
    r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-6})
    D = r.x
    eta, RVm = pars(D)
    F, CT, sg = FZ(Nc, Nf, D, eta=eta, need_p2=True, h=0.01)
    edge = (abs(D - hi) < 1e-4) or (abs(D - lo) < 1e-4)
    return D, eta, RVm, F, CT, edge


if __name__ == "__main__":
    # --- calibration of the eta (topological-mixing) machinery: BBP one-monopole identity at N_f = N_c + 1
    say("== calibration: U(1), N_f=2, W=V_+  <->  4 mesons (R=2R_Q) + singlet S (R=2R0-2), no gauge group ==")
    for D in [0.3, 0.4]:
        R0 = 2 * (1 - D)
        eta = R0 - 2
        F, CT, _ = FZ(1, 2, D, eta=eta, need_p2=True, h=0.01)
        Fd = 4 * Fch(2 * D) + Fch(2 * R0 - 2)
        CTd = 4 * CTch(2 * D) + CTch(2 * R0 - 2)
        say(f"  R_Q={D}: F_el={F:.10f} F_WZ={Fd:.10f} | C_T el={CT:.8f} WZ={CTd:.8f}")
    say("== calibration: U(2), N_f=4, W=V_+  <->  U(1) N_f=4 dual not used; U(2) N_f=8 two-sided C_T=92.0149439 ==")
    D = 5 / 8
    say(f"  U(2) N_f=8: C_T = {FZ(2, 8, D, need_p2=True, h=0.01)[1]:.8f}")

    say("\n== two-sided W = V_+ + V_- ==")
    say(f"{'Nc':>3} {'Nf':>3} | {'R_Q*(T0)':>9} {'R[V]T0':>7} {'relev':>5} {'flags':>12} | {'R_Q(TM)':>8} "
        f"{'C_T(TM)':>9} {'7/8thr':>8} {'ratio':>7} {'F(TM)':>8} {'F28-type':>8} | note")
    rows = []
    for Nc in [1, 2, 3, 4, 5]:
        Nfmax = max(3 * Nc + 6, Nc + 6)
        for Nf in range(Nc + 2, Nfmax + 1):
            t = time.time()
            Ds, Fs, RV, flags = T0(Nc, Nf)
            D, CT, F, dimH, note = TM(Nc, Nf)
            thr = 3.5 * dimH
            ratio = CT / (4 * dimH)
            unit = (2 * D >= 0.5) or ("free" in note) or ("WZ" in note)
            rel = RV < 2
            viol = ratio < 7 / 8
            bbpchain = Nf <= 3 * Nc + 3
            rows.append((Nc, Nf, Ds, RV, rel, D, CT, thr, ratio, F, dimH * Fch(0.5), viol, bbpchain, note))
            say(f"{Nc:3d} {Nf:3d} | {Ds:9.5f} {RV:7.3f} {('yes' if rel else 'no'):>5} {','.join(flags) or '-':>12} | "
                f"{D:8.4f} {CT:9.3f} {thr:8.1f} {ratio:7.4f} {F:8.4f} {dimH*Fch(0.5):8.4f} | "
                f"{'VIOL ' if viol else ''}{'(i)ok ' if rel else ''}{'BBPchain ' if bbpchain else ''}{note}  [{time.time()-t:.0f}s]")
    say("\npoints passing (i) relevance AND (iii) violation:")
    hits = [r for r in rows if r[4] and r[11]]
    say(str([(r[0], r[1]) for r in hits]) if hits else "  none")
    say("points with violation (iii) and N_f <= 3N_c+3 (BBP Ising-chain existence, weaker):")
    say(str([(r[0], r[1], round(r[8], 4)) for r in rows if r[11] and r[12]]))

    say("\n== one-sided W = V_+ (R_Q, eta F-maximised along R[V_+]=2) ==")
    say(f"{'Nc':>3} {'Nf':>3} | {'R_Q':>8} {'eta':>8} {'R[V_-]':>7} {'mesonR':>7} {'C_T':>9} {'thr':>8} {'ratio':>7} {'F':>8} {'edge':>5}")
    for Nc in [1, 2, 3, 4]:
        for Nf in range(Nc + 2, 3 * Nc + 5):
            res = onesided(Nc, Nf)
            if res is None:
                continue
            D, eta, RVm, F, CT, edge = res
            dimH = 2 * Nc * Nf - Nc * Nc
            ratio = CT / (4 * dimH)
            say(f"{Nc:3d} {Nf:3d} | {D:8.4f} {eta:8.4f} {RVm:7.3f} {2*D:7.3f} {CT:9.3f} {3.5*dimH:8.1f} {ratio:7.4f} {F:8.4f} {str(edge):>5}"
                f" {'VIOL' if ratio < 7/8 else ''} {'unitarity!' if (2*D < 0.5 or RVm < 0.5) else ''}")

    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "scan_monopole_output.txt"), "w") as f:
        f.write("\n".join(out) + "\n")
