#!/usr/bin/env python3
"""Scans of R = (c_S/C_T)/(16 pi zeta3) over the planar CS-vector families. Uses csm_ratio.py."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import brentq
import csm_ratio as c

PI, Z3 = np.pi, c.Z3
lams = np.concatenate([[1e-4, 1e-3, 1e-2], np.linspace(0.05, 0.95, 19), [0.99, 0.999]])


def line(*a):
    print(" ".join(f"{x:>12.6g}" if isinstance(x, (float, np.floating)) else f"{str(x):>12s}" for x in a))


print("=== RF / CB (bosonization pair), N=2 ===")
line("lam", "mu_RF", "R_RF", "R_CB(1-lam)", "mu_N2", "R_N2", "R_RFnaive")
for l in lams:
    line(l, c.mu_RF(l), c.R_RF(l), c.R_CB(1 - l), c.mu_N2(l), c.R_N2(l), c.R_RF_naive(l))

print("\n=== RB: stable domain lamhat in [2lam/(2-lam), 2]; also lamhat outside ===")
line("lam", "lh_min", "R(lh_min)", "R(x6D=0)", "R(lh=2)", "R(lh=1e-3)", "R(lh=10)")
for l in lams:
    lo, hi = c.RB_domain(l)
    line(l, lo, c.R_RB(l, lo), c.R_RB(l, 2 * l), c.R_RB(l, hi), c.R_RB(l, 1e-3), c.R_RB(l, 10.0))

# monotonicity of R_RB in lamhat at fixed lam, and dense 2D max over stable domain
print("\n=== RB dense scan of stable domain ===")
best = (-1, None)
viol = 0
for l in np.linspace(0.002, 0.998, 250):
    lo, hi = c.RB_domain(l)
    vals = [c.R_RB(l, h) for h in np.linspace(lo, hi, 40)]
    if np.any(np.diff(vals) > 1e-14):
        print("  non-monotone in lamhat at lam", l)
    if max(vals) > 1:
        viol += 1
    if max(vals) > best[0]:
        best = (max(vals), l)
print(f"  max R on sampled domain = {best[0]:.10f} at lam={best[1]:.4f} (sup -> 1 only as lam->0, lamhat->0); points with R>1: {viol}")

# critical fermion: ghat in (0, inf); dual-stable range ghat in [1-lam_f, 1+lam_f]
print("\n=== CF: R vs ghat (ghat=1 <-> lambda6^f=0; GN at lam_f->0) ===")
line("lam_f", "R(g=1-l)", "R(g=1)", "R(g=1+l)", "R(g=1e-3)", "R(g=1e3)")
for l in lams:
    line(l, c.R_CF(l, 1 - l), c.R_CF(l, 1.0), c.R_CF(l, 1 + l), c.R_CF(l, 1e-3), c.R_CF(l, 1e3))

# small-lambda asymptotics (the danger region)
print("\n=== small-lambda asymptotics ===")
line("lam", "1-R_RB(lo)", "pred", "1-R_RB(lh=lam)", "1-R_RF-1/4", "1-R_N2-1/8")
for l in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
    lo, _ = c.RB_domain(l)
    m = c.mu_RB(l, lo)
    phi = PI * l
    # leading-log prediction: 1-R ~ [ (mu^2+phi^2) log(1/sqrt(mu^2+phi^2)) ]/(6 zeta3)
    rho = np.hypot(m, phi)
    pred = (m * m + phi * phi) * np.log(1 / rho) / (6 * Z3)
    line(l, 1 - c.R_RB(l, lo), pred, 1 - c.R_RB(l, l), 1 - c.R_RF(l) - 0.25, 1 - c.R_N2(l) - 0.125)

print("\n=== N=1 'w' family (1305.7235), conformal point ===")


def solve_w(l, w, ngrid=300):
    """All real solutions (c_B>0, c_F>=0); returns list of (R, cB, cF)."""
    x4, x6 = (1 + w) / 2, w * w
    def cF_of(cB):
        sB = c.S_tad(cB, l)
        g = lambda x: x - abs(2 * l * (c.C_tad(x, l) - x4 * sB))
        xs = np.linspace(1e-9, 40, 200)
        gs = [g(x) for x in xs]
        roots = [brentq(g, xs[i], xs[i + 1]) for i in range(len(xs) - 1) if gs[i] * gs[i + 1] < 0]
        return roots
    sols = []
    cBs = np.linspace(1e-6, 25, ngrid)
    for branch in range(2):
        def h(cB):
            r = cF_of(cB)
            if len(r) <= branch:
                return np.nan
            cF = r[branch]
            C, S = c.C_tad(cF, l), c.S_tad(cB, l)
            return cB ** 2 - (l * l * (1 + 3 * x6) * S * S + 4 * x4 * l * l * (C * C - 2 * x4 * S * C))
        hs = np.array([h(x) for x in cBs])
        for i in range(len(cBs) - 1):
            if np.isfinite(hs[i]) and np.isfinite(hs[i + 1]) and hs[i] * hs[i + 1] < 0:
                cB = brentq(h, cBs[i], cBs[i + 1])
                cF = cF_of(cB)[branch]
                R = 0.5 * (c.G_poly(cB, PI * l) + c.G_poly(cF, PI * (1 - l))) / Z3
                sols.append((R, cB, cF))
    return sols


line("lam", "w", "nsol", "R_dominant", "c_B", "c_F")
wmax = (-1, None)
for l in [0.05, 0.25, 0.5, 0.75, 0.95]:
    for w in [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 4.0]:
        s = solve_w(l, w, ngrid=120)
        if s:
            Rd, cB, cF = max(s)          # lowest F = largest sum of G
            line(l, w, len(s), Rd, cB, cF)
            wmax = max(wmax, (Rd, (l, w)))
        else:
            line(l, w, 0, "none", "", "")
print("  max R over sampled w-family points:", wmax)

# bound for ANY boson+fermion planar theory: R <= [G(0,pi lam)+G(0,pi(1-lam))]/(2 zeta3)
print("\n=== envelope for any one-boson+one-fermion planar theory: [G(0,pl)+G(0,p(1-l))]/2z3 ===")
env = [(0.5 * (c.G_int(0, PI * l) + c.G_int(0, PI * (1 - l))) / Z3, l) for l in np.linspace(1e-4, 0.5, 60)]
print("  max", max(env), " (7/8 =", 7 / 8, ")")

# figure
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    L = np.linspace(0.002, 0.998, 200)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(L, [c.R_RF(l) for l in L], label="regular fermion = critical boson(1-λ)")
    ax.plot(L, [c.R_CB(l) for l in L], label="critical boson", ls="--")
    lo = [c.R_RB(l, c.RB_domain(l)[0]) for l in L]
    hi = [c.R_RB(l, 2.0) for l in L]
    ax.fill_between(L, hi, lo, alpha=0.3, label="regular boson, stable λ6 band")
    ax.plot(L, [c.R_N2(l) for l in L], label="N=2, one chiral")
    ax.plot(L[L < 0.9], [c.R_RF_naive(l) for l in L[L < 0.9]], ls=":", c="k", label="RF, naive zero holonomy (wrong)")
    ax.axhline(1, c="r", lw=0.8)
    ax.set_xlabel("λ = N/k (bosonic λ for RB/CB, fermionic for RF/CF)")
    ax.set_ylabel("R = (c_S/C_T)/(16π ζ(3))")
    ax.set_ylim(0.6, 1.25)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(os.path.dirname(os.path.abspath(__file__)), "R_vs_lambda.png"), dpi=140)
    print("figure written")
except Exception as e:
    print("no figure:", e)
