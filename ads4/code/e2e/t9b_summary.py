"""Compact markdown tables for SPINNING_E2E.md sect. 5: measured/predicted for the helicity-preserving and helicity-flip components
(two-mode extraction of t7_corrections.py) versus m, for l = 1, 2, k = 4, 6, nu = 0, 1, 3, beta ~ 0.4, 1, 2.
Entries whose estimated numerical error (dropped small-u region) exceeds 1e-4 are marked with *; those above 1e-3 with **."""
import os, numpy as np, mpmath as mp
from t4_spinning_e2e import main_cases, get_moments, moments_of, fname, DPHI
from t7_corrections import measured_components
from cases import info
mp.mp.dps = 25
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for l in (1, 2):
    data = {}
    for c in main_cases():
        if not os.path.exists(fname(l, *c)): continue
        m, b_, eta = info(*c)
        mom = get_moments(l, *c)
        res = moments_of(mom, (2, 4, 6), (0.0, 1.0, 3.0))
        for k in (2, 4, 6):
            for nu in (0.0, 1.0, 3.0):
                Rp, Rf, worst, bb = measured_components(l, c[0], c[1], DPHI, k, nu, res)
                tail = max(mom.tail_estimate(k, nu, l, -l), mom.tail_estimate(k, nu, l, l))
                data[(round(float(b_), 1), k, nu, c)] = (float(m), float(mp.re(Rp)), float(mp.re(Rf)), tail)
    for k in (6, 4, 2):
        pr(f"### l = {l}, k = {k}:  R_pres / R_flip  (measured/predicted)")
        for beta in (0.4, 1.0, 2.0):
            cs = sorted([key[3] for key in data if key[0] == beta and key[1] == k and key[2] == 0.0], key=lambda c: c[0])
            pr("| β ≈ %.1f, (Δ,J) → | " % beta + " | ".join(f"({c[0]:.2f}, {c[1]})" for c in cs) + " |")
            pr("|---|" + "---|"*len(cs))
            pr("| m | " + " | ".join(f"{data[(beta, k, 0.0, c)][0]:.1f}" for c in cs) + " |")
            for nu in (0.0, 1.0, 3.0):
                row = []
                for c in cs:
                    m, Rp, Rf, tail = data[(beta, k, nu, c)]
                    mark = "**" if tail > 1e-3 else ("*" if tail > 1e-4 else "")
                    row.append(f"{Rp:.4f} / {Rf:.4f}{mark}")
                pr(f"| ν = {nu:.0f} | " + " | ".join(row) + " |")
            pr("")
open('out_t9b_summary.txt', 'w').write("\n".join(out) + "\n")

# --- sensitivity controls: what the ratio measured/predicted would be with plausible wrong ingredients (largest m of each beta, k = 6, nu = 0)
from formulas import log_b
from predict import prediction, components
pr2 = []
pr2.append("### controls (k = 6, nu = 0, mode (l,-l)): ratio measured/predicted with the derived formula, and with one ingredient altered")
pr2.append("| l | (Δ,J) | m | derived | b^(0) instead of b^(1−Δφ) | (Δ,J) instead of (Δ−3/2, J+1/2) in T_l | M_pres only (flip term dropped) |")
pr2.append("|---|---|---|---|---|---|---|")
for l in (1, 2):
    for c in main_cases():
        if not os.path.exists(fname(l, *c)): continue
        m, b_, eta = info(*c)
        if float(m) < 60: continue
        mom = get_moments(l, *c)
        res = moments_of(mom, (6,), (0.0,))
        ex = res[(6, 0.0, l, -l)]; pd = prediction(l, c[0], c[1], DPHI, 6, 0.0, l, -l)
        r = mp.re(ex/pd)
        rb0 = r*mp.exp(log_b(c[0], c[1], 0.0) - log_b(c[0], c[1], 1 - DPHI))
        Ap, Am = mp.mpf(c[0]) - 1.5, c[1] + 0.5
        rDJ = r*((Ap*Am)/(mp.mpf(c[0])*c[1]))**(2*l)
        A, B = components(l, 6, 0.0, l, -l)
        from formulas import Fnu_mp
        X0 = Fnu_mp(0.0, 0, b_); Y0 = Fnu_mp(0.0, 2*l, b_)
        rP = r*(A*X0 + B*Y0)/(A*X0)
        pr2.append(f"| {l} | ({c[0]:.2f}, {c[1]}) | {float(m):.1f} | {float(r):.5f} | {float(rb0):.5f} | {float(rDJ):.5f} | {float(mp.re(rP)):.5f} |")
for s in pr2: print(s)
open('out_t9b_summary.txt', 'a').write("\n".join(pr2) + "\n")
