"""Tables for SPINNING_E2E.md: measured (exact block) and predicted (derived formula) values, separately, with their ratio,
for the modes (l,-l) [mostly helicity preserving] and (l,l) [mostly helicity flip], and the extracted component ratios."""
import os, sys, numpy as np, mpmath as mp
from t4_spinning_e2e import main_cases, get_moments, moments_of, fname, DPHI
from t7_corrections import measured_components
from predict import prediction
from cases import info
mp.mp.dps = 25
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
for l in (1, 2):
    pr(f"## l = {l}:  (-1)^(J+l) Pihat^(q1,q3)_(k,nu)[stripped W_diff], Dphi = {DPHI}")
    pr("| Delta | J | m | beta | k | nu | mode | measured | predicted | meas/pred | num. err. | R_pres | R_flip |")
    pr("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for c in main_cases():
        if not os.path.exists(fname(l, *c)): continue
        m, b_, eta = info(*c)
        mom = get_moments(l, *c)
        res = moments_of(mom, (4, 6), (0.0, 1.0, 3.0))
        sgn = (-1)**l
        for k in (4, 6):
            for nu in (0.0, 1.0, 3.0):
                Rp, Rf, worst, bb = measured_components(l, c[0], c[1], DPHI, k, nu, res)
                for (q1, q3) in ((l, -l), (l, l)):
                    ex = sgn*res[(k, nu, q1, q3)]; pd = sgn*prediction(l, c[0], c[1], DPHI, k, nu, q1, q3)
                    pr(f"| {c[0]:.4f} | {c[1]} | {float(m):.2f} | {float(b_):.3f} | {k} | {nu:.0f} | ({q1:+d},{q3:+d}) | {mp.nstr(mp.re(ex), 9)} | {mp.nstr(mp.re(pd), 9)} | {float(mp.re(ex/pd)):.6f} | {mom.tail_estimate(k, nu, q1, q3):.0e} | {float(mp.re(Rp)):.6f} | {float(mp.re(Rf)):.6f} |")
open('out_t9_tables.txt', 'w').write("\n".join(out) + "\n")
