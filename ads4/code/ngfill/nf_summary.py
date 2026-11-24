"""
nf_summary (NEARGAP_FILL.md §0-§2): the table of the J <= 16 fill from cache/nf_state_J<J>.json (+ the J = 23 / 27 / 31 table of
nf_reprsum).  Usage: python3 nf_summary.py -> cache/nf_summary.txt
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_guard  # noqa: F401
C = os.path.join(HERE, 'cache')
L = ["| J | m/Λ | published model | K₂ own row (L sign-free / λ-level) | **direct own row** λ ± jk | budget: low / edge / trunc / interp / mc | lower bound | **L** (NU_LOW = 25) | certified-floor route |",
     "|---|---|---|---|---|---|---|---|---|"]
for J in (9, 10, 11, 13, 14, 15):
    f = os.path.join(C, 'nf_state_J%d.json' % J)
    if not os.path.exists(f):
        L.append(f"| {J} | (not run) |"); continue
    for o in json.load(open(f)):
        k, d = o['k2'], o['direct']
        b = d['base'][list(d['base'])[0]]['rig'] if J % 2 else min(v['rig'] for v in d['base'].values())
        L.append(f"| {J} | {o['mfl']:.4f} | {k['lam_model']:+.4f} | {k['lam_K2']:+.4f} ({k['L']:.2f} / {k['L_lam']:.2f}) | **{d['lam_central']:+.4f}** ± {o['jk_sigma']:.4f} | "
                 f"{d['E_low']:.3f} / {d['E_edge']:.3f} / {d['E_trunc']:.3f} / {d['E_interp']:.4f} / {o['E_mc']:.4f} | {o['lower']:+.3f} | **{o['L']:.2f}** ({o['L25']:.2f}) | {b:+.3f} |")
open(os.path.join(C, 'nf_summary.txt'), 'w').write("\n".join(L) + "\n")
print("\n".join(L))

# held own row: odd J = 11, 13, 15 have no K2 row ABOVE the own row (published rows 400/800, 336/400), so states with m > m_row are held at
# the own row (ng_bind's 1/m model holds outside the row range).  Estimate of that hold (NOT in the budgets above): the 1/m step from the
# own row to the state, as a fraction f of the step from the next-lower row to the own row, times the lambda shift that step produced
# (direct own row vs the published chain, which holds the lower row), x3.
import math
def mph(D, J):
    a = D - J - 2.0; return math.sqrt(a*(a + 2*J + 1))
Lm = math.sqrt(998.0*999.0)
H2 = ["", "| J | m/Λ | f = (1/m_row − 1/m)/(1/m_low − 1/m_row) | λ(direct own row) − λ(published, held at the lower row) | E_hold = 3 f |Δλ| | L with E_hold |", "|---|---|---|---|---|---|"]
for J, Dlow in ((11, 800.0), (13, 800.0), (15, 400.0)):
    f_ = os.path.join(C, 'nf_state_J%d.json' % J)
    if not os.path.exists(f_):
        continue
    mrow, mlow = mph(1000.0 + J, J), mph(Dlow, J)
    for o in json.load(open(f_)):
        m = o['mfl']*Lm
        if m <= mrow*(1 + 1e-9):
            continue
        fr = (1/mrow - 1/m)/(1/mlow - 1/mrow); dl = o['direct']['lam_central'] - o['k2']['lam_model']; Eh = 3*fr*abs(dl)
        H2.append(f"| {J} | {o['mfl']:.4f} | {fr:.3f} | {dl:+.4f} | {Eh:.4f} | {o['direct']['lam_central']/(o['E_tot'] + Eh):.2f} |")
# K2 calibration on the own rows (coarse 10-node K2 rows; J = 15: NEARGAP2's 36-node row): actual K2 error vs its last Laplace step
H3 = ["", "| J | λ_direct − λ_K2 (own row, twist edge) | λ_K2 − λ_K1 (last Laplace step) | ratio |", "|---|---|---|---|"]
for J in (9, 10, 11, 13, 14, 15):
    f_ = os.path.join(C, 'nf_state_J%d.json' % J)
    if not os.path.exists(f_):
        continue
    o = json.load(open(f_))[0]; k = o['k2']
    a, s = o['direct']['lam_K2'] - k['lam_K2'], k['lam_K2'] - k['lam_K1']
    H3.append(f"| {J} | {a:+.4f} | {s:+.4f} | {abs(a/s)*100:.0f} % |")
open(os.path.join(C, 'nf_summary.txt'), 'a').write("\n".join(H2 + H3) + "\n")
print("\n".join(H2 + H3))
