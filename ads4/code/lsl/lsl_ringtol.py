"""
lsl/lsl_ringtol.py (LARGESPIN_L.md §4.3): Lemma-I tolerance of ring-located (AG) corrections at large impact parameter.
A beyond-LR correction whose H2 density (Lemma I of LARGESPIN_ANALYTIC) is supported on t >= (1 - eps_r) beta with total variation w
(in units of the LR ring weight) changes the k = 2 action by at most w * 2 max_{t >= (1-eps_r) beta} |Omega_2(t)| (+ k = 3 likewise);
with the two-term form (Lemma LR-inf, valid with the stated remainder for t >= 1) Omega_2 is decreasing there, so the tolerated w is
   w_tol = F_LR(beta) / [2 Omega_2((1-eps_r) beta) + u |Omega_3((1-eps_r) beta)|].
Output: cache/ringtol.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import lsl_core as C
lines = []
for b in (1.2, 1.5, 2.0, 3.0, 4.0, 5.3, 7.6):
    row = f"beta = {b}:"
    for er in (0.03, 0.1, 0.3):
        bt = (1 - er)*b
        for u in (1.0, 0.59):
            F = C.LR_action(b, u, -1.0)
            den = 2*C.Omega_form(2, bt) + u*abs(C.Omega_form(3, bt)) if bt >= 1 else 2*C.Omega(2, bt) + u*abs(C.Omega(3, bt))
            row += f"  eps_r={er},u={u}: w_tol={F/den:.3f}"
    print(row, flush=True); lines.append(row)
mono = [C.M1(b) for b in np.linspace(0.05, 12, 120)]
lines.append(f"M1 decreasing on [0.05, 12] (120 samples): {bool(np.all(np.diff(mono) < 0))}")
print(lines[-1])
open(os.path.join(C.CACHE, 'ringtol.txt'), 'w').write(__doc__ + "\n" + "\n".join(lines) + "\n")


# tolerance profile at large impact parameter (Lemma I): weight tolerated at radius t, in units of the ring weight, worse helicity
def Om(k, t):
    return C.Omega(k, t) if t < 2 else C.Omega_form(k, t)
prof = ["\nTolerance profile tol(t) = F_LR(beta)/(2|Omega_2(t)| + u|Omega_3(t)|) at the gap (u from m_gap):"]
for J in (857, 1000, 3000, 10000, 100000, 1000000):
    m = C.mgap(1000, J); b = C.beta_of(J, m); u = (C.lam(1000)/m)**2
    F = C.LR_action(b, u, -1.0)
    row = f"J = {J:7d} (beta = {b:.3f}):"
    for t in (0.003, 0.03, 0.1, 0.3, 1.0, b/2, b - 1, b):
        if t <= 0:
            continue
        row += f"  t={t:.3g}: {F/(2*abs(Om(2, t)) + u*abs(Om(3, t))):.2e}"
    prof.append(row); print(row, flush=True)
open(os.path.join(C.CACHE, 'ringtol.txt'), 'a').write("\n".join(prof) + "\n")
