"""Result Q, Stage 6 (RESULT_Q_STAGE6.md §2.4): how far in t the stored de-shift weights satisfy W(t) = Phi(t) (direct quadrature).
The heavy image at beta is the Abel transform of W on [0, beta]; a double-precision residual eps of the cosh(jt)-normalised rows becomes
|W - Phi| ~ eps cosh(j t), so the pure image is realised only up to a finite beta.  Prints |W - Phi| at t = 4 ... 16 for the reference gauges
(C4 min-norm value/derivative, dnu = 0.02; C5 min-norm, dnu = 0.04) and the SB gauge.  Output: cache/stage6_wcheck.txt"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import stage4_deshift as DS4, stage5_deshift as DS3
G = pickle.load(open(os.path.join(HERE, 'cache', 'stage6_gauge.pkl'), 'rb'))
C = pickle.load(open(os.path.join(HERE, 'cache', 'stage6_canon.pkl'), 'rb'))
T = np.array([4., 6., 8., 10., 12., 14., 16.])
lines = []
S4 = DS4.Deshift4('EE3', 0.02, 2048)
for lab, key, Phi in (('C4 value', 'C4v', lambda t: np.ones_like(t)), ('C4 deriv', 'C4d', lambda t: -t**2/2)):
    for g, chi in (('H0 min-norm', G[(key, 0.02)]['H0 min-norm']['chi']), ('SB', C[key]['chi'])):
        W = S4.W(chi, T); lines.append(f"{lab:8s} {g:11s}: |W - Phi|/|Phi| at t = {T.tolist()}: " + " ".join(f"{v:.1e}" for v in np.abs(W/Phi(T) - 1)))
S5 = DS3.Deshift3S('C5_pp', 0.04, 1024, general=True)
for g, c in (('G0 min-norm', G[('C5', 0.04)]['G0 min-norm']['c']), ('SB', C['C5']['c'])):
    W = S5.W(c, T); lines.append(f"C5       {g:11s}: |W - 1| at t = {T.tolist()}: " + " ".join(f"{v:.1e}" for v in np.abs(W - 1)))
print("\n".join(lines)); open(os.path.join(HERE, 'cache', 'stage6_wcheck.txt'), 'w').write("\n".join(lines) + "\n")
