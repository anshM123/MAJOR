"""rv12_chi_envelope.py -- numerical check of the note's Sec.4 envelopes |chi R_+-| <= 2^90 e^{-(x-L)_+}, |chi| <= 2^90 (1+x)^-2 e^{-(x-L)_+},
and of the actual decay rate of the uncut chi beyond Lambda, on the rv5 double tables (Lambda^2 = 900, 997002)."""
import numpy as np, math
for lam2 in (900, 997002):
    d = np.load(f'out/chi_{lam2}.npz'); nu, chi = d['nu'], d['chi']; L = float(d['L'])
    z = nu.astype(complex); I = 1j
    Rm = -2*(2*z-I)*(2*z+I)**2*(2*z+3*I)*(2*z+5*I)*(2*z+11*I)/(z*(2*z-7*I)*(2*z-3*I)*(2*z+7*I))
    ex = np.maximum(nu-L, 0)
    sel = np.abs(chi) > 1e-280
    a = np.max(np.log2(np.abs(chi*Rm)[sel]) + ex[sel]/math.log(2))
    b = np.max(np.log2(np.abs(chi)[sel]*(1+nu[sel])**2) + ex[sel]/math.log(2))
    # decay rate beyond L: slope of log|chi| on [L+20, L+120]
    s2 = (nu > L+20) & (nu < L+120) & sel
    slope = np.polyfit(nu[s2], np.log(np.abs(chi[s2])), 1)[0]
    print(f"Lambda^2={lam2}: max log2(|chi R_-| e^(x-L)_+) = {a:.2f} (claim <= 90); max log2(|chi|(1+x)^2 e^(x-L)_+) = {b:.2f} (claim <= 90); "
          f"fitted decay rate d log|chi|/dx on [L+20,L+120] = {slope:.4f} (pi/2 = {math.pi/2:.4f})")
