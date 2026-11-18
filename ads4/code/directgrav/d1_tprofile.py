"""
d1: where in the (u, t) moment grid do the exact <TTTT> kernel moments come from?  Cumulative kernel moment over t-columns (t ascending),
for W[X_s,X_s] and W[X_o,X_o], kernel O (odd) and K* (even), at nu = 0; also the moment with the smallest-t columns removed.
Diagnoses the small-t (near-lightcone / dropped-region) contamination of the low-k moments.
Usage: python3 d1_tprofile.py <modes.npz>
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TT = os.path.join(HERE, '..', 'fourpt', 'tttt')
for p in (TT, os.path.join(TT, '..'), os.path.join(TT, '..', 'oddk'), os.path.join(TT, '..', '..', 'e2e')):
    sys.path.insert(0, p)
import numpy as np
from t13_analysis_l2 import Modes, design, Tpred, CLASSES
from moments4 import kernel_np
sys.argv, fn = sys.argv[:1], sys.argv[1]
from a1_moments_ker import KER
mo = Modes(fn); g = mo.g; U, S = g['U'], g['S']; xi = -np.log(g['us'])
D, J = mo.Delta, mo.J
ts = g['ts']
def density(k, nu, wave, kern):
    K = KER[kern]; tot = 0; sh0 = None; pred = 0j
    A = design(D, J, k, nu); Tp = Tpred(wave[0], wave[1], D, J, 8.0)
    rho, rhob = S**2, (S*U)**2
    with np.errstate(divide='ignore', invalid='ignore'):
        logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + mo.logscale
    logw[~(g['ok'] & mo.conv)] = -np.inf
    sh = np.max(logw); base = np.where(np.isfinite(logw), np.exp(logw - sh), 0.0)*g['W']
    for (q1, q3), w13 in K.items():
        for (q2, q4), w24 in K.items():
            q = (q1, q3, q2, q4)
            F = kernel_np(nu, [(q1 + q3, q2 + q4)], xi)[(q1 + q3, q2 + q4)]
            tot = tot + w13*w24*base*F[:, None]*mo.modes[wave][q]
            pred += w13*w24*sum(A[q][c]*Tp[c] for c in CLASSES)
    return tot*(-1)**J, pred, sh
for k in (2, 4, 5, 6):
    for wave in (('same', 'same'), ('opp', 'opp')):
        for kern in ('O', 'KS'):
            dens, pred, sh = density(k, 0.0, wave, kern)
            col = dens.sum(axis=0)              # per t column
            R = lambda c: (c*np.exp(sh)/pred).real
            full = R(col.sum())
            cuts = [R(col[ts > t0].sum()) for t0 in (0.5, 1.0, 2.0, 4.0, 6.0)]
            # first t where the conv mask is complete (all u converged)
            print(f"k={k} {wave[0][0]}{wave[1][0]} {kern:2s}: R(full) = {full:+.4e} | R(t>0.5,1,2,4,6) = " + " ".join(f"{c:+.4e}" for c in cuts), flush=True)
