"""Contribution of large-r grid points (where float64 Taylor mode loses all digits, see ref_spot) to the moments and fits."""
import os, sys
import numpy as np, mpmath as mp
sys.path.insert(0, '/home/user/MAJOR/ads4/code/fourpt/tttt')
import t13_analysis_l2 as A13
from t13_analysis_l2 import Modes, CLASSES, QS, Tpred, lam
from tmoments_l2 import COMBOS
from moments import frame_data, series_converged
fn = sys.argv[1]; k, nu, cT = 6, 0.0, 8.0
mo = Modes(fn); D, J = mo.Delta, mo.J
g = mo.g; ok = g['ok']
fd = frame_data(g['U'][ok], g['S'][ok])
R = np.full(g['S'].shape, np.nan); R[ok] = fd['rs']
conv0 = mo.conv.copy()
print(f"# {os.path.basename(fn)}; converged points {conv0.sum()}, r range {np.nanmin(R[conv0]):.3f}..{np.nanmax(R[conv0]):.4f}; "
      f"#(r>0.8) {(R[conv0]>0.8).sum()}, #(r>0.9) {(R[conv0]>0.9).sum()}, #(r>0.95) {(R[conv0]>0.95).sum()}")
A = A13.design(D, J, k, nu)
lo, ls = float(lam('opp', mp.mpf(D), J)), float(lam('same', mp.mpf(D), J)); w = np.sqrt(2*ls/lo)
waves = {c: c for c in COMBOS}; waves['mix'] = {('same','same'):1.0, ('same','opp'):w, ('opp','same'):w, ('opp','opp'):w*w}
def run(mask, tag):
    mo.conv = conv0 & mask; mo._K = {}
    logref = None; out = []
    for name, cmb in waves.items():
        meas = {}
        for q in QS:
            v, sh = mo.moment(k, nu, cmb, q)
            if logref is None: logref = sh
            meas[q] = (-1)**J*complex(v)*np.exp(sh - logref)
        Tp = ({c: sum(wt*Tpred(a, b, D, J, cT)[c] for (a, b), wt in cmb.items()) for c in CLASSES} if name == 'mix' else Tpred(*name, D, J, cT))
        sol, resid, s1, r1 = A13.fit(A, meas, Tp, logref)
        ex = ''
        if name == 'mix':
            oo, ss = (1,1,-1,-1), (1,-1,-1,1); ex = f" Sig- {((sol[oo]-sol[ss])/(Tp[oo]-Tp[ss])).real:.4f}"
        out.append(f"{str(name)}: {s1.real:.4f} ({r1:.3f}){ex}")
    print(f"  {tag:28s}", " | ".join(out))
run(np.ones_like(conv0), 'all points (baseline)')
for rc in (0.95, 0.9, 0.85, 0.8):
    run(~(R > rc), f'drop r>{rc}')
# per-point weight of the dominant ss mode by r band
mo.conv = conv0; mo._K = {}
q = (1, -1, -1, 1)
U, S = g['U'], g['S']; xi = -np.log(g['us'])
F = A13.kernel_np(nu, [(0, 0)], xi)[(0, 0)]
rho, rhob = S**2, (S*U)**2
with np.errstate(divide='ignore', invalid='ignore'):
    logw = np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k - 4)*np.log(rho*rhob) + mo.logscale
logw[~(ok & conv0)] = -np.inf
for c in COMBOS:
    terms = np.where(np.isfinite(logw), np.exp(logw - np.max(logw)), 0.0)*F[:, None]*mo.modes[c][q]*g['W']
    tot = abs(terms.sum())
    print(f"  {c} mode {q}: sum|terms| r>0.95: {np.abs(terms[R>0.95]).sum()/tot:.1e}, r>0.9: {np.abs(terms[R>0.9]).sum()/tot:.1e}, r>0.8: {np.abs(terms[R>0.8]).sum()/tot:.1e}")
