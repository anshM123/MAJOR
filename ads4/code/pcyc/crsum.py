"""crsum: per-thimble diagnostics of a pcyc_crys run: <Q> per thimble (jackknife), |Y+|+|Y-| vs |Y+ + Y-| (cancellation), relative errors."""
import os, sys, pickle
import numpy as np
r = pickle.load(open(sys.argv[1], 'rb'))
R = r['R']
def jk(a, b):
    v = np.array([(a.sum() - a[i])/(b.sum() - b[i]) for i in range(len(a))]); return a.sum()/b.sum(), np.sqrt((len(a) - 1)/len(a)*np.sum(np.abs(v - v.mean())**2))
l0 = r['logN0']; Ys = []
for k, t in r['th'].items():
    q, e = jk(t['S1'], t['S0'])
    Y = np.mean(t['S0'])*np.exp(t['logN'] - l0); Ys.append(Y)
    print(f"{k}: <Q> = {q:.5e} +- {e:.1e} (rel {e/abs(q):.1e});  S0 rel err {np.std(t['S0'], ddof=1)/np.sqrt(R)/abs(np.mean(t['S0'])):.1e};  Y/Y0 = {Y:.4e}; diag {t['diag']}")
print(f"combined: <Q> = {r['Qcrys']:.6e} +- {r['Qerr']:.1e} (rel {r['Qerr']/abs(r['Qcrys']):.1e});  cancellation sum|Y|/|sum Y| = {sum(abs(y) for y in Ys)/abs(sum(Ys)):.2f}; time {r['time']:.0f}s")
