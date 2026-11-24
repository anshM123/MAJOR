"""
gn_stepnorm (GAP_NEARGAP.md §2.4): split the chain-level K2 - K1 step of a row, Gamma_K2(x) - Gamma_K1(x) (c2), into a row
renormalisation component c D_h(x) (a change of the fitted row constant a between the two orders: Gamma = m^2 (G/(a env) - D_h)) and
the remainder; least squares over the row's nodes.  The renormalisation is positivity-neutral to first order (SMALLNU §4.2).
Usage: python3 gn_stepnorm.py J D[,D...]
"""
import os, sys
sys.dont_write_bytecode = True
import gn_bind as G
import numpy as np
J = int(sys.argv[1]); Ds = [float(v) for v in sys.argv[2].split(',')]
g2, _ = G.NB.grid_for(J, 2, False); g1, _ = G.NB.grid_for(J, 1, False)
h = '+-' if J % 2 else '++'
keys = ['C2', 'C3H'] if J % 2 else [('C2', '++'), ('C2', '+-'), ('C3H', '++'), ('C3H', '+-')]
for D in Ds:
    r2, r1 = g2[D], g1[D]
    xc = np.intersect1d(np.round(r2['xs'], 6), np.round(r1['xs'], 6))
    for key in keys:
        hh = key[1] if isinstance(key, tuple) else '+-'
        v2 = np.array([r2[key][np.argmin(abs(r2['xs'] - x))] for x in xc]); v1 = np.array([r1[key][np.argmin(abs(r1['xs'] - x))] for x in xc])
        st = v2 - v1; Dh = G.E.Dc_fixed(hh, J, xc)
        c = float(st @ Dh/(Dh @ Dh)); res = st - c*Dh
        sel = xc <= 0.8 + 1e-9
        print(f"J={J} row {D:g} {str(key):14s}: max|step| {np.abs(st[sel]).max():7.0f}  renorm coeff c = {c:+8.0f} c2  "
              f"max|step - c D_h| {np.abs(res[sel]).max():7.0f}  (x <= 0.8; rms {np.sqrt(np.mean(st[sel]**2)):6.0f} -> {np.sqrt(np.mean(res[sel]**2)):6.0f})")
