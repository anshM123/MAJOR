"""
m1 (M6, scalar calibration, empirical): how the EXACT scalar-skeleton remainder rho_k(x; m, J) (exactskel e1, Mellin/Barnes, Delta_phi = 3)
scales with m at fixed x = p/(m/Lam) and fixed (m/Lam, J), across Delta_g = 100, 300, 1000 (the 64 p-nodes are shared, so x is identical).
For each (m/Lam, J) present at all three gaps and nodes with x in [0.3, 0.7]: m rho, m^2 rho, and the local exponent
alpha = -dlog|rho|/dlog m between successive gaps.  rho is defined by  eps_k m^{2k} C_k = e_k D + 2 eps_k (P - J0) + env_k rho_k  (e3/e7 model),
i.e. it is what remains after the bulk-point density and the Regge curvature term.  Also the size of env_k rho_k relative to the density.
Usage: python3 m1_scalar_scaling.py -> cache/m1_scalar_scaling.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
C = os.path.join(HERE, '..', 'exactskel', 'cache')
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
db = {Dg: pickle.load(open(os.path.join(C, 'e1_%d_Df3.pkl' % Dg), 'rb')) for Dg in (100, 300, 1000)}
common = sorted(set(db[100]) & set(db[300]) & set(db[1000]))
pr(f"# m1: exact scalar remainder rho_k at fixed x, (m/Lam, J) across Delta_g = 100/300/1000 (Df = 3); {len(common)} common states")
sel = [kk for kk in common if kk[1] in (0, 2, 4, 6) and kk[0] in (1.0, 1.01, 1.05, 1.1, 1.2, 1.3, 1.5)]
for (mfl, J) in sel:
    recs = [db[Dg][(mfl, J)] for Dg in (100, 300, 1000)]
    ms = np.array([r['m'] for r in recs]); x = recs[0]['x']
    win = np.where((x >= 0.3) & (x <= 0.7))[0]
    if len(win) == 0:
        continue
    for k in (2, 3, 4):
        rh = np.array([r['rho'][k][win] for r in recs])            # 3 x nwin
        a1 = -np.log(np.abs(rh[1]/rh[0]))/np.log(ms[1]/ms[0]); a2 = -np.log(np.abs(rh[2]/rh[1]))/np.log(ms[2]/ms[1])
        mr2 = rh*ms[:, None]**2
        pr(f"  m/Lam={mfl} J={J} k={k}: x in [{x[win[0]]:.2f},{x[win[-1]]:.2f}] ({len(win)} nodes); |rho| max at m~{ms[0]:.0f}/{ms[1]:.0f}/{ms[2]:.0f}: "
           f"{np.max(np.abs(rh[0])):.2e}/{np.max(np.abs(rh[1])):.2e}/{np.max(np.abs(rh[2])):.2e}; local exponent alpha (median, range) "
           f"{np.median(a1):.2f} [{a1.min():.2f},{a1.max():.2f}] / {np.median(a2):.2f} [{a2.min():.2f},{a2.max():.2f}]; "
           f"m^2 rho at x={x[win[len(win)//2]]:.2f}: {mr2[0, len(win)//2]:+.1f}, {mr2[1, len(win)//2]:+.1f}, {mr2[2, len(win)//2]:+.1f}")
open(os.path.join(HERE, 'cache', 'm1_scalar_scaling.txt'), 'w').write("\n".join(out) + "\n")
