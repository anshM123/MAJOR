"""
dc_raw (DIRCV.md §2): raw-level comparison at (1100, 100) (no a_h, no fit, no kappa): on each node, per common-seed replica,
  scalar : d_w = log(S_waves/S_ref) - exact,  d_P = log(S_PJ/S_ref) - exact            (exact = z3 CV mean, fw_scalar.install)
  graviton (EE, C3H; 2x2 record components): g = log(G_waves/G_PJ) at the STATE (the reference cancels; same CRN)
  CV hypothesis: E[g] = E[log(S_waves/S_PJ)] =: s  (common relative bias);  DIR hypothesis: E[g] = 0.
Jackknife over the 8 replicas.  Usage: python3 dc_raw.py [J Delta]  -> cache/dc_raw.txt
"""
import dc_guard, os, sys, pickle
import numpy as np
import z3_nodes as Z3, fw_scalar as FS, pw_gamma as PG

J = int(sys.argv[1]) if len(sys.argv) > 1 else 100
D = float(sys.argv[2]) if len(sys.argv) > 2 else 1100.0
SC = os.path.join(dc_guard.SPIN, 'cache')
w = pickle.load(open(os.path.join(SC, 'pw_pts_J%d_D%g_w.pkl' % (J, D)), 'rb'))
P = pickle.load(open(os.path.join(SC, 'pw_pts_J%d_D%g_wP.pkl' % (J, D)), 'rb'))
FS.install(J, D)
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

def jk(f, R):
    v0 = f(np.ones(R, bool)); j = np.array([f(np.arange(R) != i) for i in range(R)])
    return v0, j, np.sqrt((R - 1)/R*np.sum(np.abs(j - j.mean(0))**2, axis=0))

m = Z3.mphys(D, J)
pr(f"# dc_raw J={J} Delta={D:g} m={m:.2f}")
for nu in sorted(set(k[1] for k in w if k[0] == D)):
    sw, sp, rf = PG.to_z2(w[(D, nu)], 'EE'), PG.to_z2(P[(D, nu)], 'EE'), PG.to_z2(w[(Z3.DREF, nu)], 'EE')
    rfP = PG.to_z2(P[(Z3.DREF, nu)], 'EE')
    same_ref = np.allclose(rf['repsS'], rfP['repsS'])
    R = len(sw['repsS'])
    ow = Z3.scalar_check(J, nu, {D: sw, Z3.DREF: rf}, D); op = Z3.scalar_check(J, nu, {D: sp, Z3.DREF: rf}, D)
    fac = ow['dc2']/ow['d'].real if ow['d'].real != 0 else np.nan
    # scalar waves vs PJ at the state (same CRN)
    s0, sj, se = jk(lambda s: np.log(sw['repsS'][s].mean()/sp['repsS'][s].mean()) + sw['logNs'] - sp['logNs'], R)
    pr(f"nu={nu:6.2f} x={nu/m:.4f} (same ref in _w/_wP: {same_ref}); c2-factor m^2 ex2/E2 = {fac:.3g}")
    pr(f"   scalar  d_waves = {ow['d'].real:+.5f} +- {ow['err'].real:.5f} (Dc2 {ow['dc2']:+.0f} +- {ow['ec2']:.0f});  "
       f"d_PJ = {op['d'].real:+.5f} +- {op['err'].real:.5f} (Dc2 {op['dc2']:+.0f} +- {op['ec2']:.0f});  "
       f"log(S_w/S_P) = {s0.real:+.5f} +- {se.real:.5f}  Im {s0.imag:+.5f} +- {se.imag:.5f}")
    for kind in ('EE', 'C3H'):
        gw, gp = PG.to_z2(w[(D, nu)], kind), PG.to_z2(P[(D, nu)], kind)
        def f(s):
            return np.log(gw['reps'][s].mean(0)/gp['reps'][s].mean(0)) + gw['logN'] - gp['logN']
        g0, gj, ge = jk(f, R)
        def fd(s):
            return (np.log(gw['reps'][s].mean(0)/gp['reps'][s].mean(0)) + gw['logN'] - gp['logN']
                    - (np.log(sw['repsS'][s].mean()/sp['repsS'][s].mean()) + sw['logNs'] - sp['logNs']))
        d0, dj, de = jk(fd, R)
        relw = np.abs(gw['err']/gw['E']); relp = np.abs(gp['err']/gp['E'])
        for i in range(2):
            for k in range(2):
                pr(f"   {kind:3s}[{i}{k}] log(G_w/G_P) = {g0[i,k].real:+.5f} +- {ge[i,k].real:.5f} (Im {g0[i,k].imag:+.5f} +- {ge[i,k].imag:.5f});"
                   f"  minus scalar: {d0[i,k].real:+.5f} +- {de[i,k].real:.5f};  rel err w/P {relw[i,k]:.4f}/{relp[i,k]:.4f}")
open(os.path.join(dc_guard.CACHE, 'dc_raw_J%d_D%g.txt' % (J, D)), 'w').write('\n'.join(lines) + '\n')
