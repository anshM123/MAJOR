"""
ef_j7 (EXT_FOLD.md §3.4): the J = 7 near-gap state at Delta_g = 1000 (odd, 1x1 parity-odd +-) with the external DIRECT own-mass row
(job A: 34 n1 nodes on Delta = 1007, x = 0.05 ... 0.775; no K2 row exists at 1007 -- the published chain holds J = 7 between its
K2 rows 800 and 1600 by 1/m interpolation).  Exactly highspin/hs_even.py's protocol for J = 12, 16 (neargap2/ng_bind.run_state,
read-only, with the direct row injected through hs_core / ng_val's loader wrappers; E_mc = 3 x jackknife sigma over the 8 replicas),
here for odd J.  The direct records need the parity-odd <TTO> shadow of (1007, 7), computed by ef_shadow.py (p2_shadow_odd, defaults;
it reproduces cached shadows bit for bit).
Budget (NEARGAP2 §1.2, x3 on model terms): E_trunc (own-row nodes direct: 0 at the twist edge; at m > m_row the 1/m partner row 1600
is K2), E_interp, E_low (nu < 50: 1e4 c2 sign-free with the gauge factor; no direct node below x = 0.05 here), E_edge (x > 0.775:
G_edge = max(1e4, 10 max|Gamma| of the last 4 own-row nodes), sign-free), E_mc.
Usage: python3 ef_j7.py [mfl,...] -> cache/ef_j7.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C
import numpy as np
H = C.H
V = H.V
import ng_bind as NB          # read-only import

if __name__ == '__main__':
    J, Down, Dg = 7, 1007.0, 1000.0
    mtw = H.HN.mphys(Down, J)/H.E.lam_of(Dg)
    mfls = [mtw] + ([float(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1.01, 1.02, 1.05])
    lines, res = [], []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    xd = H.direct_x(J, Down); bad = H.failed_x(J, Down)
    pr(f"# ef_j7 J={J} own row {Down:g} (m = {H.HN.mphys(Down, J):.3f}, twist edge m/Lam = {mtw:.7f}): direct nodes {len(xd)} "
       f"(x {xd.min():.3f} .. {xd.max():.3f}); failed: {sorted(bad)}; NU_LOW = {NB.NU_LOW}")
    for mfl in mfls:
        pr(f"## m/Lam = {mfl:.7f}")
        V.INJ.pop(J, None); V.JK[0] = None
        pr("#  published K2 rows only (200/400/800/1600; the state 1/m-interpolated between 800 and 1600):")
        o0 = NB.run_state(Dg, J, Down, mfl, pr)
        V.INJ[J] = [Down]
        pr("#  own-mass DIRECT row 1007 added (all replicas):")
        o1 = NB.run_state(Dg, J, Down, mfl, pr)
        lj = []
        for jk in range(8):
            V.JK[0] = jk
            lj.append(NB.run_state(Dg, J, Down, mfl, lambda s: None)['lam_central'])
        V.JK[0] = None; V.INJ.pop(J, None)
        lj = np.array(lj); sj = float(np.sqrt(7/8*((lj - lj.mean())**2).sum()))
        E_mc = 3*sj
        tot = o1['E_tot'] + E_mc; totl = o1['E_tot_lam'] + E_mc
        o = dict(J=J, mfl=mfl, Down=Down, k2=o0, direct=o1, jk=lj.tolist(), jk_sigma=sj, E_mc=E_mc, E_tot=tot, lower=o1['lam_central'] - tot,
                 L=o1['lam_central']/tot, E_tot_lam=totl, lower_lam=o1['lam_central'] - totl, L_lam=o1['lam_central']/totl)
        res.append(o)
        pr(f"  => J={J} m/Lam={mfl:.7f}: published K2 chain {o0['lam_K2']:+.5f} (its NEARGAP2-style budget: L {o0['L']:.2f} sign-free / "
           f"{o0['L_lam']:.2f} lambda-level) | own direct row {o1['lam_K2']:+.5f}, central {o1['lam_central']:+.5f} +- {sj:.5f} (jk) | "
           f"budget sign-free {tot:.3e} (trunc {o1['E_trunc']:.2e}, interp {o1['E_interp']:.2e}, low {o1['E_low']:.2e}, edge {o1['E_edge']:.2e} "
           f"[G_edge {o1['G_edge']:.0f}], mc {E_mc:.2e}) -> lower {o['lower']:+.4f}, L = {o['L']:.1f}; lambda-level {totl:.3e} -> L = {o['L_lam']:.1f}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C.EFC, 'ef_j7.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.EFC, 'ef_j7.json'), 'w'), indent=1, default=float)
