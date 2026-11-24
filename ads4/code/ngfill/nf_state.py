"""
nf_state (NEARGAP_FILL.md §1-§2): near-gap states J = 9, 10, 11, 13, 14, 15 of Omega^EE at Delta_g = 1000 (odd J: the parity-odd 1x1 +-;
even J: the 2x2 helicity block) on their OWN-MASS rows Delta = 1000 + J, by the HIGHSPIN §4 protocol (highspin/hs_even.py, here for
both parities, as EXT_FOLD's ef_j7 did for J = 7): neargap2/ng_bind.run_state UNCHANGED (read-only import; the EE_NEARGAP w2 chain in
EE mode), with
  (a) the own-mass K2 row (nf_rows.py; a coarse 10-node row x = 0.05, 0.1, ..., 0.9 -- J = 15: NEARGAP2's 36-node row) -> "K2 own row"
  (b) the full DIRECT n1 row (nf_points.py, 38 nodes x = 0.025 ... 0.95) injected at every node (direct records M2 = M4 = 0: no Laplace
      truncation on the own row; RQMC-failed nodes, replica rel. err > 0.2, dropped) -> the value and the budget
Budget: NEARGAP2 §1.2, sign-free, x3 on every model term (E_trunc, E_interp), E_low = 1e4 c2 on nu < 50 (gauge factor), E_edge beyond
the last node at max(1e4, 10 max|Gamma| of the last 4 nodes); plus E_mc = 3 x jackknife sigma of lam (8 RQMC replicas).
Variant NU_LOW = 25 (HIGHSPIN's 'NG_NULOW=25'): the direct nodes at nu ~ 25 and 50 stand for the window [25, 50); 1e4 c2 below 25.
Usage: python3 nf_state.py J [mfl,...]   -> cache/nf_state_J<J>.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import nf_core as C
import numpy as np
import ng_bind as NB          # read-only
H, V = C.H, C.V

if __name__ == '__main__':
    J = int(sys.argv[1]); Down = 1000.0 + J; Dg = 1000.0
    mtw = C.twist(J, Dg)
    mfls = [mtw] + ([float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.02, 1.03, 1.05])
    lines, res = [], []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    xd = H.direct_x(J, Down); bad = H.failed_x(J, Down)
    pr(f"# nf_state J={J} own row {Down:g} (m = {C.mphys(Down, J):.3f}, J^2/m = {J*J/C.mphys(Down, J):.3f}, twist edge m/Lam = {mtw:.7f}; "
       f"reference zero x0 = {H.ref_zero_x(J, Down):.3f}): direct nodes {len(xd)} (x {xd.min():.3f} .. {xd.max():.3f}); failed (dropped): {sorted(bad)}")
    _rec0 = H.HN.records

    def _records(kind, J_, D, jk=None, **kw):
        out = _rec0(kind, J_, D, jk, **kw)
        m = H.HN.mphys(D, J_)
        return {k: v for k, v in out.items() if round(k[1]/m, 6) not in bad}
    H.HN.records = _records
    for mfl in mfls:
        pr(f"## m/Lam = {mfl:.7f}")
        NB.NU_LOW = 50.0
        V.INJ.pop(J, None); V.JK[0] = None
        pr("#  K2 rows only (published rows + the own-mass K2 row):")
        o0 = NB.run_state(Dg, J, Down, mfl, pr)
        V.INJ[J] = [Down]
        pr("#  own-mass DIRECT row injected (all replicas):")
        o1 = NB.run_state(Dg, J, Down, mfl, pr)
        lj = []
        for jk in range(8):
            V.JK[0] = jk
            lj.append(NB.run_state(Dg, J, Down, mfl, lambda s: None)['lam_central'])
        V.JK[0] = None
        NB.NU_LOW = 25.0
        pr("#  variant NU_LOW = 25 (direct nodes stand for nu in [25, 50)):")
        o25 = NB.run_state(Dg, J, Down, mfl, pr)
        V.INJ.pop(J, None); NB.NU_LOW = 50.0
        lj = np.array(lj); sj = float(np.sqrt(7/8*((lj - lj.mean())**2).sum()))
        E_mc = 3*sj
        tot = o1['E_tot'] + E_mc; totl = o1['E_tot_lam'] + E_mc; t25 = o25['E_tot'] + E_mc
        o = dict(J=J, mfl=mfl, Down=Down, k2=o0, direct=o1, nu25=o25, jk=lj.tolist(), jk_sigma=sj, E_mc=E_mc, E_tot=tot,
                 lower=o1['lam_central'] - tot, L=o1['lam_central']/tot, E_tot_lam=totl, L_lam=o1['lam_central']/totl,
                 E_tot25=t25, L25=o25['lam_central']/t25, lower25=o25['lam_central'] - t25)
        res.append(o)
        pr(f"  => J={J} m/Lam={mfl:.7f}: published model {o0['lam_model']:+.5f} | K2 own row {o0['lam_K2']:+.5f} (NEARGAP2 budget L {o0['L']:.2f} / "
           f"lambda-level {o0['L_lam']:.2f}) | DIRECT {o1['lam_K2']:+.5f}, central {o1['lam_central']:+.5f} +- {sj:.5f} (jk) | budget sign-free "
           f"{tot:.3e} (trunc {o1['E_trunc']:.2e}, interp {o1['E_interp']:.2e}, low {o1['E_low']:.2e}, edge {o1['E_edge']:.2e} [G {o1['G_edge']:.0f}], "
           f"mc {E_mc:.2e}) -> lower {o['lower']:+.4f}, L = {o['L']:.2f}; lambda-level L = {o['L_lam']:.2f}; NU_LOW=25: low {o25['E_low']:.2e} "
           f"-> L = {o['L25']:.2f}")
        b = o1['base'][list(o1['base'])[0]]
        pr(f"     certified-floor route ({list(o1['base'])[0]}): (Nf/5 - LemmaC)/Ns = {b['rig']:+.4f}")
    H.HN.records = _rec0
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C.NFC, 'nf_state_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.NFC, 'nf_state_J%d.json' % J), 'w'), indent=1, default=float)
