"""
hs_even (HIGHSPIN.md §4): near-gap even-J states of Omega^EE at Delta_g = 1000 (J = 12, 16, ...) on their own-mass K2 rows
(neargap2/ng_rows, Delta = 1000 + J), with DIRECT n1 nodes (hs_points.py) injected at every node where they exist (hc_nodes records,
M2 = M4 = 0, so the Laplace truncation is absent there) -- through neargap2/ng_bind.run_state, unchanged (the injection is done by the
loader wrappers of neargap2/ng_val, installed by hs_core).  Same budget as NEARGAP2 §1.2 (x3 on the model terms), plus
  E_mc = 3 x jackknife sigma of lambda over the 8 RQMC replicas (every replica of every injected node left out together).
NEARGAP2's K2-only values are reproduced first (no injection).
Usage: [NG_NULOW=25] python3 hs_even.py J Delta_own mfl[,mfl...] [tag]   -> cache/hs_even_<tag>.txt/.json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_core as H
import numpy as np
import ng_bind as NB
V = H.V

if __name__ == '__main__':
    J, Down = int(sys.argv[1]), float(sys.argv[2])
    mfls = [float(v) for v in sys.argv[3].split(',')]
    tag = sys.argv[4] if len(sys.argv) > 4 else 'J%d' % J
    Dg = 1000.0
    lines, res = [], []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    xd = H.direct_x(J, Down); bad = H.failed_x(J, Down)
    pr(f"# hs_even J={J} own row {Down:g}: direct nodes {len(xd)} (x {xd.min() if len(xd) else 0:.3f} .. {xd.max() if len(xd) else 0:.3f}); "
       f"failed (dropped from the injection: rel. err > 0.2): {sorted(bad)}; NU_LOW = {NB.NU_LOW}")
    # failed nodes: hide them from hc_nodes.records by a numin/numax-free filter -> patch HN.records to skip them
    _rec0 = H.HN.records

    def _records(kind, J_, D, jk=None, **kw):
        out = _rec0(kind, J_, D, jk, **kw)
        m = H.HN.mphys(D, J_)
        return {k: v for k, v in out.items() if round(k[1]/m, 6) not in bad}
    H.HN.records = _records
    for mfl in mfls:
        pr(f"## m/Lam = {mfl}")
        V.INJ.pop(J, None); V.JK[0] = None
        pr("#  K2 only (NEARGAP2 reproduction):")
        o0 = NB.run_state(Dg, J, Down, mfl, pr)
        V.INJ[J] = [Down]
        pr("#  direct nodes injected (all replicas):")
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
        pr(f"  => J={J} m/Lam={mfl:.7f}: K2 row {o0['lam_K2']:+.5f} (NEARGAP2 budget L {o0['L']:.1f} / lambda-level {o0['L_lam']:.1f}) | "
           f"direct-injected {o1['lam_K2']:+.5f}, central {o1['lam_central']:+.5f} +- {sj:.5f} (jk) | budget sign-free {tot:.3e} "
           f"(trunc {o1['E_trunc']:.2e}, interp {o1['E_interp']:.2e}, low {o1['E_low']:.2e}, edge {o1['E_edge']:.2e}, mc {E_mc:.2e}) "
           f"-> lower {o['lower']:+.4f}, L = {o['L']:.1f}; lambda-level {totl:.3e} -> L = {o['L_lam']:.1f}")
    H.HN.records = _rec0
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(H.HSC, 'hs_even_%s.txt' % tag), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(H.HSC, 'hs_even_%s.json' % tag), 'w'), indent=1, default=float)
