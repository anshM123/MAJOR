"""
hs_edgeval (HIGHSPIN.md §3): H_edge on the rows that carry an even-J near-gap state at Delta_g = 1000 (rows not at the state's mass;
1/m interpolation between them, as published) -- exactly NEARGAP2's ng_val protocol (direct n1 edge nodes injected as hc_chain's
HC_INJ does, exempt from the Hampel rule; w2 spline chain, held beyond the last node), with the direct points read from highspin/cache
and every RQMC-failed node (rel. err > 0.2) left out.  Prints Gamma at the direct nodes (K2 where it exists) and the shift of each
state, with its jackknife.
Usage: [HS_IXMIN=0.7] python3 hs_edgeval.py J Delta1[,Delta2] mfl[,mfl...]   -> cache/hs_edgeval_J<J>.txt
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_core as H
import numpy as np
V, NC, HN = H.V, H.NC, H.HN

if __name__ == '__main__':
    J = int(sys.argv[1]); Downs = [float(v) for v in sys.argv[2].split(',')]
    mfls = [float(v) for v in sys.argv[3].split(',')]
    Dg = 1000
    IXMIN = float(os.environ.get('HS_IXMIN', '0'))     # inject only direct nodes with x >= IXMIN (edge only)
    bad = {D: H.failed_x(J, D) for D in Downs}
    _rec0 = HN.records

    def _records(kind, J_, D, jk=None, **kw):
        out = _rec0(kind, J_, D, jk, **kw)
        m = HN.mphys(D, J_)
        return {k: v for k, v in out.items() if round(k[1]/m, 6) not in bad.get(D, set()) and k[1]/m >= IXMIN - 1e-9}
    HN.records = _records
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    V.INJ.pop(J, None)
    gK2 = V.grid(J)
    V.INJ[J] = Downs
    reps = []
    for jk in [None] + list(range(8)):
        V.JK[0] = jk
        reps.append(V.grid(J))
    V.JK[0] = None; V.INJ.pop(J, None)
    keys = ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]
    pr(f"# hs_edgeval J={J} rows {Downs} (model rows {sorted(gK2)}); failed direct nodes dropped: {bad}; reference zero at x = " +
       ", ".join(f"{H.ref_zero_x(J, D):.3f} (row {D:g})" for D in Downs))
    for Down in Downs:
        g0 = gK2[Down]; gd = reps[0][Down]; m = HN.mphys(Down, J)
        P = HN.pts('EE', J, Down)
        xdir = sorted(round(k[1]/m, 6) for k in P if k[0] == Down and round(k[1]/m, 6) not in bad[Down])
        pr(f"#  row {Down:g} (m={m:.2f}, J^2/m={J*J/m:.2f}): Gamma (c2) at the direct nodes (K2 where it exists)")
        for x in xdir:
            i1 = int(np.argmin(np.abs(gd['xs'] - x)))
            if abs(gd['xs'][i1] - x) > 1e-5:
                pr(f"   x={x:.4f}: dropped"); continue
            i0 = int(np.argmin(np.abs(g0['xs'] - x))); k2 = abs(g0['xs'][i0] - x) <= 1e-5
            s = f"   x={x:.4f}:"
            for k in keys:
                vj = np.array([r[Down][k][i1] for r in reps[1:]]); err = np.sqrt(7/8*((vj - vj.mean())**2).sum())
                s += f" {k[0] if isinstance(k, tuple) else k}{'' if not isinstance(k, tuple) else ' ' + k[1]}: " + \
                     (f"{g0[k][i0]:+.0f}/" if k2 else "") + f"{gd[k][i1]:+.0f}+-{err:.0f}"
            pr(s)
    for mfl in mfls:
        lam = []
        for r in reps:
            grids = {'even': {}, 'odd': {}}; grids['odd' if J % 2 else 'even'][J] = r
            em, om = NC.models(grids); f = (om if J % 2 else em)[J][0]
            lam.append(NC.lam_of(NC.dens(Dg, J, mfl, f, 256))[0])
        grids = {'even': {}, 'odd': {}}; grids['odd' if J % 2 else 'even'][J] = gK2
        em, om = NC.models(grids); f = (om if J % 2 else em)[J][0]
        lk = NC.lam_of(NC.dens(Dg, J, mfl, f, 256))[0]
        lj = np.array(lam[1:]); e = np.sqrt(7/8*((lj - lj.mean())**2).sum())
        pr(f"   state m/Lam={mfl:.7f}: lam/N_s K2 rows {lk:+.5f} | with direct edge nodes {lam[0]:+.5f} +- {e:.5f} (jk) | shift {lam[0] - lk:+.5f}")
    pr(f"# Hampel drops (K2 nodes only): {sorted(set((a, b, tuple(c)) for a, b, c in V.DESP_LOG))}")
    pr(f"[{time.time() - t0:.0f}s]")
    HN.records = _rec0
    open(os.path.join(H.HSC, 'hs_edgeval_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
