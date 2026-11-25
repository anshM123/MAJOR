"""
b4_state (BATCH4_FOLD.md §2; a copy of b5fold/b5_state.py with Delta_g = 2000): near-gap states at Delta_g = 2000 on their OWN-MASS rows Delta = 2000 + J by the own-row protocol of
ngfill/nf_state.py (= HIGHSPIN §4 / EXT_FOLD's J = 7; neargap2/ng_bind.run_state UNCHANGED, read-only), for both parities:
  even J = 12, 16 (and the seam J = 20; 2x2 helicity block): the batch-4 direct n1 row (x = 0.025 ... 0.95, 38 nodes), the reference-
  zero flanks and the n1 small-nu nodes (nu = 20 ... 47), injected into the own row (J = 12, 16: on top of gapnear's K2 row 2000 + J,
  every one of whose nodes is overridden; J = 20: no K2 row, the shadow b(2020, 20) of b4_shadow.py and 5 K2 EE state records of
  b4_k2.py for the l8 helicity basis);
  odd J (cross-check only; parity-odd 1x1 +-): the direct bulk row + blo + GK15 edge + small-nu nodes injected into the own row;
with:
  * nodes dropped: RQMC-failed (hs_core.failed_x: replica rel. err > 0.2) and nodes closer to the Delta_ref = 2e4 reference zero x0
    (leading ++ wave, the one the batch's flanking 'rz' nodes were placed for) than the flanks, |x - x0| < 0.012 (rule fixed by the job
    design, BATCH5_FOLD §2.2 item 3, before looking at Gamma);
  * the representation (interpolation) term ALSO evaluated on lam_c (NEARGAP_FILL §3 consistency correction).
Budget (NEARGAP2 §1.2 sign-free, x3 on model terms) + E_mc = 3 x jackknife sigma (8 replicas).  The nu < 50 window is then replaced by
the measured window (b4_window.py: smallnu/sn_even.main re-typed with Delta_g = 2000, FULL normalisation).
Env B4_SET: point set (default 'fix' = batch-4 even records mapped to the local structure basis, b4_basis.py; 'all' = unmapped).
Usage: python3 b4_state.py J [mfl,...]   -> cache/state/hs_even_J<J>.json (nf_state / hs_even format), cache/state/b4_state_J<J>.{txt,json}
"""
import os, sys, json, time, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_core as C
import numpy as np
import ng_bind as NB          # read-only
H, V = C.H, C.V
NC = H.NC; W = NC.W
C.SET[0] = os.environ.get('B4_SET', 'fix')
OUT = os.path.join(C.B4C, 'state' if C.SET[0] == 'fix' else 'state_' + C.SET[0]); os.makedirs(OUT, exist_ok=True)
H.HSC = OUT; C.S.HSC = OUT; C.S.SNC = OUT
RZW = 0.012          # strictly closer than the flanking rz nodes (placed at x0 -+ 0.0125, rounded to 4 decimals)


def interp_c(Dg, J, Down, mfl, npn=256):
    """ng_bind.run_state's E_interp probes (linear x-interpolation; every-other-node subsamples of the own row), evaluated on lam_c."""
    Lm = H.E.lam_of(Dg); m = mfl*Lm
    NC.PK = None
    g2, _ = NB.grid_for(J, 2, False)
    own = g2[float(Down)]
    x0, xmax = own['xs'][0], own['xs'][-1]
    d2 = NC.dens(Dg, J, mfl, NB.model_of(J, g2), npn)
    sel = ~((d2['nu'] < max(NB.NU_LOW, x0*m)) | (d2['x'] > xmax))
    lc = NC.lam_of(d2, sel=sel)[0]
    lcl = NC.lam_of(NC.dens(Dg, J, mfl, NB.lin_model(J, g2), npn), sel=sel)[0]
    ls = []
    for keep in (0, 1):
        gs = dict(g2); gs[float(Down)] = NB.sub({float(Down): own}, keep)[float(Down)]
        ls.append(NC.lam_of(NC.dens(Dg, J, mfl, NB.model_of(J, gs), npn), sel=sel)[0])
    return dict(lam_c=float(lc), lam_c_lin=float(lcl), lam_c_sub=[float(v) for v in ls],
                E_interp_c=float(3*max(abs(lc - lcl), *[abs(lc - v) for v in ls])))


if __name__ == '__main__':
    J = int(sys.argv[1]); Down = 2000.0 + J; Dg = 2000.0
    mtw = C.twist(J, Dg)
    mfls = [mtw] + [float(v) for v in (sys.argv[2].split(',') if len(sys.argv) > 2 else ['1.02', '1.03', '1.05']) if float(v) > mtw + 1e-6]
    lines, res = [], []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    m = C.mphys(Down, J)
    xd = H.direct_x(J, Down); failed = H.failed_x(J, Down)
    z = {h: H.ref_zero_x(J, Down, h) for h in (('+-',) if J % 2 else ('++', '+-'))}
    if J % 2:
        z['++'] = z['+-']
    rz = {round(float(x), 6) for x in xd if abs(x - z['++']) < RZW} if z['++'] < 0.98 else set()
    bad = set(failed) | rz
    C.DROPX[(J, Down)] = set(bad)
    pr(f"# b4_state J={J} set {C.SET[0]} own row {Down:g} (m = {m:.3f}, J^2/m = {J*J/m:.3f}, twist edge m/Lam = {mtw:.7f}); reference zeros x0(++) = "
       f"{z['++']:.4f}, x0(+-) = {z['+-']:.4f}; direct nodes {len(xd)} (x {xd.min():.4f} .. {xd.max():.4f}); failed (rel. err > 0.2): "
       f"{sorted(failed)}; dropped by the reference-zero rule |x - x0(++)| < {RZW}: {sorted(rz)}")
    # node table (all direct nodes, before dropping)
    g = H.grid(J, inject=[Down])[Down]
    reps = [H.grid(J, inject=[Down], jk=j)[Down] for j in range(8)]
    keys = ['C2', 'C3H'] if J % 2 else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]
    offk = [] if J % 2 else [('C2', 'off'), ('C3H', 'off')]
    nt = []
    pr("#  node table (Gamma, c2; jackknife sigma over 8 replicas; c_off = m Gamma_off):")
    for i, x in enumerate(g['xs']):
        rec = dict(x=float(x), nu=float(x*m), flag=('FAILED' if round(float(x), 6) in failed else ('RZ' if round(float(x), 6) in rz else '')))
        s = f"   x={x:.4f} nu={x*m:7.1f} {rec['flag']:6s}"
        for k in keys + offk:
            sc = m if isinstance(k, tuple) and k[1] == 'off' else 1.0
            vj = np.array([r[k][int(np.argmin(np.abs(r['xs'] - x)))] for r in reps])*sc
            e = float(np.sqrt(7/8*((vj - vj.mean())**2).sum())); v = float(g[k][i])*sc
            kk = k if isinstance(k, str) else k[0] + k[1]
            rec[kk] = (v, e)
            s += f" | {kk} {v:+10.1f} +- {e:7.1f}"
        nt.append(rec); pr(s)
    _rec0 = H.HN.records

    def _records(kind, J_, D, jk=None, **kw):
        out = _rec0(kind, J_, D, jk, **kw)
        mm = H.HN.mphys(D, J_)
        return {k: v for k, v in out.items() if round(k[1]/mm, 6) not in bad}
    H.HN.records = _records
    g0 = H.grid(J)
    for mfl in mfls:
        pr(f"## m/Lam = {mfl:.7f}")
        NB.NU_LOW = 50.0
        V.INJ.pop(J, None); V.JK[0] = None
        o0 = None
        if g0:
            pr("#  published K2 rows only (no own row):")
            o0 = NB.run_state(Dg, J, Down, mfl, pr)
        V.INJ[J] = [Down]
        pr("#  own-mass DIRECT row injected (all replicas):")
        o1 = NB.run_state(Dg, J, Down, mfl, pr)
        ic = interp_c(Dg, J, Down, mfl)
        lj = []
        for jk in range(8):
            V.JK[0] = jk
            lj.append(NB.run_state(Dg, J, Down, mfl, lambda s: None)['lam_central'])
        V.JK[0] = None
        NB.NU_LOW = 25.0
        o25 = NB.run_state(Dg, J, Down, mfl, lambda s: None)
        V.INJ.pop(J, None); NB.NU_LOW = 50.0
        lj = np.array(lj); sj = float(np.sqrt(7/8*((lj - lj.mean())**2).sum()))
        E_mc = 3*sj
        tot = o1['E_tot'] + E_mc; totl = o1['E_tot_lam'] + E_mc
        totc = tot - o1['E_interp'] + ic['E_interp_c']
        o = dict(J=J, mfl=mfl, Down=Down, k2=o0, direct=o1, nu25=o25, interp_c=ic, jk=lj.tolist(), jk_sigma=sj, E_mc=E_mc, E_tot=tot,
                 lower=o1['lam_central'] - tot, L=o1['lam_central']/tot, E_tot_lam=totl, L_lam=o1['lam_central']/totl,
                 E_tot_c=totc, L_c=o1['lam_central']/totc, held=(len(o1['rows']) == 1), z=z, rz=sorted(rz), failed=sorted(failed))
        res.append(o)
        pr(f"  => J={J} m/Lam={mfl:.7f}: published {o0['lam_model'] if o0 else float('nan'):+.5f} | DIRECT {o1['lam_K2']:+.5f}, central "
           f"{o1['lam_central']:+.5f} +- {sj:.5f} (jk) | sign-free {tot:.3e} (trunc {o1['E_trunc']:.2e}, interp {o1['E_interp']:.2e} "
           f"[on lam_c {ic['E_interp_c']:.2e}], low {o1['E_low']:.2e}, edge {o1['E_edge']:.2e} [G {o1['G_edge']:.0f}, x > {o1['xmax']:.4f}], "
           f"mc {E_mc:.2e}) -> L = {o['L']:.2f} (interp on lam_c: {o['L_c']:.2f}); lambda-level L = {o['L_lam']:.2f}; rows {o1['rows']}"
           + (" [HELD: own row only]" if o['held'] else ""))
        b = o1['base'][list(o1['base'])[0]]
        pr(f"     certified-floor route: (Nf/5 - LemmaC)/Ns = {b['rig']:+.4f}")
    H.HN.records = _rec0
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(OUT, 'b4_state_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    js = dict(J=J, nodes=nt, states=res)
    json.dump(js, open(os.path.join(OUT, 'b4_state_J%d.json' % J), 'w'), indent=1, default=lambda v: float(v) if isinstance(v, (np.floating, np.integer)) else str(v))
    json.dump(res, open(os.path.join(OUT, 'hs_even_J%d.json' % J), 'w'), indent=1, default=lambda v: float(v) if isinstance(v, (np.floating, np.integer)) else str(v))
