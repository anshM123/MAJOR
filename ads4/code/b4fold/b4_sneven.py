"""
b4_sneven (BATCH4_FOLD.md §2): a copy of smallnu/sn_even.py with Delta_g = 2000 (Lam = lam_of(DG), sn_window.DG = DG, the
twist edge from b4_core), otherwise unchanged; outputs to sn_core.SNC (b4_core redirects it).

sn_even (SMALLNU_WINDOW.md): even-J near-gap states with the nu < 50 window from MEASURED small-nu direct values.
  * own-row states J = 12, 16 (rows 1012, 1016; HIGHSPIN §4 budgets read from highspin/cache/hs_even_J<J>.json, whose central value
    already excludes the window, ng_bind's lam_central): E_low replaced, everything else unchanged.
  * J = 20 (rows 800, 1600), 32 (row 1040): no budget exists; node tables and window costs only.
Node values: sn_core.small_nodes on each row (baseline normalisation, jackknife); diagonal ++, +- and c_off = m Gamma_off.
Window treatments (sn_window; NEARGAP2's E_low structure: |weights| x gauge factor, normalised 2x2, Frobenius norm):
  conv    : |Gamma| <= 1e4 for every entry (the old convention; c_off included with the same G, as ng_bind)
  env-k   : per-ENTRY piecewise measured envelope (++, +-: max over C2/C3H of |value| + 3 sigma; off: |c_off| + 3 sigma), k = 1, 3, 10,
            1e4 below each row's lowest reliable direct node; the state's envelope at x = max over its rows (rows used at fixed x)
  signed  : the measured window matrix dF_win (direct values, natural spline in x per row and entry, held beyond the last node, rows
            1/m-interpolated as the chain does) -> its spectral norm |dF_win|_2 (largest possible first-order effect on lambda, either
            sign) + 3 x MC (sum |w_i| sigma_i, fully correlated) + 3 x |spline - linear|, Frobenius; 1e4 below nu_min.
Usage: python3 sn_even.py J row[,row...] [nu_min]   -> cache/sn_even_J<J>.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_core as C
import sn_core as S
import sn_window as SW
DG = C.DG
SW.DG = DG
import numpy as np
from scipy.interpolate import CubicSpline
H = S.H
G_LO = 1e4
RELMAX = 0.2
TW = {12: 1.0119407, 16: 1.0158933, 20: 1.0198198, 32: 1.0315315}
KEYS = [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-', 'off')]
ENT = ('++', '+-', 'off')


def m_row(D, J):
    return H.HN.mphys(D, J)


def main(J, rows, numin_req=None):
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    Lm = H.E.lam_of(DG); mtw = C.twist(J)*Lm
    res = dict(J=J, rows=rows, nodes={}, states={})
    pr(f"# sn_even J={J} rows {rows}: small-nu nodes (c2, baseline normalisation; c_off = m Gamma_off)")
    good = {}
    for D in rows:
        nodes, chk, da, m = S.small_nodes(J, D, numax=80.0*m_row(D, J)/mtw)
        res['nodes'][D] = [{(k if isinstance(k, str) else k[0] + k[1]): v for k, v in r.items()} for r in nodes]
        pr(f"#  row {D:g} (m={m:.2f}): transform check {chk:.2e}; normalisation change r-1 {da}")
        good[D] = []
        for r in nodes:
            rel = S.rel_err(J, D, r['nu']) if r['src'] == 'sn' else 0.0
            ok = r['src'] != 'K2' and rel <= RELMAX
            pr(f"   nu={r['nu']:7.2f} x={r['x']:.5f} (state nu {r['x']*mtw:6.1f}) [{r['src']}] " +
               "  ".join(f"{k[0]}{k[1]}: {r[k][0]:+8.1f}+-{r[k][1]:5.1f}" for k in KEYS if k in r) +
               (f"  (rel {rel:.1e}{'' if ok else ' NOT USED'})" if r['src'] == 'sn' else ''))
            if ok and (numin_req is None or r['x']*mtw >= float(numin_req) - 1e-9):
                good[D].append(r)
    numin = min(r['x']*mtw for D in rows for r in good[D])
    bnd = {D: {e: np.array([max(abs(r[(kd, e)][0]) + 3*r[(kd, e)][1] for kd in ('C2', 'C3H')) for r in good[D]]) for e in ENT}
           for D in rows}
    gx = {D: np.array([r['x'] for r in good[D]]) for D in rows}
    win = [(D, i) for D in rows for i, r in enumerate(good[D]) if r['x']*mtw <= 52.0]
    Gd = max(max(bnd[D]['++'][i], bnd[D]['+-'][i]) for D, i in win)
    Go = max(bnd[D]['off'][i] for D, i in win)
    Gc = max(abs(good[D][i][k][0]) for D, i in win for k in KEYS if k[1] != 'off')
    res.update(numin=numin, G_meas=Gd, G_meas_off=Go, G_cent=Gc)
    pr(f"#  nu_min (state) = {numin:.2f}; window nodes {len(win)}; max (|Gamma| + 3 sigma) diag {Gd:.0f} (max |Gamma| {Gc:.0f}), "
       f"|c_off| + 3 sigma {Go:.1f} c2")
    for D in rows:
        pr(f"#  node bounds row {D:g} (state nu: ++ / +- / off): " +
           ", ".join(f"{x*mtw:.0f}:{a:.0f}/{b:.0f}/{c:.0f}" for x, a, b, c in zip(gx[D], bnd[D]['++'], bnd[D]['+-'], bnd[D]['off'])))

    def gfun(k):
        fs = {D: {e: SW.envelope(gx[D], bnd[D][e], k, G_LO) for e in ENT} for D in rows}
        return lambda e, x: np.max([fs[D][e](x) for D in rows], axis=0)
    xmin = max(gx[D].min() for D in rows)
    # signed: splines per row and key, 1/m interpolation between rows
    ms = np.array([m_row(D, J) for D in rows])
    spl = {D: {k: CubicSpline(gx[D], [r[k][0] for r in good[D]], bc_type='natural') for k in KEYS} for D in rows}
    sig = {D: {k: np.array([r[k][1] for r in good[D]]) for k in KEYS} for D in rows}

    def tw(mfl):
        m = mfl*Lm
        if len(rows) == 1:
            return {rows[0]: 1.0}
        o = np.argsort(ms)
        if m <= ms[o[0]]:
            return {rows[o[0]]: 1.0}
        if m >= ms[o[-1]]:
            return {rows[o[-1]]: 1.0}
        i1 = int(np.searchsorted(ms[o], m)); a, b = rows[o[i1 - 1]], rows[o[i1]]
        t = (1/m_row(a, J) - 1/m)/(1/m_row(a, J) - 1/m_row(b, J))
        return {a: 1 - t, b: t}

    fj = os.path.join(H.HSC, 'hs_even_J%d.json' % J)
    old = {o['mfl']: o for o in json.load(open(fj))} if os.path.exists(fj) else {}
    mfls = list(old) or [C.twist(J), 1.03, 1.05]
    for mfl in mfls:
        mst = mfl*Lm; T = tw(mfl)
        E_conv = SW.cost(J, mfl, G_LO); E_lo = SW.cost(J, mfl, G_LO, 0.0, xmin*mst)
        o = dict(mfl=mfl, E_conv=E_conv, E_lo=E_lo, W=SW.weight(J, mfl), rows_t={str(k): v for k, v in T.items()})
        for k in (1, 3, 10):
            o[f'E_env{k}'] = SW.cost_entries(J, mfl, gfun(k))

        def gam(key, x, mode='spl'):
            out = 0.0
            for D, t in T.items():
                xc = np.clip(x, gx[D].min(), gx[D].max())
                v = spl[D][key](xc) if mode == 'spl' else np.interp(xc, gx[D], spl[D][key](gx[D]))
                out = out + t*v
            return out
        dF, base, Ns = SW.signed_dF(J, mfl, gam, xmin*mst, SW.NU_W)
        dFl, _, _ = SW.signed_dF(J, mfl, lambda k, x: gam(k, x, 'lin'), xmin*mst, SW.NU_W)
        nn = 1/np.sqrt(Ns)
        dFn = nn[:, None]*dF*nn[None, :]
        Em = 3*np.abs(nn[:, None]*(dF - dFl)*nn[None, :])
        for D, t in T.items():                       # MC: each node of each row, fully correlated
            for i in range(len(gx[D])):
                e_i = np.zeros(len(gx[D])); e_i[i] = 1.0
                for key in KEYS:
                    if sig[D][key][i] == 0:
                        continue

                    def g1(k2, x, key=key, D=D, e_i=e_i, t=t):
                        return t*np.interp(np.clip(x, gx[D].min(), gx[D].max()), gx[D], e_i) if k2 == key else 0.0*x
                    d1, _, _ = SW.signed_dF(J, mfl, g1, xmin*mst, SW.NU_W)
                    Em += 3*sig[D][key][i]*np.abs(nn[:, None]*d1*nn[None, :])
        o['dFwin_norm'] = dFn.tolist(); o['dFwin_spec'] = float(np.max(np.abs(np.linalg.eigvalsh(dFn))))
        o['E_signed'] = float(np.sqrt((Em*Em).sum()))
        o['lam_hyb'] = float(H.E.lam2(base, Ns, 0*dF)); o['dlam_win_hyb'] = float(H.E.lam2(base, Ns, dF) - o['lam_hyb'])
        o['E_meas'] = o['dFwin_spec'] + o['E_signed'] + E_lo
        s = (f"  m/Lam={mfl:.7f} (rows {', '.join(f'{D:g}:{t:.3f}' for D, t in T.items())}): window conv 1e4 {E_conv:.3f} | below nu_min "
             f"{E_lo:.4f} | env x1 {o['E_env1']:.3f} x3 {o['E_env3']:.3f} x10 {o['E_env10']:.3f} | measured window: |dF_win|_2 "
             f"{o['dFwin_spec']:.4f} (first-order on the hybrid {o['dlam_win_hyb']:+.4f}) + error {o['E_signed']:.4f} + below nu_min -> {o['E_meas']:.3f}")
        if mfl in old:
            ob = old[mfl]; di = ob['direct']
            lamc = di['lam_central']
            rest = ob['E_tot'] - di['E_low']; rest_l = ob['E_tot_lam'] - di['E_low']
            o.update(lamc=lamc, E_low_hs=di['E_low'], rest=rest, rest_lam=rest_l, L_old=ob['L'], L_old_lam=ob['L_lam'])
            for k in (1, 3, 10):
                o[f'L_env{k}'] = lamc/(rest + o[f'E_env{k}']); o[f'Llam_env{k}'] = lamc/(rest_l + o[f'E_env{k}'])
            o['L_meas'] = lamc/(rest + o['E_meas']); o['Llam_meas'] = lamc/(rest_l + o['E_meas'])
            o['G_tol'] = (lamc - rest)/o['W']
            s += (f"\n        lam_c {lamc:+.4f}; HIGHSPIN E_low {di['E_low']:.3f}, other terms {rest:.3f} (lambda-level {rest_l:.3f}); "
                  f"L sign-free old {ob['L']:.2f} -> env3 {o['L_env3']:.2f} / env10 {o['L_env10']:.2f} / measured {o['L_meas']:.2f}; "
                  f"lambda-level old {ob['L_lam']:.2f} -> {o['Llam_env3']:.2f} / {o['Llam_env10']:.2f} / {o['Llam_meas']:.2f}; "
                  f"tolerated uniform |Gamma| on nu<50 {o['G_tol']:.3g}")
        res['states'][f'{mfl:.7f}'] = o
        pr(s)
    open(os.path.join(S.SNC, 'sn_even_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(S.SNC, 'sn_even_J%d.json' % J), 'w'), indent=1,
              default=lambda v: float(v) if isinstance(v, (np.floating, np.integer)) else str(v))


if __name__ == '__main__':
    t0 = time.time()
    main(int(sys.argv[1]), [float(v) for v in sys.argv[2].split(',')], sys.argv[3] if len(sys.argv) > 3 else None)
    print(f"[{time.time() - t0:.0f}s]")
