"""
sn_even (SMALLNU_WINDOW.md): even-J near-gap states with the nu < 50 window from MEASURED small-nu direct values.
  * own-row states J = 12, 16 (rows 1012, 1016; HIGHSPIN §4 budgets read from highspin/cache/hs_even_J<J>.json, whose central value
    already excludes the window, ng_bind's lam_central): E_low replaced, everything else unchanged.
  * J = 20, 32 (rows 800/1600, 1040/...; no budget exists): node table and window cost only.
Node values: sn_core.small_nodes on each row (baseline normalisation, jackknife); diagonal ++, +- and c_off = m Gamma_off.
Envelope: piecewise measured (sn_window.envelope) per row, node bound |value| + 3 sigma (max over kinds/entries), k = 1, 3, 10;
  1e4 below each row's lowest reliable direct node; the state's envelope at x = max over its rows (the chain interpolates rows at
  fixed x).  Window cost: sn_window.cost (NEARGAP2's E_low convention incl. the gauge factor).
Usage: python3 sn_even.py J row[,row...] [nu_min]   -> cache/sn_even_J<J>.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sn_core as S
import sn_window as SW
import numpy as np
H = S.H
G_LO = 1e4
RELMAX = 0.2
TW = {12: 1.0119407, 16: 1.0158933, 20: 1.0198198, 32: 1.0315315}


def main(J, rows, numin_req=None):
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    res = dict(J=J, rows=rows, nodes={}, states={})
    inwin = []
    pr(f"# sn_even J={J} rows {rows}: small-nu nodes (c2, baseline normalisation; c_off = m Gamma_off)")
    Lm = H.E.lam_of(1000.0)
    keys = [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-', 'off')]
    allgood = []
    for D in rows:
        nodes, chk, da, m = S.small_nodes(J, D, numax=80.0*m_row(D, J)/(TW.get(J, 1.02)*Lm))
        res['nodes'][D] = [{(k if isinstance(k, str) else k[0] + k[1]): v for k, v in r.items()} for r in nodes]
        pr(f"#  row {D:g} (m={m:.2f}): transform check {chk:.2e}; normalisation change r-1 {da}")
        for r in nodes:
            rel = S.rel_err(J, D, r['nu']) if r['src'] == 'sn' else 0.0
            ok = r['src'] != 'K2' and rel <= RELMAX
            nus = r['x']*TW.get(J, 1.02)*Lm               # the state's nu at this x (the chain interpolates rows at fixed x)
            pr(f"   nu={r['nu']:7.2f} x={r['x']:.5f} (state nu {nus:6.1f}) [{r['src']}] " +
               "  ".join(f"{k[0]}{k[1]}: {r[k][0]:+8.1f}+-{r[k][1]:5.1f}" for k in keys if k in r) +
               (f"  (rel {rel:.1e}{'' if ok else ' NOT USED'})" if r['src'] == 'sn' else ''))
            if ok:
                allgood.append((nus, r, D))
    numin = float(numin_req) if numin_req else min(a[0] for a in allgood)
    allgood = [(v, r, D) for v, r, D in allgood if v >= numin - 1e-9]
    inwin = [r for v, r, D in allgood if v <= 52.0]
    Gm = max(abs(r[k][0]) + 3*r[k][1] for r in inwin for k in keys if k in r)
    Gc = max(abs(r[k][0]) for r in inwin for k in keys if k in r)
    res.update(numin=numin, G_meas=Gm, G_cent=Gc)
    pr(f"#  nu_min (state) = {numin:.2f}; window nodes {len(inwin)}; max (|value| + 3 sigma) over the window = {Gm:.0f} (max |value| {Gc:.0f}) c2")
    byrow = {}
    for v, r, D in allgood:
        byrow.setdefault(D, []).append((r['x'], max(abs(r[k][0]) + 3*r[k][1] for k in keys if k in r)))
    pr("#  node bounds (state nu: |value| + 3 sigma): " + "; ".join(f"row {D:g}: " + ", ".join(f"{x*TW.get(J, 1.02)*Lm:.0f}:{g:.0f}" for x, g in sorted(v)) for D, v in byrow.items()))

    def env(k):
        fs = [SW.envelope([x for x, g in v], [g for x, g in v], k, G_LO) for v in byrow.values()]
        return lambda x: np.max([f(x) for f in fs], axis=0)
    xmin = max(min(x for x, g in v) for v in byrow.values())
    fj = os.path.join(H.HSC, 'hs_even_J%d.json' % J)
    old = {o['mfl']: o for o in json.load(open(fj))} if os.path.exists(fj) else {}
    mfls = [o for o in old] or [TW[J]]
    for mfl in mfls:
        mst = mfl*Lm
        E_conv = SW.cost(J, mfl, G_LO); E_lo = SW.cost(J, mfl, G_LO, 0.0, xmin*mst)
        o = dict(mfl=mfl, E_conv=E_conv, E_lo=E_lo, W=SW.weight(J, mfl))
        for k in (1, 3, 10):
            o[f'E_env{k}'] = SW.cost_x(J, mfl, env(k))
        s = f"  m/Lam={mfl:.7f}: window conv 1e4 {E_conv:.3f} | below nu_min {E_lo:.4f} | env x3 {o['E_env3']:.3f} | env x10 {o['E_env10']:.3f}"
        if mfl in old:
            ob = old[mfl]; di = ob['direct']
            lamc = di['lam_central']
            rest = ob['E_tot'] - di['E_low']; rest_l = ob['E_tot_lam'] - di['E_low']
            o.update(lamc=lamc, E_low_hs=di['E_low'], rest=rest, rest_lam=rest_l, L_old=ob['L'], L_old_lam=ob['L_lam'])
            for k in (3, 10):
                o[f'L_env{k}'] = lamc/(rest + o[f'E_env{k}']); o[f'Llam_env{k}'] = lamc/(rest_l + o[f'E_env{k}'])
            o['G_tol'] = (lamc - rest)/o['W']
            s += (f" || lam_c {lamc:+.4f}; HIGHSPIN E_low {di['E_low']:.3f}, other terms {rest:.3f} (lambda-level {rest_l:.3f}); "
                  f"L sign-free old {ob['L']:.2f} -> env3 {o['L_env3']:.2f} / env10 {o['L_env10']:.2f}; lambda-level old {ob['L_lam']:.2f} -> "
                  f"{o['Llam_env3']:.2f} / {o['Llam_env10']:.2f}; tolerated uniform |Gamma| on nu<50 {o['G_tol']:.3g}")
        res['states'][f'{mfl:.7f}'] = o
        pr(s)
    open(os.path.join(S.SNC, 'sn_even_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(S.SNC, 'sn_even_J%d.json' % J), 'w'), indent=1,
              default=lambda v: float(v) if isinstance(v, (np.floating, np.integer)) else str(v))


def m_row(D, J):
    return H.HN.mphys(D, J)


if __name__ == '__main__':
    t0 = time.time()
    main(int(sys.argv[1]), [float(v) for v in sys.argv[2].split(',')], sys.argv[3] if len(sys.argv) > 3 else None)
    print(f"[{time.time() - t0:.0f}s]")
