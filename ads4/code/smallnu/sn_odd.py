"""
sn_odd (SMALLNU_WINDOW.md): odd-J near-gap states (J = 21, 23, 27, 29, 31; rows Delta = 1000 + J) with the nu < 50 window from
MEASURED small-nu direct values instead of the |Gamma| <= 1e4 c2 convention.

Inputs (read-only): HIGHSPIN's hs_odd results (highspin/cache/hs_odd_J<J>[_final].json: main value, lambda-level budget terms
E_mc, E_quad, E_miss, E_trunc, E_repr), recomputed here only to (a) check the main value is reproduced and (b) get X_low, the bulk
model's own contribution on nu < 50 (hs_odd keeps the shape-fit continuation there in its central value; here it is REMOVED and
the window is covered by the envelope, as in ng_bind's lam_central).
Small-nu values: sn_core.small_nodes (direct sn points + HIGHSPIN direct nodes, baseline normalisation, jackknife).
Window treatments (all with NEARGAP2's gauge factor, sn_window.cost):
  conv   : G = 1e4 on nu < 50                                           (the old convention, for comparison)
  env-k  : piecewise MEASURED envelope (sn_window.envelope): node bound g_i = max over kinds of |Gamma_i| + 3 sigma_i; on [x_i, x_i+1]
           k x max(g_i, g_i+1) (k = 1 bare, 3 main, 10 conservative); below the lowest reliable direct node (nu_min)
           G_lo = 1e4 (the old convention, kept ONLY there)
  signed : the direct values on [nu_min, 50) integrated (natural spline in nu through the direct nodes up to nu = 80) in the central
           value; error 3 x (sum |w_i| sigma_i, fully correlated) + 3 x |spline - linear|; G_lo below nu_min.
Usage: python3 sn_odd.py J [nu_min]   -> cache/sn_odd_J<J>.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sn_core as S
import sn_window as SW
import numpy as np
from numpy.polynomial.legendre import leggauss
H = S.H
import hs_odd as HO          # read-only import (functions only; its __main__ is not run)

DG = 1000.0
G_LO = 1e4
RELMAX = 0.2                 # direct nodes with replica error > 20 % (state or reference) are not used


def main(J, numin_req=None):
    D = 1000.0 + J
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    fj = os.path.join(S.SNC, 'hs_odd_J%d.json' % J)          # this note's own hs_odd run (direct bulk row in smallnu/cache), if any
    if os.path.exists(fj):
        S.FULL[0] = True
    else:
        fj = os.path.join(H.HSC, 'hs_odd_J%d_final.json' % J)
        if not os.path.exists(fj):
            fj = os.path.join(H.HSC, 'hs_odd_J%d.json' % J)
    old = json.load(open(fj))
    nodes, chk, da, m = S.small_nodes(J, D, numax=80.0)
    pr(f"# sn_odd J={J} row {D:g} (m={m:.3f}); HIGHSPIN results from {os.path.basename(fj)}; normalisation transform check on old "
       f"nodes {chk:.2e} c2; a_new/a_base - 1 = {da:+.2e}")
    pr("#  small-nu nodes (Gamma - kappa, c2, baseline normalisation; src: sn = this note, hs = HIGHSPIN direct, K2):")
    good = []
    for r in nodes:
        rel = S.rel_err(J, D, r['nu']) if r['src'] == 'sn' else 0.0
        ok = r['src'] != 'K2' and rel <= RELMAX
        pr(f"   nu={r['nu']:7.2f} x={r['x']:.5f} [{r['src']}] " + "  ".join(f"{k}: {r[k][0]:+9.1f} +- {r[k][1]:6.1f}" for k in ('C2', 'C3H'))
           + (f"  (rel err {rel:.1e}{'' if ok else ' -> NOT USED'})" if r['src'] == 'sn' else ''))
        if ok:
            good.append(r)
    gx = np.array([r['x'] for r in good])
    gb = np.array([max(abs(r[k][0]) + 3*r[k][1] for k in ('C2', 'C3H')) for r in good])
    numin = float(numin_req) if numin_req else float(min(r['nu'] for r in good))
    good = [r for r in good if r['nu'] >= numin - 1e-9]
    gx = np.array([r['x'] for r in good]); gb = np.array([max(abs(r[k][0]) + 3*r[k][1] for k in ('C2', 'C3H')) for r in good])
    inwin = [r for r in good if r['nu'] <= 52.0]
    Gm = max(abs(r[k][0]) + 3*r[k][1] for r in inwin for k in ('C2', 'C3H'))
    Gc = max(abs(r[k][0]) for r in inwin for k in ('C2', 'C3H'))
    pr(f"#  nu_min = {numin:.2f}; window nodes {len(inwin)}; max over the window of (|Gamma| + 3 sigma) = {Gm:.0f} c2 (max |Gamma| {Gc:.0f}); "
       f"node bounds: " + ", ".join(f"{r['nu']:.0f}:{g:.0f}" for r, g in zip(good, gb)))
    # bulk model (baseline points) for X_low and the reproduction check
    S.BASE[0] = not S.FULL[0]
    try:
        rK, rI, isdir = HO.rows(J, D)
        bad = H.refbad(J, D, HO.RTOL)
        nbk = int(((rK['xs'] >= 0.05 - 1e-9) & (rK['xs'] < HO.XCUT - 1e-9)).sum())
        nbd = int(((rI['xs'] >= 0.05 - 1e-9) & (rI['xs'] < HO.XCUT - 1e-9) & isdir).sum())
        bulk_main = 'dir' if nbd >= 0.9*nbk else 'K2'
        M = HO.Model(J, rK, rI, isdir, bad, bulk=bulk_main)
    finally:
        S.BASE[0] = False
    Lm = H.E.lam_of(DG); rowm = H.HN.mphys(D, J); HO.XROW[0] = Lm/rowm
    res = dict(J=J, D=D, numin=numin, G_meas=Gm, G_cent=Gc, nodes=nodes, states={}, chk=chk, da=da)
    # spline of the direct values in nu (signed variant)
    from scipy.interpolate import CubicSpline
    sx = np.array([r['nu'] for r in good]); o = np.argsort(sx); sx = sx[o]
    sp = {k: CubicSpline(sx, np.array([good[i][k][0] for i in o]), bc_type='natural') for k in ('C2', 'C3H')}
    se = {k: np.array([good[i][k][1] for i in o]) for k in ('C2', 'C3H')}
    for key, s in old['states'].items():
        mfl = float(key)
        edge = s['main_key'].split('/')[1]
        lam_re = HO.value(J, mfl, M, edge)[0]
        hs, base, Ns, st = H.base_of(DG, J, mfl)
        mst = mfl*Lm; x50 = SW.NU_W/mst
        t, w = leggauss(400)
        xx = 0.5*x50*(t + 1); wx = 0.5*x50*w
        w2, w3 = H.state_weights(DG, J, mfl, '+-', xx)
        X_low = float(((wx*w2*M.k['C2']['bulk'](xx)).sum() + (wx*w3*M.k['C3H']['bulk'](xx)).sum())/Ns[0])
        lamc = s['main'] - X_low
        Eb = s['E_tot']
        nuq = xx*m                       # ROW nu at the state's x (the chain uses the row at fixed x)
        xmin = gx.min(); mid = xx >= xmin
        # signed: spline of the direct values on [numin, 50)
        Dmid = float(((wx*w2*sp['C2'](nuq))[mid].sum() + (wx*w3*sp['C3H'](nuq))[mid].sum())/Ns[0])
        lin = float(((wx*w2*np.interp(nuq, sx, sp['C2'](sx)))[mid].sum() + (wx*w3*np.interp(nuq, sx, sp['C3H'](sx)))[mid].sum())/Ns[0])
        # MC: effective node weights of linear interpolation, fully correlated (sum |w_i| sigma_i)
        emc = 0.0
        for k, wk in (('C2', w2), ('C3H', w3)):
            for i in range(len(sx)):
                e_i = np.zeros(len(sx)); e_i[i] = 1.0
                emc += abs(float((wx*wk*np.interp(nuq, sx, e_i))[mid].sum()))/Ns[0]*se[k][i]
        E_lo = SW.cost(J, mfl, G_LO, 0.0, xmin*mst)
        E_conv = SW.cost(J, mfl, G_LO)
        E_conv_ng = SW.cost(J, mfl, G_LO, gauge=False)
        W_mid = SW.weight(J, mfl, xmin*mst, SW.NU_W)
        W_all = SW.weight(J, mfl)
        o = dict(mfl=mfl, main=s['main'], lam_re=lam_re, X_low=X_low, lamc=lamc, E_base=Eb, E_trunc=s['E_trunc'],
                 E_conv=E_conv, E_conv_nogauge=E_conv_ng, W_all=W_all, W_mid=W_mid, E_lo=E_lo, D_mid=Dmid, E_signed=3*emc + 3*abs(Dmid - lin))
        for k in (1, 3, 10):
            Ek = SW.cost_x(J, mfl, SW.envelope(gx, gb, k, G_LO))
            o[f'E_env{k}'] = Ek
            o[f'L_env{k}'] = lamc/(Eb + Ek)
            o[f'Lcal_env{k}'] = lamc/(Eb - s['E_trunc'] + 0.05*s['E_trunc'] + Ek)
        o['L_conv'] = lamc/(Eb + E_conv)
        o['lam_signed'] = lamc + Dmid
        o['L_signed'] = o['lam_signed']/(Eb + E_lo + o['E_signed'])
        o['Lcal_signed'] = o['lam_signed']/(Eb - s['E_trunc'] + 0.05*s['E_trunc'] + E_lo + o['E_signed'])
        o['G_tol'] = (lamc - Eb)/W_all
        o['G_tol_mid'] = (lamc - Eb - E_lo)/W_mid
        o['G_lo_tol3'] = (lamc - Eb - (o['E_env3'] - E_lo))/max(SW.weight(J, mfl, 0.0, xmin*mst), 1e-30)
        res['states'][key] = o
        pr(f"  m/Lam={mfl:.7f}: main {s['main']:+.4f} (recomputed {lam_re:+.4f}) | X_low(model, nu<50) {X_low:+.4f} -> lam_c {lamc:+.4f} | "
           f"budget w/o window {Eb:.3f} (of which trunc {s['E_trunc']:.3f})")
        pr(f"        window: conv 1e4 {E_conv:.3f} (no gauge {E_conv_ng:.3f}) -> L {o['L_conv']:.2f} | below nu_min (1e4) {E_lo:.4f} | "
           f"env x1 {o['E_env1']:.3f} | env x3 {o['E_env3']:.3f} -> L {o['L_env3']:.2f} | env x10 {o['E_env10']:.3f} -> L {o['L_env10']:.2f} | signed: D_mid "
           f"{Dmid:+.4f} +- {o['E_signed']:.4f} -> lam {o['lam_signed']:+.4f}, L {o['L_signed']:.2f}")
        pr(f"        W(nu<50) {W_all:.2e}, W[nu_min,50) {W_mid:.2e} per c2; tolerated uniform |Gamma| on nu<50 {o['G_tol']:.3g}, on "
           f"[nu_min,50) with 1e4 below {o['G_tol_mid']:.3g} c2; tolerated below nu_min with env x3 above: {o['G_lo_tol3']:.3g} | calibrated-trunc (5 % of step) L: env3 {o['Lcal_env3']:.2f}, "
           f"env10 {o['Lcal_env10']:.2f}, signed {o['Lcal_signed']:.2f}")
    open(os.path.join(S.SNC, 'sn_odd_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(S.SNC, 'sn_odd_J%d.json' % J), 'w'), indent=1, default=lambda v: v if not isinstance(v, (np.floating, np.integer)) else float(v))


if __name__ == '__main__':
    t0 = time.time()
    main(int(sys.argv[1]), sys.argv[2] if len(sys.argv) > 2 else None)
    print(f"[{time.time() - t0:.0f}s]")
