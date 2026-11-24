"""
b5_window (BATCH5_FOLD.md §2): the measured nu < 50 window for the own-row (b5_state.py) states, replacing ng_bind's E_low = 1e4 c2:
  even J: smallnu/sn_even.main UNCHANGED (read-only import; per-entry envelope x1/x3/x10 with 1e4 below the lowest reliable node, and
          the 'measured' window = spectral norm of the signed window matrix + its error), on b5_state's hs_even_J<J>.json, FULL
          normalisation (the row exists only with its direct points: there is no baseline grid);
  odd J (own-row route, J = 17, 19, 25): smallnu/sn_odd's envelope and signed-window formulas (copied below, with sn_core.small_nodes /
          sn_window unchanged) applied to the ng_bind central value lam_c (which already excludes nu < 50).
Env B5_SET as b5_state.py.
Usage: python3 b5_window.py J   -> cache/state/sn_even_J<J>.{txt,json} (even) | cache/state/b5_window_J<J>.{txt,json} (odd)
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_core as C
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
H, S = C.H, C.S
C.SET[0] = os.environ.get('B5_SET', 'fix')
OUT = os.path.join(C.B5C, 'state' if C.SET[0] == 'fix' else 'state_' + C.SET[0])
H.HSC = OUT; S.HSC = OUT; S.SNC = OUT
S.FULL[0] = True
import sn_window as SW
G_LO = 1e4
RELMAX = 0.2


def odd(J):
    D = 1000.0 + J; Lm = H.E.lam_of(1000.0)
    st = json.load(open(os.path.join(OUT, 'hs_even_J%d.json' % J)))
    lines, res = [], dict(J=J, states={})

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    nodes, chk, da, m = S.small_nodes(J, D, numax=80.0)
    good = []
    pr(f"# b5_window J={J} (odd, own-row route) row {D:g}: small-nu nodes (Gamma, c2, FULL normalisation)")
    for r in nodes:
        rel = S.rel_err(J, D, r['nu']) if r['src'] == 'sn' else 0.0
        ok = r['src'] != 'K2' and rel <= RELMAX
        pr(f"   nu={r['nu']:7.2f} x={r['x']:.5f} [{r['src']}] " + "  ".join(f"{k}: {r[k][0]:+9.1f} +- {r[k][1]:6.1f}" for k in ('C2', 'C3H'))
           + ('' if ok else '  NOT USED'))
        if ok:
            good.append(r)
    gx = np.array([r['x'] for r in good]); gb = np.array([max(abs(r[k][0]) + 3*r[k][1] for k in ('C2', 'C3H')) for r in good])
    sx = np.array([r['nu'] for r in good]); o = np.argsort(sx); sx = sx[o]
    sp = {k: CubicSpline(sx, np.array([good[i][k][0] for i in o]), bc_type='natural') for k in ('C2', 'C3H')}
    se = {k: np.array([good[i][k][1] for i in o]) for k in ('C2', 'C3H')}
    pr(f"#  nu_min = {sx.min():.2f}; node bounds |Gamma| + 3 sigma: " + ", ".join(f"{r['nu']:.0f}:{g:.0f}" for r, g in zip(good, gb)))
    for s in st:
        mfl = s['mfl']; di = s['direct']
        lamc = di['lam_central']; rest = s['E_tot'] - di['E_low']; rest_c = s['E_tot_c'] - di['E_low']
        hs, base, Ns, _ = H.base_of(1000.0, J, mfl)
        mst = mfl*Lm; x50 = SW.NU_W/mst
        t, w = leggauss(400); xx = 0.5*x50*(t + 1); wx = 0.5*x50*w
        w2, w3 = H.state_weights(1000.0, J, mfl, '+-', xx)
        nuq = xx*m; mid = xx >= gx.min()
        Dmid = float(((wx*w2*sp['C2'](nuq))[mid].sum() + (wx*w3*sp['C3H'](nuq))[mid].sum())/Ns[0])
        lin = float(((wx*w2*np.interp(nuq, sx, sp['C2'](sx)))[mid].sum() + (wx*w3*np.interp(nuq, sx, sp['C3H'](sx)))[mid].sum())/Ns[0])
        emc = 0.0
        for k, wk in (('C2', w2), ('C3H', w3)):
            for i in range(len(sx)):
                e_i = np.zeros(len(sx)); e_i[i] = 1.0
                emc += abs(float((wx*wk*np.interp(nuq, sx, e_i))[mid].sum()))/Ns[0]*se[k][i]
        E_lo = SW.cost(J, mfl, G_LO, 0.0, gx.min()*mst)
        oo = dict(mfl=mfl, lamc=lamc, rest=rest, rest_c=rest_c, E_low_ng=di['E_low'], E_lo=E_lo, D_mid=Dmid, E_signed=3*emc + 3*abs(Dmid - lin))
        for k in (1, 3, 10):
            oo[f'E_env{k}'] = SW.cost_x(J, mfl, SW.envelope(gx, gb, k, G_LO))
            oo[f'L_env{k}'] = lamc/(rest + oo[f'E_env{k}']); oo[f'Lc_env{k}'] = lamc/(rest_c + oo[f'E_env{k}'])
        oo['lam_signed'] = lamc + Dmid
        oo['L_signed'] = oo['lam_signed']/(rest + E_lo + oo['E_signed']); oo['Lc_signed'] = oo['lam_signed']/(rest_c + E_lo + oo['E_signed'])
        res['states'][f'{mfl:.7f}'] = oo
        pr(f"  m/Lam={mfl:.7f}: lam_c {lamc:+.4f}; other terms {rest:.3f} (interp on lam_c: {rest_c:.3f}); window 1e4 {di['E_low']:.3f} -> env x3 "
           f"{oo['E_env3']:.3f} (x10 {oo['E_env10']:.3f}); signed {Dmid:+.4f} +- {oo['E_signed']:.4f} (+ below nu_min {E_lo:.4f}) | L env3 "
           f"{oo['L_env3']:.2f} (interp on lam_c {oo['Lc_env3']:.2f}), env10 {oo['Lc_env10']:.2f}, signed {oo['Lc_signed']:.2f}")
    open(os.path.join(OUT, 'b5_window_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(OUT, 'b5_window_J%d.json' % J), 'w'), indent=1, default=lambda v: float(v) if isinstance(v, (np.floating, np.integer)) else str(v))


if __name__ == '__main__':
    J = int(sys.argv[1])
    if J % 2:
        odd(J)
    else:
        import sn_even as SE            # read-only import
        SE.TW[J] = C.twist(J)
        SE.main(J, [1000.0 + J])
