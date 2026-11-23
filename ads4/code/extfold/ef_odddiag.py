"""
ef_odddiag (EXT_FOLD.md §3.2-§3.3): diagnostics of hs_odd's budget terms on the direct-bulk rows J = 23, 27, 31 (and J = 21, SMALLNU's
row, for comparison when its points are present), with highspin/hs_odd's Model / value (read-only import):
  (a) representation: the degree-KB+1 shape fit that sets E_repr -- its own jackknife sigma (is the deviation Monte-Carlo noise
      amplified by the extra parameters, or a stable representation ambiguity?), and a model-free bulk (natural spline through the
      direct nodes, held beyond the last one) as an independent representation;
  (b) J = 31 edge (no GK15 window: 8 equispaced HC nodes 0.8-0.9, spacing 0.0125 against a period ~0.05): the published spline edge
      HOLDS Gamma(0.9) ~ 1.4e4 c2 over (0.9, Lam/m]; here that tail is set to 0 in the central value and covered sign-free with the
      light-cone envelope 3 A ((1 - x_l)/(1 - x))^2 of hs_odd's E_miss (A = max |Gamma| of the last 3 nodes, x_l = 0.9); the
      [0.8, 0.9] spline is checked against linear interpolation and the two every-other-node subsamples (x3: E_interp).
Usage: python3 ef_odddiag.py J [mfl,...] -> cache/ef_odddiag_J<J>.txt / .json
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C
import numpy as np
from numpy.polynomial.legendre import leggauss
H = C.H
import hs_odd as HO          # read-only import (functions only)
Dg = 1000


def edge_custom(J, mfl, M, f2, f3, xcut=0.8, nb=800):
    """lam/N_s with the bulk model of M on [0, xcut] and given edge functions on [xcut, 1/mfl]."""
    hs, base, Ns, st = H.base_of(Dg, J, mfl)
    t, w = leggauss(nb); xb = 0.5*xcut*(t + 1); wb = 0.5*xcut*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xb)
    d = (wb*w2*M.k['C2']['bulk'](xb)).sum() + (wb*w3*M.k['C3H']['bulk'](xb)).sum()
    xmax = 1.0/mfl
    t, w = leggauss(600); xx = xcut + 0.5*(xmax - xcut)*(t + 1); wx = 0.5*(xmax - xcut)*w
    w2e, w3e = H.state_weights(Dg, J, mfl, '+-', xx)
    d += (wx*w2e*f2(xx)).sum() + (wx*w3e*f3(xx)).sum()
    return float((base[0] + d)/Ns[0])


def tail_cost(J, mfl, A, xl, xcut=0.8):
    _, _, Ns, _ = H.base_of(Dg, J, mfl)
    xmax = 1.0/mfl
    if xmax <= xl:
        return 0.0
    t, w = leggauss(600); xx = xl + 0.5*(xmax - xl)*(t + 1); wx = 0.5*(xmax - xl)*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xx)
    g = 3*A*((1 - xl)/(1 - xx))**2
    return float((wx*g*(np.abs(w2) + np.abs(w3))).sum()/Ns[0])


if __name__ == '__main__':
    J = int(sys.argv[1]); D = 1000.0 + J
    rowm = H.HN.mphys(D, J); Lm = H.E.lam_of(Dg); mtw = rowm/Lm
    HO.XROW[0] = 1.0/mtw
    mfls = [mtw] + ([float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [])
    lines, res = [], dict(J=J, states={})

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    rK, rI, isdir = HO.rows(J, D)
    bad = H.refbad(J, D, HO.RTOL)
    reps = [HO.rows(J, D, jk=j)[1] for j in range(8)]
    dset = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
    isd = lambda rr: np.array([round(x, 6) in dset for x in rr['xs']])
    M = HO.Model(J, rK, rI, isdir, bad, bulk='dir')
    MK = {k: HO.Model(J, rK, rI, isdir, bad, bulk='dir', kb=k) for k in (4, 5, 6, 7, 8)}
    MKj = {k: [HO.Model(J, rK, rr, isd(rr), bad, bulk='dir', kb=k) for rr in reps] for k in (6, 7)}
    # model-free bulk: natural spline through the direct bulk nodes (x in [0.05, 0.8)), held outside
    sel = (rI['xs'] >= 0.05 - 1e-9) & (rI['xs'] < 0.8 - 1e-9) & isdir
    Ms = HO.Model(J, rK, rI, isdir, bad, bulk='dir')
    for kd in ('C2', 'C3H'):
        Ms.k[kd]['bulk'] = H.spline_of(rI['xs'][sel], np.asarray(rI[kd])[sel])
    gk_ok = len(M.ex) >= 10 and J != 31
    edge = 'gk' if gk_ok else 'espl'
    pr(f"# ef_odddiag J={J} row {D:g}: edge treatment {edge}; bulk direct nodes {int(sel.sum())}; shape-fit rms/rms by degree: " +
       ", ".join(f"K{k}: C2 {m.k['C2']['bulk'].rms_rel:.3f} C3H {m.k['C3H']['bulk'].rms_rel:.3f}" for k, m in MK.items()))
    for mfl in mfls:
        o = {}
        for k, m in MK.items():
            o[f'K{k}'] = HO.value(J, mfl, m, edge)[0]
        for k in (6, 7):
            lj = [HO.value(J, mfl, mj, edge)[0] for mj in MKj[k]]
            o[f'K{k}_jk'] = HO.jk_sigma(lj)
        o['spline_bulk'] = HO.value(J, mfl, Ms, edge)[0]
        pr(f"  m/Lam={mfl:.7f}: shape degree " + ", ".join(f"{k} {o[k]:+.4f}" for k in ('K4', 'K5', 'K6', 'K7', 'K8')) +
           f" | jk sigma K6 {o['K6_jk']:.4f}, K7 {o['K7_jk']:.4f} | (K7 - K6)/sigma_jk(K7) = {(o['K7'] - o['K6'])/o['K7_jk']:+.1f} | "
           f"model-free spline bulk {o['spline_bulk']:+.4f}")
        if J == 31:
            ex = M.ex; xl = float(ex.max())
            last = np.argsort(ex)[-3:]
            A = max(np.abs(M.k[kd]['ey'][last]).max() for kd in ('C2', 'C3H'))
            cut = lambda f: (lambda x: np.where(np.asarray(x) <= xl, f(x), 0.0))
            sp = {kd: M.k[kd]['espl'] for kd in ('C2', 'C3H')}
            lam_spl0 = edge_custom(J, mfl, M, cut(sp['C2']), cut(sp['C3H']))
            lin = {kd: (lambda kd_: (lambda x: np.interp(x, ex, M.k[kd_]['ey'])))(kd) for kd in ('C2', 'C3H')}
            lam_lin0 = edge_custom(J, mfl, M, cut(lin['C2']), cut(lin['C3H']))
            subs = []
            for s0 in (0, 1):
                o_ = np.argsort(ex); keep = o_[s0::2].tolist()
                for e_ in (o_[0], o_[-1]):
                    if e_ not in keep:
                        keep.append(e_)
                keep = np.array(sorted(keep, key=lambda i: ex[i]))
                fs = {kd: H.spline_of(ex[keep], M.k[kd]['ey'][keep]) for kd in ('C2', 'C3H')}
                subs.append(edge_custom(J, mfl, M, cut(fs['C2']), cut(fs['C3H'])))
            lam_held = HO.value(J, mfl, M, 'espl')[0]
            Etail = tail_cost(J, mfl, A, xl)
            Eint = 3*max(abs(lam_spl0 - lam_lin0), *[abs(lam_spl0 - v) for v in subs])
            o.update(lam_espl_held=lam_held, lam_espl_cut=lam_spl0, lam_lin_cut=lam_lin0, lam_sub=subs, E_tail=Etail, E_edge_interp=Eint, A=A, xl=xl)
            pr(f"        J=31 edge: spline held beyond x={xl:.3f} {lam_held:+.4f} -> tail set to 0: {lam_spl0:+.4f} (linear {lam_lin0:+.4f}; "
               f"subsamples {subs[0]:+.4f} / {subs[1]:+.4f}) | E_edge_interp (x3) {Eint:.3f} | tail (x > {xl:.2f}) sign-free, envelope "
               f"3 A ((1 - x_l)/(1 - x))^2, A = {A:.0f} c2: {Etail:.3f}")
        res['states'][f'{mfl:.7f}'] = o
    open(os.path.join(C.EFC, 'ef_odddiag_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.EFC, 'ef_odddiag_J%d.json' % J), 'w'), indent=1, default=float)
