"""
edger/er_existing.py (EDGE_R.md §3): the graviton (EDGE-R) quantity on the EXISTING direct odd-J near-gap rows (HIGHSPIN / EXT_FOLD:
J = 21, 23, 27, 29 own-mass rows Delta = 1000 + J: direct bulk row + GK15 direct edge on [0.8, Lam/m]; RQMC fails at nu/Delta >~ 0.94).

Gamma(x) exactly as hs_odd (kappa subtracted; bulk x in [0.05, 0.8): degree-6 shape fit of the direct nodes; edge: GK15 quadrature of
the direct values, no interpolation).  For each state (twist edge of the row) and each cutoff (p_e, sigma) of AG_PROOF's grid with
kappa(p_e) sigma >= 10 (and the nearest Theorem W state's own cutoff), reports
    cost_E = |cost_2| + |cost_3| over [p_w, 1]  with the failed GK nodes (x >~ 0.94)  (a) set to 0,  (b) covered sign-free (and the reference-zero GK node, HIGHSPIN §2, by 3 x max|Gamma| of its 4 nearest valid nodes) by the
    HIGHSPIN envelope 3 A ((1 - x_l)/(1 - x))^2,
in units of H (P4W), S_E (plain) and S_E,t (+- transferred); and the jackknife (8 RQMC replicas) sigma of cost_E.
Usage: python3 er_existing.py J[,J...]   -> cache/existing.txt, cache/existing.json
"""
import er_core as R
import os, sys, json, time
import numpy as np
import hs_core as HC
import hs_odd as HO

OUT = []
pr = lambda s='': (print(s, flush=True), OUT.append(s))
GRID_PE = (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)
GRID_T = (1.0, 0.8, 0.65, 0.5)
GRID_S = (0.35, 0.45, 0.6)


def model_costs(J, mfl, M, xrow, cuts, miss_mode):
    L_ = R.lam(); m = mfl*L_
    pw = 3.0/(R.beta_of(J, mfl)*L_); xw = pw*L_/m
    xb, wb = R.gl(xw, 0.8, 600)
    xe, wk, _ = HC.gk15_nodes(0.8, xrow)
    have = np.array([np.any(np.abs(M.ex - x) < 5e-6) for x in xe])
    x0 = HC.ref_zero_x(J, 1000.0 + J)
    inz = np.abs(xe - x0) < HO.ZWIN            # reference-zero node: direct ratio unusable (HIGHSPIN §2) -> 0 + sign-free cover
    okx = np.abs(M.ex - x0) >= HO.ZWIN
    near = np.argsort(np.abs(M.ex - x0) + 1e9*(~okx))[:4]
    gz = 3*max(np.abs(M.k[k_]['ey'][near]).max() for k_ in ('C2', 'C3H'))
    out = []
    for (pe, sg) in cuts:
        c = []
        for kd in ('C2', 'C3H'):
            gb = M.k[kd]['bulk'](xb)
            ge = np.array([M.k[kd]['ey'][np.argmin(np.abs(M.ex - x))] if h else 0.0 for x, h in zip(xe, have)])
            ge[inz] = 0.0
            cb = R.edge_cost(J, mfl, '+-', xb, wb, gb, gb, pe, sg)
            ce = R.edge_cost(J, mfl, '+-', xe, wk, ge, ge, pe, sg)
            k = 0 if kd == 'C2' else 1
            val = cb[k] + ce[k]
            sf = 0.0
            if miss_mode == 'cover':
                w2, w3 = R.weights(J, mfl, '+-', xe)
                p = xe*m/L_
                wkk = (w2 if k == 0 else w3)*R.chi(p, pe, sg)*(p <= 1)*wk
                sf += float((gz*np.abs(wkk))[inz].sum())
                last = np.argsort(M.ex + -9*(~okx))[-3:]
                A = max(np.abs(M.k[k_]['ey'][last]).max() for k_ in ('C2', 'C3H')); xl = M.ex[last].max()
                genv = 3*A*np.maximum(1.0, ((1 - xl)/(1 - xe))**2)
                sf += float((genv*np.abs(wkk))[~have & ~inz].sum())
            c.append((val, sf))
        out.append(c)
    return out


def main(Js):
    t0 = time.time()
    pr(__doc__.split('Usage')[0])
    res = {}
    for J in Js:
        D = 1000.0 + J
        rK, rI, isdir = HO.rows(J, D)
        bad = HC.refbad(J, D, HO.RTOL)
        rowm = HC.HN.mphys(D, J); L_ = R.lam(); mtw = rowm/L_; xrow = 1/mtw
        M = HO.Model(J, rK, rI, isdir, bad, bulk='dir')
        Mj = []
        dset = set(np.round(HC.direct_x(J, D), 6)) - HC.failed_x(J, D)
        for j in range(8):
            rr = HO.rows(J, D, jk=j)[1]
            Mj.append(HO.Model(J, rK, rr, np.array([round(x, 6) in dset for x in rr['xs']]), bad, bulk='dir'))
        H, f = R.H_of(J, mtw)
        SE, SEt, _ = R.SE_units(J, mtw)
        beta = R.beta_of(J, mtw); T0 = 25/(0.97*beta*L_)
        kl = lambda pe: (2*J + 1)/mtw/np.sqrt(1 - (pe/mtw)**2)
        cuts = [(pe, sf*tf*T0) for pe in GRID_PE for tf in GRID_T for sf in GRID_S if kl(pe)*sf*tf*T0 >= 10]
        wr = R.nearest_cut(J, mtw)
        cuts = [(wr['pe'], wr['sigma'])] + cuts
        c0 = model_costs(J, mtw, M, xrow, cuts, 'cover')
        cj = [model_costs(J, mtw, mj, xrow, cuts[:1], 'zero') for mj in Mj]
        pr(f"J = {J}, row {D:g}, twist edge m/Lam = {mtw:.5f} (x1 = {xrow:.4f}); H = {H:.3e}; S_E = {SE:.3e}, S_E,t = {SEt:.3e}; "
           f"failed GK nodes: {sorted(HC.failed_x(J, D))}; nearest W cutoff (J = {wr['J']}, m/L = {wr['mfl']:.3f}): p_e = {wr['pe']}, "
           f"sigma = {wr['sigma']:.3f}")
        rows = []
        for (pe, sg), c in zip(cuts, c0):
            e2, e3 = c[0][0], c[1][0]; s2, s3 = c[0][1], c[1][1]
            cost = abs(e2) + abs(e3); cov = s2 + s3
            rows.append(dict(pe=pe, sigma=sg, cost2=e2, cost3=e3, cost=cost, cover=cov))
        v = np.array([[abs(cj_[0][0][0]) + abs(cj_[0][1][0]), cj_[0][0][0], cj_[0][1][0]] for cj_ in cj])
        sj = float(np.sqrt(7/8*((v[:, 0] - v[:, 0].mean())**2).sum()))
        r0 = rows[0]
        pr(f"   own-W cutoff: cost_2 = {r0['cost2']:+.3e}, cost_3 = {r0['cost3']:+.3e} -> |.|+|.| = {r0['cost']:.3e} = {r0['cost']/H:.3e} H "
           f"= {r0['cost']/SE:.3g} S_E = {r0['cost']/SEt:.3g} S_E,t (jk sigma {sj/SEt:.2g} S_E,t); failed nodes sign-free {r0['cover']/SEt:.3g} S_E,t "
           f"({r0['cover']/H:.2e} H)")
        cs = np.array([r['cost'] for r in rows[1:]])/SEt; cv = np.array([r['cover'] for r in rows[1:]])/SEt
        pr(f"   AG grid ({len(rows) - 1} cutoffs): cost_E/S_E,t min {cs.min():.3g} median {np.median(cs):.3g} max {cs.max():.3g}; "
           f"sign-free failed-node cover max {cv.max():.3g} S_E,t, min {cv.min():.3g}")
        res[J] = dict(J=J, mfl=mtw, H=H, f=f, SE=SE, SEt=SEt, rows=rows, jk=sj, failed=sorted(HC.failed_x(J, D)))
    pr(f'[{time.time() - t0:.0f}s]')
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(R.CACHE, 'existing.json'), 'w'), default=float, indent=1)
    open(os.path.join(R.CACHE, 'existing.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main([int(v) for v in sys.argv[1].split(',')])
