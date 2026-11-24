"""
edger/er_row.py (EDGE_R.md §3, §5): the graviton (EDGE-R) quantity on an odd-J own-mass row Delta = 1000 + J with the direct n1 points of
HIGHSPIN / EXT_FOLD (highspin/cache, read-only) and the n3 points of this note (edger/cache, er_points.py) MERGED (n3 fills the GK15 edge
nodes where n1 fails, nu/Delta >~ 0.94; at the overlap nodes n1 is kept and n3 is the cross-check).  Merged point files are written to
edger/cache/merged/ (copies; nothing in highspin/ is touched) and injected through hs_core exactly as HIGHSPIN does.

Gamma(x) as hs_odd (kappa subtracted; bulk [0.05, 0.8): degree-KB shape fit of the direct nodes; edge [0.8, Lam/m]: GK15 quadrature of
the direct values, no interpolation).  Per cutoff (p_e, sigma): cost_E = |cost_2| + |cost_3| over [p_w, 1] (er_core), and the budget
    mc    : 3 x jackknife sigma over the 8 RQMC replicas (all injected nodes, n1 and n3),
    quad  : 3 x |K15 - G7| on the edge window,
    repr  : 3 x max |cost_E(KB) - cost_E(KB +- 1, KB - 2)| (bulk shape degree),
    zero  : the reference-zero GK node (HIGHSPIN §2) covered sign-free by 3 x max|Gamma| of its 4 nearest valid edge nodes,
    miss  : any edge node still without a valid value, covered sign-free by 3 A ((1 - x_l)/(1 - x))^2 (HIGHSPIN),
in units of H, S_E (plain) and S_E,t (+- transferred).
Usage: python3 er_row.py J [KB]   -> cache/row_J<J>.txt / .json
"""
import er_core as R
import os, sys, json, time, pickle, glob
import numpy as np
import hs_core as HC
import hs_odd as HO

MERGED = os.path.join(R.CACHE, 'merged')
OUT = []
pr = lambda s='': (print(s, flush=True), OUT.append(s))
GRID_PE = (0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.75)
GRID_T = (1.0, 0.8, 0.65, 0.5)
GRID_S = (0.35, 0.45, 0.6)


def relerr(r):
    return float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300)))


def merge(J, D, use_n3_overlap=False):
    """highspin n1 points + edger n3 points -> merged/ (n3 replaces a node only where n1 failed, or everywhere if use_n3_overlap)."""
    os.makedirs(MERGED, exist_ok=True)
    m = HC.HN.mphys(D, J); info = {}
    for kind in ('EE', 'C3H'):
        f = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
        a = pickle.load(open(os.path.join(R.CODE, 'highspin', 'cache', f), 'rb'))
        bp = os.path.join(R.CACHE, f)
        b = pickle.load(open(bp, 'rb')) if os.path.exists(bp) and not os.environ.get('ER_NO_N3') else {}   # ER_NO_N3=1: n1 only
        out = dict(a)
        for (DD, nu), r in b.items():
            if DD != D or (20000.0, nu) not in b:
                continue
            old = a.get((DD, nu)); oldr = a.get((20000.0, nu))
            n1bad = old is None or oldr is None or relerr(old) > 0.2 or relerr(oldr) > 0.2
            if n1bad or use_n3_overlap:
                for kk in [k for k in out if abs(k[1]/m - nu/m) < 5e-6 and k[1] != nu]:   # n1 record at the same (rounded) x
                    del out[kk]
                out[(DD, nu)] = r; out[(20000.0, nu)] = b[(20000.0, nu)]
                info[round(nu/m, 6)] = 'n3' + ('' if n1bad else ' (overlap)')
        pickle.dump(out, open(os.path.join(MERGED, f), 'wb'))
    return info


def setup(J, D, use_n3_overlap=False):
    info = merge(J, D, use_n3_overlap)
    HC.HN.C = MERGED; HC.HN._SC.clear()
    return info


def costs(J, mfl, M, xrow, cuts, kb_models=None, edge='gk'):
    L_ = R.lam(); m = mfl*L_
    pw = 3.0/(R.beta_of(J, mfl)*L_); xw = pw*L_/m
    xb, wb = R.gl(xw, 0.8, 600)
    xe, wk, wg = HC.gk15_nodes(0.8, xrow)
    ww = wk if edge == 'gk' else wg
    have = np.array([np.any(np.abs(M.ex - x) < 5e-6) for x in xe])
    x0 = HC.ref_zero_x(J, 1000.0 + J)
    inz = np.abs(xe - x0) < HO.ZWIN
    okx = np.abs(M.ex - x0) >= HO.ZWIN
    near = np.argsort(np.abs(M.ex - x0) + 1e9*(~okx))[:4]
    gz = 3*max(np.abs(M.k[k_]['ey'][near]).max() for k_ in ('C2', 'C3H'))
    p = xe*m/L_
    w2e, w3e = R.weights(J, mfl, '+-', xe)
    out = []
    for (pe, sg) in cuts:
        rec = {}
        for kd, k in (('C2', 0), ('C3H', 1)):
            gb = M.k[kd]['bulk'](xb)
            ge = np.array([M.k[kd]['ey'][np.argmin(np.abs(M.ex - x))] if h else 0.0 for x, h in zip(xe, have)])
            ge[inz] = 0.0
            cb = R.edge_cost(J, mfl, '+-', xb, wb, gb, gb, pe, sg)[k]
            ce = R.edge_cost(J, mfl, '+-', xe, ww, ge, ge, pe, sg)[k]
            wkk = (w2e if k == 0 else w3e)*R.chi(p, pe, sg)*(p <= 1)*ww
            zero = float((gz*np.abs(wkk))[inz].sum())
            miss = 0.0
            if (~have & ~inz).any():
                last = np.argsort(M.ex + -9*(~okx))[-3:]
                A = max(np.abs(M.k[k_]['ey'][last]).max() for k_ in ('C2', 'C3H')); xl = M.ex[last].max()
                genv = 3*A*np.maximum(1.0, ((1 - xl)/(1 - xe))**2)
                miss = float((genv*np.abs(wkk))[~have & ~inz].sum())
            rec[kd] = dict(bulk=cb, edge=ce, val=cb + ce, zero=zero, miss=miss)
        out.append(rec)
    return out


def tot(c):
    return abs(c['C2']['val']) + abs(c['C3H']['val'])


def main(J, KB=6):
    t0 = time.time()
    pr(__doc__.split('Usage')[0])
    D = 1000.0 + J
    rowm = HC.HN.mphys(D, J); L_ = R.lam(); mtw = rowm/L_; xrow = 1/mtw
    HO.XROW[0] = xrow
    res = dict(J=J, D=D, mfl=mtw)
    # (a) n1 + n3 where n1 failed;  (b) n3 also at the overlap nodes
    tabs = {}
    for var in ('a', 'b'):
        info = setup(J, D, use_n3_overlap=(var == 'b'))
        rK, rI, isdir = HO.rows(J, D)
        tabs[var] = (rI, isdir, info)
    rIa, isdir, info = tabs['a'][0], tabs['a'][1], tabs['a'][2]
    rIb = tabs['b'][0]
    pr(f"J = {J}, row {D:g}, twist edge m/Lam = {mtw:.5f}, x1 = {xrow:.5f}; n3 nodes: {info}; overlap n3: "
       f"{ {k: v for k, v in tabs['b'][2].items() if 'overlap' in v} }")
    pr("  edge node table (Gamma - kappa, c2): x | C2 n1-or-merged | C2 with n3 | C3H merged | C3H with n3")
    for i, x in enumerate(rIa['xs']):
        if x < 0.895:
            continue
        j = int(np.argmin(np.abs(rIb['xs'] - x)))
        pr(f"   x = {x:.6f}: C2 {rIa['C2'][i]:+12.1f} | {rIb['C2'][j]:+12.1f} | C3H {rIa['C3H'][i]:+12.1f} | {rIb['C3H'][j]:+12.1f}"
           f"   {'direct' if isdir[i] else 'K2'} {info.get(round(x, 6), '')}")
    # main evaluation on variant (a)
    setup(J, D, use_n3_overlap=False)
    rK, rI, isdir = HO.rows(J, D)
    bad = HC.refbad(J, D, HO.RTOL)
    fx = HC.failed_x(J, D)
    M = HO.Model(J, rK, rI, isdir, bad, bulk='dir', kb=KB)
    dset = set(np.round(HC.direct_x(J, D), 6)) - fx
    Mj = []
    for j in range(8):
        rr = HO.rows(J, D, jk=j)[1]
        Mj.append(HO.Model(J, rK, rr, np.array([round(x, 6) in dset for x in rr['xs']]), bad, bulk='dir', kb=KB))
    Mr = {k: HO.Model(J, rK, rI, isdir, bad, bulk='dir', kb=k) for k in (KB - 2, KB - 1, KB + 1)}
    H, f = R.H_of(J, mtw)
    SE, SEt, _ = R.SE_units(J, mtw)
    beta = R.beta_of(J, mtw); T0 = 25/(0.97*beta*L_)
    kl = lambda pe: (2*J + 1)/mtw/np.sqrt(1 - (pe/mtw)**2)
    wr = R.nearest_cut(J, mtw)
    cuts = [(wr['pe'], wr['sigma'])] + [(pe, sf*tf*T0) for pe in GRID_PE for tf in GRID_T for sf in GRID_S if kl(pe)*sf*tf*T0 >= 10]
    c0 = costs(J, mtw, M, xrow, cuts)
    c7 = costs(J, mtw, M, xrow, cuts, edge='g7')
    cj = [costs(J, mtw, mj, xrow, cuts) for mj in Mj]
    cr = {k: costs(J, mtw, mm, xrow, cuts) for k, mm in Mr.items()}
    pr(f"  failed n1 nodes (rel err > 0.2) after merging: {sorted(fx)}; H = {H:.3e}, S_E = {SE:.3e}, S_E,t = {SEt:.3e} (S_E,t/H = {SEt/H:.3f})")
    pr(f"  nearest W cutoff: J = {wr['J']}, m/L = {wr['mfl']:.3f}: p_e = {wr['pe']}, sigma = {wr['sigma']:.3f}; its tolerance "
       f"{wr['tolE_tr']:.3g} S_E,t (even-J H; this odd state's own H gives (1 - d_chi) H / S_E,t = {(1 - wr['d_tr'])*H/SEt:.3g})")
    rows = []
    for i, (pe, sg) in enumerate(cuts):
        v = tot(c0[i])
        js = np.array([tot(c[i]) for c in cj]); mc = 3*float(np.sqrt(7/8*((js - js.mean())**2).sum()))
        quad = 3*abs(v - tot(c7[i]))
        rep = 3*max(abs(v - tot(c[i])) for c in cr.values())
        zero = c0[i]['C2']['zero'] + c0[i]['C3H']['zero']; miss = c0[i]['C2']['miss'] + c0[i]['C3H']['miss']
        b = mc + quad + rep + zero + miss
        rows.append(dict(pe=pe, sigma=sg, cost=v, cost2=c0[i]['C2']['val'], cost3=c0[i]['C3H']['val'], mc=mc, quad=quad, repr=rep,
                         zero=zero, miss=miss, upper=v + b, edge2=c0[i]['C2']['edge'], bulk2=c0[i]['C2']['bulk']))
    r = rows[0]
    u = lambda v: f"{v/SEt:.3f}"
    pr(f"  own-W cutoff: cost_2 = {r['cost2']:+.3e} (bulk {r['bulk2']:+.2e}, edge {r['edge2']:+.2e}), cost_3 = {r['cost3']:+.3e}; "
       f"cost_E = {u(r['cost'])} S_E,t; budget: mc {u(r['mc'])}, quad {u(r['quad'])}, repr {u(r['repr'])}, zero-node {u(r['zero'])}, "
       f"missing {u(r['miss'])} -> upper {u(r['upper'])} S_E,t = {r['upper']/H:.3f} H = {r['upper']/SE:.1f} S_E")
    up = np.array([q['upper'] for q in rows[1:]])/SEt; cs = np.array([q['cost'] for q in rows[1:]])/SEt
    pr(f"  AG grid ({len(rows) - 1} cutoffs): cost_E/S_E,t median {np.median(cs):.3f} max {cs.max():.3f}; upper (cost + budget) median "
       f"{np.median(up):.3f} max {up.max():.3f}")
    res.update(H=H, SE=SE, SEt=SEt, rows=rows, info={str(k): v for k, v in info.items()}, wcut=dict(J=wr['J'], mfl=wr['mfl']))
    pr(f'[{time.time() - t0:.0f}s]')
    json.dump(res, open(os.path.join(R.CACHE, 'row_J%d.json' % J), 'w'), default=float, indent=1)
    open(os.path.join(R.CACHE, 'row_J%d.txt' % J), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 6)
