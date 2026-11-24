"""
edger2/e2_gamma.py (EDGE_R2.md §2-§3): the graviton (EDGE-R) quantity E_chi[G_e] at Theorem W's gap-mass states from the batch-3 n1
rows (ext, partial), i.e. what is DETERMINED by them and what is not.

Gamma_{k,e}(x) (c2 units, kappa-curvature image subtracted) is produced by the project's own chain, run unchanged in-process
(hs_core -> neargap2/ng_val -> ng_chain -> w2_chain/ee1000 in EE mode -> q3/l9 grid, odd: p4_odd_analysis grid), exactly as HIGHSPIN /
EDGE_R do: the direct records (hc_nodes.to_record, M2 = M4 = 0) are injected over the row; the only redirections are
  * the direct point directory (hc_nodes.C) -> edger2/cache/pts (e2_merge.py),
  * the K2-record loader (w2_chain._pk) and the shadow candidates (hc_nodes._cands) also read edger2/cache/k2 (e2_k2.py): the K2 EE
    records that the chain's l8 helicity basis needs (state row x = 0.2-0.8; reference x_ref = 0.3-0.8 for J >= 34) and the <TTO>
    shadows b(Delta, J).  These K2 records are NOT row nodes (no reference at their nu); every node is a direct n1 node.
Per state (J, m = twist edge of the row Delta = 1000 + J) and entry e (++, +-, off for even J; +- for odd J):
  cost_k = int chi(p) w_k(x) Gamma_k(x) dx over p in [p_w, 1]   (er_core.edge_cost; w_k = hs_core.state_weights; off: ng_chain.dens
           convention, density psi_k env_off Gamma_off / m),
  bulk [x_w, a): degree-KB shape fit of the direct 0.025-grid nodes (hs_odd 'dir' bulk); edge [a, x1]: the job list's two GK15
  windows (er_jobs.state_nodes), K15 weights on the measured nodes ("determined"), the rest = REMAINDER (n3-only nodes not run,
  and any failed node).
Budget of the determined part: mc = 3 x jackknife sigma (8 RQMC replicas, everything rebuilt), repr = 3 x max|KB -> KB-2, KB-1, KB+1|,
quad = 3 x |K15 - G7| on the fully measured lower window.  Remainder: sign-free, per node |K15 w| chi |w_k| G_env(x) with
G_env = c A ((1 - x_l)/(1 - x))^2 (A = max|Gamma| of the last three valid nodes of the key, x_l the last valid node; c from e2_calib).
Usage: python3 e2_gamma.py J [KB]   -> cache/gamma_J<J>.json / .txt
"""
import e2_guard  # noqa: F401  (first; registers itself as er_guard / hs_guard / ng_guard)
import os, sys, json, glob, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(CODE, 'edger')); sys.path.insert(0, os.path.join(CODE, 'highspin'))
import numpy as np
import er_core as R            # read-only (its er_guard import resolves to e2_guard)
import hs_core as HC           # read-only
import er_jobs as EJ           # read-only: the job-list node layout
CACHE = os.path.join(HERE, 'cache')
PTS = os.path.join(CACHE, 'pts')
K2D = os.path.join(CACHE, 'k2')
NC, HN, V = HC.NC, HC.HN, HC.V
W = NC.W

# ---------------------------------------------------------------------------------------------------------- redirections
HN.C = PTS
HN._SC.clear()
_pk0 = W._pk


def _pk(pats):
    Rr = _pk0(pats)
    for p in pats:
        for fn in sorted(glob.glob(os.path.join(K2D, p))):
            Rr.update(pickle.load(open(fn, 'rb')))
    return Rr


W._pk = _pk
_c1 = HN._cands


def _cands(odd):
    d = _c1(odd)
    key = 'odd_e2' if odd else 'even_e2'
    if key not in HN._SC:
        dd = {k: list(v) for k, v in d.items()}
        f = os.path.join(K2D, 'w1_p2shadow_e2.pkl' if odd else 'w1_shadow_k2_e2.pkl')
        if os.path.exists(f):
            for k, v in pickle.load(open(f, 'rb')).items():
                dd.setdefault(k, []).append((os.path.basename(f), v))
        HN._SC[key] = dd
    return HN._SC[key]


HN._cands = _cands
FAILED = set()          # x (rounded 6) of nodes treated as missing
_rec0 = HN.records


def _records(kind, J, D, jk=None, **kw):
    out = _rec0(kind, J, D, jk, **kw)
    m = HN.mphys(D, J)
    return {k: v for k, v in out.items() if round(k[1]/m, 6) not in FAILED}


HN.records = _records
_clean0 = V._clean


def _clean(g, J):
    """the off-diagonal Gamma_off of these rows is NaN in the chain (l9 grid: sqrt(a_++ a_+-) with a_++ a_+- < 0 -- the helicity basis
    signs of l8 basis' realvec are arbitrary), which would make ng_val._clean drop EVERY node; the off keys are removed first (the
    off-diagonal entry is not evaluated here, as in EDGE_R §1) and the diagonal entries are cleaned as published."""
    for D, r in g.items():
        for k in [k for k in r if isinstance(k, tuple) and len(k) >= 2 and k[1] == 'off']:
            OFFNAN.append((J, D, int(np.isfinite(np.asarray(r[k], complex)).sum())))
            r.pop(k)
    return _clean0(g, J)


OFFNAN = []
V._clean = _clean
Dg = 1000
TAGS = {26: 1.026, 30: 1.030, 34: 1.033, 33: 1.033, 38: 1.037, 42: 1.041, 46: 1.045, 50: 1.049}


def keys(J):
    return [('C2', '+-'), ('C3H', '+-')] if J % 2 else [(kd, h) for h in ('++', '+-') for kd in ('C2', 'C3H')]


def gkey(J, key):
    return key[0] if J % 2 else key


def node_relerr(J, D):
    """per direct node: max over kinds of the Frobenius relative replica error ||err||/||E|| (state, reference) and the max-entry one."""
    m = HN.mphys(D, J); out = {}
    for kind in ('EE', 'C3H'):
        P = HN.pts(kind, J, D)
        for (DD, nu), r in P.items():
            if DD != D or (20000.0, nu) not in P:
                continue
            rf = P[(20000.0, nu)]
            fro = max(np.linalg.norm(r['err'])/np.linalg.norm(r['E']), np.linalg.norm(rf['err'])/np.linalg.norm(rf['E']))
            ent = max(float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300))),
                      float(np.max(np.abs(rf['err'])/np.maximum(np.abs(rf['E']), 1e-300))))
            x = round(nu/m, 6); o = out.setdefault(x, [0.0, 0.0])
            o[0] = max(o[0], fro); o[1] = max(o[1], ent)
    return out


def row(J, D, jk=None):
    g = HC.grid(J, 2, inject=[D], jk=jk)
    if D not in g:
        raise RuntimeError(('row not built', J, D, sorted(g)))
    return HC.subtract_kappa(g, J)[D]


class Model:
    def __init__(self, J, r, a, KB):
        xs = np.asarray(r['xs'], float)
        self.xs = xs
        bs = (xs >= 0.0) & (xs < a - 1e-9)
        self.bulk = {}; self.val = {}
        for key in keys(J):
            y = np.asarray(r[gkey(J, key)], float)
            self.val[key] = y
            h = key[1]
            if h == 'off':
                self.bulk[key] = HC.spline_of(xs[bs], y[bs]) if KB >= 0 else None
            else:
                self.bulk[key] = HC.ShapeFit(J, xs[bs], y[bs], h, K=KB)
        self.nb = int(bs.sum())

    def at(self, key, x):
        i = int(np.argmin(np.abs(self.xs - x)))
        return float(self.val[key][i]) if abs(self.xs[i] - x) < 5e-6 else None


def weights(J, mfl, h, x):
    """(w2, w3) per unit Gamma, cost units (er_core.weights); off-diagonal in ng_chain.dens' convention (density psi env_off Gamma/m)."""
    if h != 'off':
        return R.weights(J, mfl, h, x)
    w2, w3 = HC.state_weights(Dg, J, mfl, 'off', x)
    m = mfl*R.lam()
    return w2*m, w3*m


def costs(J, mfl, M, lay, pe, sg, edge='k15'):
    """per key: (bulk, edge_determined) signed costs, and the per-node remainder weights |K15 w| chi |w_k| at the missing edge nodes."""
    L_ = R.lam(); m = mfl*L_
    pw = 3.0/(R.beta_of(J, mfl)*L_); xw = pw*L_/m
    xb, wb = R.gl(xw, lay['a'], 800)
    out = {}
    for key in keys(J):
        kd, h = key; k = 0 if kd == 'C2' else 1
        gb = M.bulk[key](xb)
        wbk = weights(J, mfl, h, xb)[k]
        pb = xb*mfl
        cb = float((wb*R.chi(pb, pe, sg)*((pb >= pw) & (pb <= 1))*wbk*gb).sum())
        ce = 0.0; miss = []; ce7 = 0.0
        for (lo, hi) in lay['win']:
            xe, wk, wg = HC.gk15_nodes(lo, hi)
            pe_ = xe*mfl
            we = weights(J, mfl, h, xe)[k]*R.chi(pe_, pe, sg)*(pe_ <= 1)
            for xi, wki, wgi, wei in zip(xe, wk, wg, we):
                v = M.at(key, xi)
                if v is None:
                    miss.append((float(xi), float(abs(wki*wei)), float(abs(wgi*wei))))
                else:
                    ce += wki*wei*v; ce7 += wgi*wei*v
        out[key] = dict(bulk=cb, edge=ce, edge_g7=ce7, det=cb + ce, miss=miss)
    return out


def layout(J, mfl):
    wr = R.w_state(J, mfl)
    s = EJ.state_nodes(J, wr['mfl'], wr['pe'], wr['sigma'], wr['pw'])
    a, x1 = s['edge_a'], s['x1']; xm = 0.5*(a + x1)
    return dict(a=a, x1=x1, xm=xm, win=[(a, xm), (xm, x1)], x0=s['x0'], wr=wr, D=s['D'])


def main(J, KB=6, rtol=0.2):
    t0 = time.time(); OUT = []
    pr = lambda s='': (print(s, flush=True), OUT.append(s))
    mfl0 = TAGS[J]
    lay = layout(J, mfl0); wr = lay['wr']; mfl = wr['mfl']; D = 1000.0 + J
    assert abs(lay['D'] - D) < 1e-6, lay['D']
    L_ = R.lam()
    rel = node_relerr(J, D)
    FAILED.clear(); FAILED.update(x for x, (fro, ent) in rel.items() if fro > rtol)
    r0 = row(J, D)
    M = Model(J, r0, lay['a'], KB)
    reps = [Model(J, row(J, D, jk=j), lay['a'], KB) for j in range(8)]
    Mr = {kb: Model(J, r0, lay['a'], kb) for kb in (KB - 2, KB - 1, KB + 1)}
    SE, SEt, _ = R.SE_units(J, mfl)
    H = wr['H']; tol = wr['tolE_tr']*SEt
    pe, sg = wr['pe'], wr['sigma']
    c0 = costs(J, mfl, M, lay, pe, sg)
    cj = [costs(J, mfl, mj, lay, pe, sg) for mj in reps]
    cr = {kb: costs(J, mfl, mm, lay, pe, sg) for kb, mm in Mr.items()}
    xv = np.array(sorted(x for x in M.xs))
    pr(f"J = {J}, row Delta = {D:g}, m/Lam = {mfl:.5f}, x1 = {lay['x1']:.5f}; cutoff p_e = {pe}, sigma = {sg:.4f}; edge windows "
       f"[{lay['a']:.4f}, {lay['xm']:.4f}], [{lay['xm']:.4f}, {lay['x1']:.4f}]; reference zero x0 = {lay['x0']:.4f}")
    pr(f"  direct nodes in the row: {len(xv)} (x {xv.min():.4f} .. {xv.max():.4f}); bulk fit nodes {M.nb}; failed (Frobenius rel. err > "
       f"{rtol}): {sorted(FAILED)}; max-entry rel. err > {rtol} at {sorted(x for x, v in rel.items() if v[1] > rtol)}")
    pr(f"  H = {H:.4e}, S_E = {SE:.4e}, S_E,t = {SEt:.4e} (chi.json S_E,t {wr['SEt']:.4e}); tolerance {wr['tolE_tr']:.3f} S_E,t = {tol:.4e}")
    res = dict(J=J, D=D, mfl=mfl, x1=lay['x1'], a=lay['a'], xm=lay['xm'], pe=pe, sigma=sg, H=H, SE=SE, SEt=SEt, tol=tol,
               tol_tr=wr['tolE_tr'], failed=sorted(FAILED), entries={}, nodes={})
    # node table near the edge
    for key in keys(J):
        ys = M.val[key]; js = np.array([mj.val[key] for mj in reps])
        sj = np.sqrt(7/8*((js - js.mean(0))**2).sum(0))
        res['nodes'][str(key)] = dict(x=M.xs.tolist(), G=ys.tolist(), jk=sj.tolist())
        sel = M.xs >= 0.85
        pr(f"  {key}: Gamma at x >= 0.85: " + ", ".join(f"{x:.4f}: {y:+.3g}+-{e:.2g}" for x, y, e in zip(M.xs[sel], ys[sel], sj[sel])))
    for h in (('+-',) if J % 2 else ('++', '+-')):
        ent = {}
        for kd in ('C2', 'C3H'):
            key = (kd, h); c = c0[key]
            js = np.array([q[key]['det'] for q in cj]); mc = 3*float(np.sqrt(7/8*((js - js.mean())**2).sum()))
            rep = 3*max(abs(c['det'] - q[key]['det']) for q in cr.values())
            quad = 3*abs(c['edge'] - c['edge_g7']) if not c['miss'] else float('nan')
            # lower window only (fully measured): K15 - G7 there
            lo, hi = lay['win'][0]
            xe, wk, wg = HC.gk15_nodes(lo, hi)
            we = weights(J, mfl, h, xe)[0 if kd == 'C2' else 1]*R.chi(xe*mfl, pe, sg)
            vals = np.array([M.at(key, x) if M.at(key, x) is not None else np.nan for x in xe])
            quad_lo = 3*abs(np.nansum((wk - wg)*we*vals)) if np.all(np.isfinite(vals)) else float('nan')
            # last valid nodes for the envelope
            ok = np.array([M.at(key, x) is not None for x in M.xs])
            xs_ok = M.xs[ok]; order = np.argsort(xs_ok)[-3:]
            A = float(np.max(np.abs(M.val[key][ok][order]))); xl = float(xs_ok.max())
            # same with A over both kinds (HIGHSPIN's rule)
            ent[kd] = dict(det=c['det'], bulk=c['bulk'], edge=c['edge'], mc=mc, repr=rep, quad_lo=quad_lo, miss=c['miss'], A=A, xl=xl)
        ent['A_both'] = max(ent['C2']['A'], ent['C3H']['A'])
        res['entries'][h] = ent
        u = lambda v: v/SEt
        pr(f"  entry {h}: determined cost_2 = {u(ent['C2']['det']):+.4f} (bulk {u(ent['C2']['bulk']):+.4f}, edge {u(ent['C2']['edge']):+.4f}) "
           f"+- mc {u(ent['C2']['mc']):.4f}, repr {u(ent['C2']['repr']):.4f}, quad(lower window) {u(ent['C2']['quad_lo']):.4f} | "
           f"cost_3 = {u(ent['C3H']['det']):+.4f} (bulk {u(ent['C3H']['bulk']):+.4f}, edge {u(ent['C3H']['edge']):+.4f}) +- mc "
           f"{u(ent['C3H']['mc']):.4f}, repr {u(ent['C3H']['repr']):.4f}, quad {u(ent['C3H']['quad_lo']):.4f}   [S_E,t]")
        for kd in ('C2', 'C3H'):
            e = ent[kd]
            w1 = sum(q[1] for q in e['miss'])
            pr(f"      {kd}: missing edge nodes {len(e['miss'])} (x {', '.join(f'{q[0]:.4f}' for q in e['miss'])}); "
               f"sum |K15 w chi w_k| = {w1:.3e} per c2; A = {e['A']:.3g} c2 (last 3 valid, x_l = {e['xl']:.4f}); "
               f"x1 sign-free with G = A ((1-x_l)/(1-x))^2: {sum(q[1]*e['A']*max(1, ((1 - e['xl'])/(1 - q[0]))**2) for q in e['miss'])/SEt:.4f} S_E,t")
    res['offnan'] = OFFNAN[:3]
    res['a'] = {str(k): float(v) for k, v in r0['a'].items()} if isinstance(r0['a'], dict) else float(r0['a'])
    pr(f"  normalisation a_h: {res['a']}; off-diagonal finite values (of n) in the chain's grid: {OFFNAN[:1]}")
    res['t'] = time.time() - t0
    pr(f'[{res["t"]:.0f}s]')
    json.dump(res, open(os.path.join(CACHE, 'gamma_J%d.json' % J), 'w'), default=float, indent=1)
    open(os.path.join(CACHE, 'gamma_J%d.txt' % J), 'w').write('\n'.join(OUT) + '\n')
    return res


if __name__ == '__main__':
    main(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 6)
