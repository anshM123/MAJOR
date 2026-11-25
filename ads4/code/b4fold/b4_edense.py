"""
b4_edense (BATCH4_FOLD.md §5, optional): dense-grid recheck of EDGE_R2's ++ warning at (46, 1.045 Lam) and (50, 1.049 Lam): the ++
central values of the (EDGE-R) cost there (EDGE_R2 table: D = -3.86 / +6.57 S_E,t against the tolerance 2.05 / 2.57) come from a
degree-6 shape fit of a bulk on the 0.025 grid, which at J >= 46 is about one oscillation of D_h at x ~ 0.6 (aliasing).
Inputs:
  * EDGE_R2's merged external batch-3 n1 rows (edger2/cache/pts, read-only; even J: external structure basis);
  * LOCAL n1 nodes (b4_points.py, same code path as ext/point.py) at the 20 mid-points x = 0.3125 ... 0.7875 of the 0.025 grid (EE +
    C3H), and one local anchor node at an existing x = 0.6 (EE, state and reference) per state.
Step 'map J': the external records are mapped into the local structure basis exactly as BATCH5_FOLD §1.2 / b4_basis.py (O(2) +
reflection fit of G, H on the anchor node, E_loc = G^-T E_ext H^-1, replica by replica), then merged with the local dense nodes ->
cache/edense/pts.  Step 'gamma J': edger2/e2_gamma.main run UNCHANGED (read-only import) with its point directory redirected to
cache/edense/pts and its outputs to cache/edense (its own joint-diagonalisation repair then acts on an all-local set), and the same on
the 0.025-grid nodes only (the local anchor replaces nothing; dense mid-points excluded) as the like-for-like reference; plus a
spline-bulk variant (natural cubic through the bulk nodes instead of the K = 6 shape fit), which the dense grid makes meaningful.
Usage: python3 b4_edense.py map J | gamma J
"""
import os, sys, glob, json, pickle, shutil
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b4_guard  # noqa: F401  (first; stands in for e2_guard / er_guard / hs_guard / ng_guard)
import numpy as np
CODE = os.path.abspath(os.path.join(HERE, '..'))
E2 = os.path.join(CODE, 'edger2')
E2PTS = os.path.join(E2, 'cache', 'pts')
B4C = os.path.join(HERE, 'cache')
OUT = os.path.join(B4C, 'edense')
LOC = os.path.join(B4C, 'loc')
DREF = 20000.0


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def do_map(J):
    import b4_basis as BB              # functions only (ubasis, fit_GH, fit_general)
    D = 1000.0 + J; m = mphys(D, J)
    dst = {s: os.path.join(OUT, s) for s in ('pts', 'pts025')}
    for d in dst.values():
        os.makedirs(d, exist_ok=True)
    PL = {k: pickle.load(open(os.path.join(LOC, 'hc_pts_%s_J%d_D%g.pkl' % (k, J, D)), 'rb')) for k in ('EE', 'C3H')
          if os.path.exists(os.path.join(LOC, 'hc_pts_%s_J%d_D%g.pkl' % (k, J, D)))}
    PE = {k: pickle.load(open(os.path.join(E2PTS, 'hc_pts_%s_J%d_D%g.pkl' % (k, J, D)), 'rb')) for k in ('EE', 'C3H')}
    ext_nus = sorted(k[1] for k in PE['EE'] if k[0] == D)
    anc = [nu for nu in sorted(k[1] for k in PL['EE'] if k[0] == D) if min(abs(nu - v) for v in ext_nus) < 1e-6]
    lines = [f"# b4_edense map J={J} D={D:g}: anchor nodes (local, at external x) {[round(v/m, 4) for v in anc]}"]
    T = {}
    for DD in (D, DREF):
        Ul, Ut = BB.ubasis(DD, J), BB.ubasis(3.0 - DD, J)
        k0 = (DD, anc[0])
        rr, G, H, par = BB.fit_GH(Ul, Ut, BB.recs(PL['EE'], k0), BB.recs(PE['EE'], k0))
        rg = BB.fit_general(BB.recs(PL['EE'], k0), BB.recs(PE['EE'], k0), G, H)
        dE = float(np.abs(np.asarray(PL['EE'][k0]['E']) - np.asarray(PE['EE'][k0]['E'])).max()/np.abs(PE['EE'][k0]['E']).max())
        T[DD] = (G, H)
        lines.append(f"  Delta={DD:g}: fit resid {rr:.2e} (general 2x2 {rg:.2e}); raw |E_loc - E_ext|/|E| {dE:.2e}; params {par}")
    out = {}
    for kind in ('EE', 'C3H'):
        P = {}
        for k, r in PE[kind].items():
            G, H = T[k[0]]
            Gi, Hi = np.linalg.inv(G).T, np.linalg.inv(H)
            r = dict(r); reps = np.array([Gi @ np.asarray(q) @ Hi for q in r['reps']])
            r['reps'] = reps; r['E'] = reps.mean(0); r['err'] = reps.std(0, ddof=1)/np.sqrt(len(reps))
            r['method'] = str(r.get('method', '')) + ' [b4_edense: mapped to the local structure basis]'
            P[k] = r
        P025 = dict(P)
        dense = 0
        for k, r in PL.get(kind, {}).items():
            if min(abs(k[1] - v) for v in ext_nus) > 1e-6:
                P[k] = r; dense += k[0] == D
        out[kind] = (P, P025, dense)
        pickle.dump(P, open(os.path.join(dst['pts'], 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'wb'))
        pickle.dump(P025, open(os.path.join(dst['pts025'], 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)), 'wb'))
        lines.append(f"  {kind}: external {sum(1 for k in PE[kind] if k[0] == D)} state records mapped; local dense state records added {dense}")
    # check at the anchor: mapped vs local (both kinds where computed)
    for kind in PL:
        for k in PL[kind]:
            if k in out[kind][1]:
                a, b = np.asarray(out[kind][1][k]['E']), np.asarray(PL[kind][k]['E'])
                lines.append(f"  check {kind} {k[0]:g} x={k[1]/m:.4f}: |E_mapped - E_local|/|E| = {np.abs(a - b).max()/np.abs(b).max():.2e}")
    open(os.path.join(OUT, 'map_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    print("\n".join(lines))


def do_gamma(J):
    sys.path.insert(0, E2)
    import e2_gamma as G               # read-only import; e2_guard / er_guard resolve to b4_guard
    res = {}
    for tag, sub, KB in (('dense', 'pts', 6), ('grid025', 'pts025', 6), ('dense_spl', 'pts', -1), ('grid025_spl', 'pts025', -1)):
        od = os.path.join(OUT, tag); os.makedirs(od, exist_ok=True)
        G.PTS = os.path.join(OUT, sub); G.FIX = os.path.join(od, 'ptsfix'); G.CACHE = od
        G.HN.C = G.PTS; G.HN._SC.clear(); G.HN.RESID.clear(); G.OFFNAN.clear(); G.FAILED.clear()
        if KB < 0:
            _init = G.Model.__init__

            def _spl_init(self, J_, r, a, KB_, _init=_init):
                _init(self, J_, r, a, 6 if KB_ < 0 else KB_)
                if KB_ < 0:                     # spline bulk for the diagonal entries (variant)
                    xs = np.asarray(r['xs'], float); bs = (xs >= 0.0) & (xs < a - 1e-9)
                    for key in G.keys(J_):
                        if key[1] != 'off':
                            self.bulk[key] = G.HC.spline_of(xs[bs], np.asarray(r[G.gkey(J_, key)], float)[bs])
            G.Model.__init__ = _spl_init
            try:
                # repr variants KB-2.. of a spline are meaningless: run with KB = -1 -> Model(-3..0) also spline; report det and mc only
                r = G.main(J, -1)
            finally:
                G.Model.__init__ = _init
        else:
            r = G.main(J, KB)
        res[tag] = {h: {kd: dict(det=e[kd]['det']/r['SEt'], bulk=e[kd]['bulk']/r['SEt'], edge=e[kd]['edge']/r['SEt'],
                                 mc=e[kd]['mc']/r['SEt'], repr=e[kd]['repr']/r['SEt']) for kd in ('C2', 'C3H')}
                    for h, e in r['entries'].items()}
        res[tag]['tol_tr'] = r['tol_tr']; res[tag]['nb'] = None
    json.dump(res, open(os.path.join(OUT, 'edense_J%d.json' % J), 'w'), indent=1, default=float)
    lines = [f"# b4_edense gamma J={J} (S_E,t units; tolerance {res['dense']['tol_tr']:.3f})"]
    for tag, r in res.items():
        for h in ('++', '+-'):
            c2, c3 = r[h]['C2'], r[h]['C3H']
            lines.append(f"  {tag:12s} {h}: |D2|+|D3| = {abs(c2['det']) + abs(c3['det']):.3f}  (D2 {c2['det']:+.3f} [bulk {c2['bulk']:+.3f}] +- mc {c2['mc']:.3f} "
                         f"repr {c2['repr']:.3f}; D3 {c3['det']:+.3f} [bulk {c3['bulk']:+.3f}] +- mc {c3['mc']:.3f} repr {c3['repr']:.3f})")
    open(os.path.join(OUT, 'edense_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    {'map': do_map, 'gamma': do_gamma}[sys.argv[1]](int(sys.argv[2]))
