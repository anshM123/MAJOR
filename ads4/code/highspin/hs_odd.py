"""
hs_odd (HIGHSPIN.md §1-§3): an odd-J (parity-odd, 1x1 +-) near-gap state of Omega^EE at Delta_g = 1000 on its own-mass row
Delta = 1000 + J, with DIRECT n1 nodes (hs_points.py) for the light-cone edge x in [x_cut, Lam/m] and, where they exist, for the bulk.

Gamma(x): the chain's rows (kappa-curvature image subtracted, SMOOTH_BOUND §4 = the published odd-J convention; the base carries the
leading-Regge curvature once).  K2 nodes at which the Delta_ref = 2e4 reference's Laplace series fails (hs_core.refbad: the reference's
leading partial wave passes through its first zero there) are dropped.
Bulk x in [0.05, x_cut): shape fit (hs_core.ShapeFit, degree KB) of the K2 nodes ('K2' bulk) or of the direct nodes ('dir' bulk, used
as the main value when >= 90 % of the bulk nodes are direct).  Below x = 0.05 the fit's polynomial continuation (as published); that
window (nu < 50) is ALSO covered sign-free (E_low).
Edge x in [x_cut, Lam/m]:
  gk    : Gauss-Kronrod-15 quadrature of the DIRECT values at the GK15 nodes of [x_cut, Lam/m_row] (twist-edge state only; no
          interpolation); GK nodes without a valid direct value (RQMC failure as nu -> Delta) are set to 0 and covered sign-free with
          G = max(1e4, 10 max|Gamma_edge|) (E_miss); 'g7' = the Gauss-7 subset (quadrature error estimate);
  espl  : natural spline through all valid direct edge nodes, held beyond the last one (any mass);
  eshape: shape fit (degree KE) of the direct edge nodes (if enough nodes);
  held  : the bulk fit continued (the published treatment beyond its last node).
Published configuration ('pub'): K2 nodes 0.05 <= x <= 0.8 minus the reference-zero node, shape K = 6 over everything (reproduces
sm_chainfix fodd_hampel_K6_row_kappa: J = 29 twist edge +0.9477 vs +0.948).
Budget, lambda level (x3 on every model term, 3 sigma on the Monte Carlo):
  E_mc    : 3 x jackknife sigma (8 RQMC replicas) of the main value (all direct nodes, bulk and edge),
  E_quad  : 3 x |gk - g7| (twist edge) or 3 x max|espl - eshape, held| (other masses),
  E_miss  : sign-free cost of the missing GK nodes with |Gamma| <= 3 A ((1 - x_l)/(1 - x))^2 (A = max|Gamma| of the last three valid
            nodes, x_l the last one: the measured growth of the edge amplitude ~ (1 - x)^-2, with x3), plus, for a GK node within ZWIN of the
            reference zero, the sign-free cost of 3 x the largest |Gamma| of its 4 nearest valid edge nodes,
  E_trunc : 3 x |lam(K2 bulk, K = 2) - lam(K = 1)|  (only when the bulk is K2),
  E_repr  : 3 x max |lam(bulk degree KB) - lam(KB - 2, KB - 1, KB + 1; two half-samples at degree 4)|,
  E_low   : sign-free |Gamma| <= G on nu < 50: G = 1e4 (NEARGAP2 convention) and G = 10 x max|Gamma| of the three lowest nodes.
Usage: python3 hs_odd.py J [mfl,...]   -> cache/hs_odd_J<J>.txt / .json
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_core as H
import numpy as np
from numpy.polynomial.legendre import leggauss

Dg = 1000
XCUT = float(os.environ.get('HS_XCUT', '0.8'))
KB = int(os.environ.get('HS_KB', '6'))
KE = int(os.environ.get('HS_KE', '3'))
RTOL = 0.1
G_LOW = 1e4


def rows(J, D, K=2, jk=None):
    rK = H.subtract_kappa(H.grid(J, K), J)[D]
    rI = H.subtract_kappa(H.grid(J, K, inject=[D], jk=jk), J)[D]
    xd = set(np.round(H.direct_x(J, D), 6)) - H.failed_x(J, D)
    isdir = np.array([round(x, 6) in xd for x in rI['xs']])
    return rK, rI, isdir


class Model:
    """bulk fit + edge data of one row (both kinds)."""
    def __init__(self, J, rK, rI, isdir, bad, bulk='K2', kb=KB, sub=None, xcut=XCUT):
        self.J = J
        if bulk == 'dir':
            xs, src = rI['xs'], rI
            bsel = (xs >= 0.05 - 1e-9) & (xs < xcut - 1e-9) & isdir
        else:
            xs, src = rK['xs'], rK
            bsel = (xs >= 0.05 - 1e-9) & (xs < xcut - 1e-9) & np.array([round(x, 6) not in bad for x in xs])
        if sub is not None:
            idx = np.where(bsel)[0]; keep = np.zeros_like(bsel); keep[idx[sub::2]] = True; keep[idx[0]] = keep[idx[-1]] = True
            bsel = keep
        self.nb = int(bsel.sum()); self.bx = xs[bsel]
        xe = rI['xs']; esel = (xe >= xcut - 1e-9) & isdir
        self.ex = xe[esel]
        self.k = {}
        for kd in ('C2', 'C3H'):
            ey = np.asarray(rI[kd], float)[esel]
            self.k[kd] = dict(bulk=H.ShapeFit(J, xs[bsel], np.asarray(src[kd], float)[bsel], '+-', K=kb), ey=ey,
                              eshape=H.ShapeFit(J, self.ex, ey, '+-', K=KE, x0=xcut) if len(ey) > 3*(KE + 1) + 2 else None,
                              espl=H.spline_of(self.ex, ey) if len(ey) >= 3 else None)
        self.gmax = max(1e4, 10*max([np.abs(self.k[kd]['ey']).max() for kd in self.k if len(self.k[kd]['ey'])] + [0]))


ZWIN = 0.002      # |x - x_refzero| < ZWIN: direct values not used (MC noise x 1/|D_h(x_ref)|)
ZCOST = [0.0]
XROW = [None]       # Lam/m_row (upper end of the row's GK window)


def value(J, mfl, M, edge, xcut=XCUT, nb=800):
    """lam/N_s of the state (J, mfl) with bulk model M and edge treatment; returns (lam, base/N_s, sign-free missing-node cost/N_s)."""
    ZCOST[0] = 0.0
    hs, base, Ns, st = H.base_of(Dg, J, mfl)
    xmax = 1.0/mfl
    t, w = leggauss(nb); xb = 0.5*xcut*(t + 1); wb = 0.5*xcut*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xb)
    d = (wb*w2*M.k['C2']['bulk'](xb)).sum() + (wb*w3*M.k['C3H']['bulk'](xb)).sum()
    miss = 0.0
    if edge in ('gk', 'g7'):
        xe, wk, wg = H.gk15_nodes(xcut, XROW[0])    # the ROW's GK nodes (= the state's own window at the twist edge)
        cut = (xe <= xmax).astype(float)             # p = x m/Lam <= 1; psi_k ~ (1-p)^5: the integrand is C^4 at the cut
        ww = (wk if edge == 'gk' else wg)*cut
        w2e, w3e = H.state_weights(Dg, J, mfl, '+-', xe)
        have = np.array([np.any(np.abs(M.ex - x) < 5e-6) for x in xe])
        x0 = H.ref_zero_x(J, 1000.0 + J)
        inz = np.abs(xe - x0) < ZWIN                    # reference-zero window: the direct ratio is unusable there too
        for kd, we in (('C2', w2e), ('C3H', w3e)):
            ey = np.array([M.k[kd]['ey'][np.argmin(np.abs(M.ex - x))] if h else 0.0 for x, h in zip(xe, have)])
            ey[inz] = 0.0
            d += (ww*we*ey).sum()
            # missing (RQMC-failed, x -> Lam/m) nodes: 3 x the light-cone envelope A ((1 - x_l)/(1 - x))^2 fitted to the last valid nodes
            okx = np.abs(M.ex - x0) >= ZWIN
            last = np.argsort(M.ex + -9*(~okx))[-3:]
            A = max(np.abs(M.k[k_]['ey'][last]).max() for k_ in ('C2', 'C3H')); xl = M.ex[last].max()
            genv = 3*A*np.maximum(1.0, ((1 - xl)/(1 - xe))**2)
            miss += (genv*ww*np.abs(we))[~have & ~inz].sum()
            if inz.any():                               # sign-free: 3 x the largest |Gamma| of the 4 nearest valid edge nodes
                okx = np.abs(M.ex - x0) >= ZWIN
                near = np.argsort(np.abs(M.ex - x0) + 1e9*(~okx))[:4]
                gl = 3*max(np.abs(M.k[k_]['ey'][near]).max() for k_ in ('C2', 'C3H'))
                ZCOST[0] += gl*(ww*np.abs(we))[inz].sum()
    else:
        t, w = leggauss(300); xx = xcut + 0.5*(xmax - xcut)*(t + 1); wx = 0.5*(xmax - xcut)*w
        w2e, w3e = H.state_weights(Dg, J, mfl, '+-', xx)
        for kd, we in (('C2', w2e), ('C3H', w3e)):
            f = M.k[kd]['bulk'] if edge == 'held' else M.k[kd][edge]
            if f is None:
                return np.nan, np.nan, np.nan
            d += (wx*we*f(xx)).sum()
    ZCOST[0] /= Ns[0]
    return float((base[0] + d)/Ns[0]), float(base[0]/Ns[0]), float(miss/Ns[0])


def pub_value(J, mfl, rK, bad):
    xs = rK['xs']
    sel = (xs >= 0.05 - 1e-9) & (xs <= 0.8 + 1e-9) & np.array([round(x, 6) not in bad for x in xs])
    f = {kd: H.ShapeFit(J, xs[sel], np.asarray(rK[kd])[sel], '+-', K=6) for kd in ('C2', 'C3H')}
    hs, base, Ns, st = H.base_of(Dg, J, mfl); s = st['+-']; x = s['x']; m = s['m']
    d = (s['p2']*H.Q.P4.env('+-', x, 2)*f['C2'](x)/m**2).sum() + (s['p3']*s['u']*H.Q.P4.env('+-', x, 3)*f['C3H'](x)/m**2).sum()
    return float((base[0] + d)/Ns[0])


def wwin(J, mfl, a, b):
    t, w = leggauss(400); xx = a + 0.5*(b - a)*(t + 1); ww = 0.5*(b - a)*w
    w2, w3 = H.state_weights(Dg, J, mfl, '+-', xx)
    _, _, Ns, _ = H.base_of(Dg, J, mfl)
    return float((ww*(np.abs(w2) + np.abs(w3))).sum()/Ns[0])


def jk_sigma(v):
    v = np.asarray(v); return float(np.sqrt(7/8*((v - v.mean())**2).sum()))


if __name__ == '__main__':
    J = int(sys.argv[1]); D = 1000.0 + J
    rowm = H.HN.mphys(D, J); Lm = H.E.lam_of(Dg); mtw = rowm/Lm
    XROW[0] = 1.0/mtw
    mfls = [mtw] + ([float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else [1.04, 1.05, 1.06])
    lines, res = [], dict(J=J, D=D, states={})

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    rK, rI, isdir = rows(J, D)
    bad = H.refbad(J, D, RTOL)
    xd = H.direct_x(J, D); fx = H.failed_x(J, D)
    pr(f"# hs_odd J={J} row {D:g} (m={rowm:.3f}, J^2/m={J*J/rowm:.3f}), twist edge m/Lam={mtw:.7f}, x_max={1/mtw:.6f}; x_cut={XCUT}; "
       f"bulk shape K={KB}, edge shape K={KE}")
    pr(f"#  direct nodes {len(xd)} (failed, RQMC rel. err > 0.2: {sorted(fx)})")
    pr(f"#  K2 nodes dropped (reference series fails, excess |r2/r0| > {RTOL}): {sorted(bad)}; reference zero predicted at "
       f"x = {H.ref_zero_x(J, D):.4f}")
    reps = [rows(J, D, jk=j)[1] for j in range(8)]
    pr("#  node table (Gamma - kappa, c2): x | C2: K2, direct +- jk | C3H: K2, direct +- jk | K2 Laplace step |K2-K1|")
    rK1 = H.subtract_kappa(H.grid(J, 1), J)[D]
    nodetab = []
    for i, x in enumerate(rI['xs']):
        if not isdir[i]:
            continue
        s = f"   x={x:.4f}"; rec = dict(x=float(x))
        for kd in ('C2', 'C3H'):
            vj = np.array([rr[kd][int(np.argmin(np.abs(rr['xs'] - x)))] for rr in reps]); e = jk_sigma(vj)
            j0 = int(np.argmin(np.abs(rK['xs'] - x))); on = abs(rK['xs'][j0] - x) < 1e-6
            k2 = f"{rK[kd][j0]:+9.1f}" if on else "     --  "
            st = f"{abs(rK[kd][j0] - rK1[kd][int(np.argmin(np.abs(rK1['xs'] - x)))]):7.0f}" if on else "     --"
            s += f" | {kd}: K2 {k2} dir {rI[kd][i]:+9.1f} +- {e:6.1f} step {st}"
            rec[kd] = dict(K2=float(rK[kd][j0]) if on else None, dir=float(rI[kd][i]), err=e)
        nodetab.append(rec); pr(s)
    res['nodes'] = nodetab
    nbk = int(((rK['xs'] >= 0.05 - 1e-9) & (rK['xs'] < XCUT - 1e-9)).sum())
    nbd = int(((rI['xs'] >= 0.05 - 1e-9) & (rI['xs'] < XCUT - 1e-9) & isdir).sum())
    bulk_main = 'dir' if nbd >= 0.9*nbk else 'K2'
    M = {b: Model(J, rK, rI, isdir, bad, bulk=b) for b in (['K2', 'dir'] if nbd >= 3*(KB + 1) + 2 else ['K2'])}
    MK1 = Model(J, rK1, rI, isdir, bad, bulk='K2')
    dset = set(np.round(xd, 6)) - fx
    Mj = [Model(J, rK, rr, np.array([round(x, 6) in dset for x in rr['xs']]), bad, bulk=bulk_main) for rr in reps]
    Mrep = {f'K{k}': Model(J, rK, rI, isdir, bad, bulk=bulk_main, kb=k) for k in (KB - 2, KB - 1, KB + 1)}
    Mrep.update({f'half{s}': Model(J, rK, rI, isdir, bad, bulk=bulk_main, kb=4, sub=s) for s in (0, 1)})
    pr(f"#  bulk: {bulk_main} ({M[bulk_main].nb} nodes; K2 nodes {nbk}, direct {nbd}); edge direct nodes {len(M[bulk_main].ex)}; fit rms/rms "
       + "; ".join(f"{b}: C2 {m_.k['C2']['bulk'].rms_rel:.3f} C3H {m_.k['C3H']['bulk'].rms_rel:.3f}" for b, m_ in M.items()))
    xe, _, _ = H.gk15_nodes(XCUT, 1/mtw)
    for mfl in mfls:
        tw = abs(mfl - mtw) < 1e-6
        Mm = M[bulk_main]
        have = np.array([np.any(np.abs(Mm.ex - x) < 5e-6) for x in xe])
        o = dict(hyb=value(J, mfl, Mm, 'held')[1], pub=pub_value(J, mfl, rK, bad))
        gk_ok = have.sum() >= 10
        for edge in (['gk', 'g7'] if gk_ok else []) + ['espl', 'eshape', 'held']:
            for b, m_ in M.items():
                lam, _, miss = value(J, mfl, m_, edge)
                o[f'{b}/{edge}'] = lam
                if edge == 'gk' and b == bulk_main:
                    o['E_miss'] = miss; o['E_zero'] = ZCOST[0]
        edge_main = 'gk' if gk_ok else ('espl' if np.isfinite(o[f'{bulk_main}/espl']) else 'held')
        main = f'{bulk_main}/{edge_main}'
        o['main'] = o[main]; o['main_key'] = main
        lj = [value(J, mfl, mj, edge_main)[0] for mj in Mj]; sj = jk_sigma(lj)
        o['jk'] = sj
        E_mc = 3*sj
        if gk_ok:
            E_quad = 3*abs(o[f'{bulk_main}/gk'] - o[f'{bulk_main}/g7'])
        else:
            alt = [o[f'{bulk_main}/{e}'] for e in ('eshape', 'held') if np.isfinite(o[f'{bulk_main}/{e}']) and e != edge_main]
            E_quad = 3*max([abs(o[main] - a) for a in alt] + [0.0]) if edge_main != 'held' else np.nan
        E_miss = o.get('E_miss', 0.0) + o.get('E_zero', 0.0)
        o['K1'] = value(J, mfl, MK1, edge_main)[0]
        E_tr = 3*abs(o['K2/' + edge_main] - o['K1']) if bulk_main == 'K2' else 0.0
        o['repr'] = {k: value(J, mfl, m_, edge_main)[0] for k, m_ in Mrep.items()}
        E_rep = 3*max(abs(o['main'] - v) for v in o['repr'].values())
        Wl = wwin(J, mfl, 0, 50.0/(mfl*Lm))
        lowx = np.sort(rI['xs'][(rI['xs'] >= 0.05 - 1e-9)])[:3]
        glow_meas = 10*max(abs(rI[kd][int(np.argmin(np.abs(rI['xs'] - x)))]) for kd in ('C2', 'C3H') for x in lowx)
        tot = E_mc + E_quad + E_miss + E_tr + E_rep
        o.update(E_mc=E_mc, E_quad=E_quad, E_miss=E_miss, E_trunc=E_tr, E_repr=E_rep, E_tot=tot, lower=o['main'] - tot,
                 L=o['main']/tot if tot > 0 else np.inf, W_low=Wl, E_low=G_LOW*Wl, G_low_meas=glow_meas, E_low_meas=glow_meas*Wl,
                 W_bulk=wwin(J, mfl, 0.05, XCUT), W_edge=wwin(J, mfl, XCUT, 1/mfl))
        o['L_low'] = o['main']/(tot + o['E_low']); o['L_lowmeas'] = o['main']/(tot + o['E_low_meas'])
        o['G_tol_low'] = (o['main'] - tot)/Wl if Wl > 0 else np.inf
        res['states'][f'{mfl:.7f}'] = o
        pr(f"  m/Lam={mfl:.7f}{' (twist edge)' if tw else ''}: hybrid {o['hyb']:+.4f} | published config {o['pub']:+.4f} | " +
           " | ".join(f"{k} {v:+.4f}" for k, v in o.items() if '/' in k and np.isfinite(v)))
        pr(f"        MAIN {main} = {o['main']:+.4f} +- {sj:.4f} (jk) | bulk K1 {o['K1']:+.4f} | repr " +
           ", ".join(f"{k} {v:+.4f}" for k, v in o['repr'].items()))
        pr(f"        budget (lambda level): mc {E_mc:.3f}  quad {E_quad:.3f}  miss {E_miss:.3f}  trunc {E_tr:.3f}  repr {E_rep:.3f} -> "
           f"{tot:.3f}; lower {o['lower']:+.4f}, L = {o['L']:.2f}")
        pr(f"        nu < 50 sign-free (W = {Wl:.2e}/c2): G = 1e4 -> {o['E_low']:.3f} (L {o['L_low']:.2f}); G = 10 x max|Gamma| of the "
           f"3 lowest nodes = {glow_meas:.0f} -> {o['E_low_meas']:.3f} (L {o['L_lowmeas']:.2f}); tolerated uniform |Gamma| there: "
           f"{o['G_tol_low']:.2e} c2 | W bulk {o['W_bulk']:.2e}, edge {o['W_edge']:.2e} per c2")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(H.HSC, 'hs_odd_J%d.txt' % J), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(H.HSC, 'hs_odd_J%d.json' % J), 'w'), indent=1, default=float)
