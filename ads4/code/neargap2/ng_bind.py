"""
ng_bind (NEARGAP2.md §1): error-bounded margin of a near-gap state of Omega^EE at Delta_g = 1000 (or 300) on its OWN-MASS K2 row.
For each requested state (J, m/Lam) it builds, in-process, through ng_chain (= the unchanged EE_NEARGAP chain, own rows added):
  * the chain model with the published rows only (NG_NOOWN)                       -> lam_model      (H_M as published)
  * the chain with the own-mass row added, Laplace order K = 2                   -> lam_K2         (central value)
  * the same at K = 1                                                             -> truncation estimate per node (the K = 2 term)
  * linear instead of spline x-interpolation; every-other-node subsample          -> interpolation estimate
and the error budget (sign-free, each with a stated safety factor), per x-window of the state:
  E_trunc  : 3 x sum_nodes |w| |Gamma_K2 - Gamma_K1|                         x in [x0, xmax] of the own row
  E_interp : 3 x max(|lam_spline - lam_linear|, |lam - lam_subsample|)
  E_low    : G_low x W(nu < nu_low) with G_low = 1e4 c2  (no data / engine only; x10 above every small-nu value in the record)
  E_edge   : G_edge x W(x > xmax) with G_edge = 10 x max |Gamma| on the row's last 4 nodes (>= 1e4 c2)
  E_base   : the hybrid base: certified flat floor N_flat/5 (CERT_CLOSE) and the Lemma-C bound on the LR curvature (POS_ERRORS §2.1)
Output: cache/ng_bind_<tag>.txt / .json.
Usage: python3 ng_bind.py Dg tag J:Delta_own:mfl[,mfl...] [J:...]
"""
import os, sys, json, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ng_chain as NC
import numpy as np
E, Q, W = NC.E, NC.Q, NC.W
NGC = NC.NGC
G_LOW = 1.0e4
NU_LOW = float(os.environ.get('NG_NULOW', '50'))


_GF = {}


def gauge_factor(Dg, nu):
    """max(1, |Re X_min/psi2|)(nu): the min-norm de-shift weight of the k = 2 part relative to the local psi2 node weight (u3/x2 caches,
    read-only); 1 above nu = 25 (W2_GAUGE convention).  Applied to the sign-free small-nu budget so that it holds in the min-norm gauge
    (the functional's own) as well as in the chain's local form.  (It multiplies the k = 3 weights too: conservative.)"""
    key = int(round(Dg))
    if key not in _GF:
        a = pickle.load(open(os.path.join(NC.SFC, 'x2_chi_EE_can_%d.pkl' % key), 'rb'))
        b = pickle.load(open(os.path.join(NC.SFC, 'u3_chi_EE_%d.pkl' % key), 'rb'))
        phi = np.interp(a['nu'], b['nu'], b['phi'].real)
        with np.errstate(divide='ignore', invalid='ignore'):
            r = np.abs(a['X_minnorm'].real/phi)
        r = np.where(np.isfinite(r), r, 1.0)
        _GF[key] = (a['nu'], r)
    g, r = _GF[key]
    # local maximum over a +-0.5 window (the ratio has poles where psi2's grid value is small)
    out = np.ones_like(nu)
    for i, v in enumerate(nu):
        if v <= 25.0:
            sel = np.abs(g - v) <= 0.5
            out[i] = max(1.0, float(r[sel].max()) if sel.any() else 1.0)
    return out


def grid_for(J, K, noown, Ds=None):
    NC.NOOWN = noown
    if J % 2:
        NC.KK, NC.PK = K, None
    else:
        NC.KK, NC.PK = 2, K
    W.DESPIKE_LOG.clear()
    g = Q.odd_grid(J) if J % 2 else Q.even_grid(J)
    if Ds is not None:
        g = {D: r for D, r in g.items() if D in Ds}
    return g, list(W.DESPIKE_LOG)


def model_of(J, g):
    grids = {'even': {}, 'odd': {}}
    grids['odd' if J % 2 else 'even'][J] = g
    em, om = NC.models(grids)
    return (om if J % 2 else em)[J][0]


def lin_model(J, g):
    """x-linear interpolation variant (natural spline -> np.interp, held outside), same 1/m interpolation."""
    from scipy.interpolate import CubicSpline as CS
    import l9_neargap as NGm
    class _Lin:
        def __init__(self, x, y, bc_type=None):
            self.x, self.y = np.asarray(x), np.asarray(y)
        def __call__(self, t):
            return np.interp(np.asarray(t, float), self.x, self.y)
    Q.NG.CubicSpline = _Lin; Q.CubicSpline = _Lin
    try:
        return model_of(J, g)
    finally:
        Q.NG.CubicSpline = CS; Q.CubicSpline = CS


def sub(g, keep):
    out = {}
    for D, r in g.items():
        n = len(r['xs']); sel = np.zeros(n, bool); sel[keep::2] = True; sel[0] = True; sel[-1] = True
        out[D] = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == (n,) else v) for k, v in r.items()}
    return out


def errmat(d, sel, G):
    """sign-free bound on the change of the normalised matrix if |Gamma| <= G (array over nodes or scalar) on nodes sel."""
    hs = d['hs']; n = len(hs); Ns = d['Ns']
    Gv = np.broadcast_to(G, d['nu'].shape)
    M = np.zeros((n, n))
    for i, h in enumerate(hs):
        a2, a3 = d['abs'][h]
        M[i, i] = ((a2 + a3)*Gv)[sel].sum()/Ns[i]
    if n == 2:
        a2, a3 = d['abs']['off']
        M[0, 1] = M[1, 0] = ((a2 + a3)*Gv)[sel].sum()/np.sqrt(Ns[0]*Ns[1])
    return np.sqrt((M*M).sum())


def absdens(Dg, J, mfl, npn):
    """|weight| per unit |Gamma| (c2) per node, per helicity and off-diagonal (as POS_ERRORS large_m_map, odd J: +- only)."""
    Lm = E.lam_of(Dg)
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    st = {h: E.hyb(Dg, J, mfl, h, npn) for h in hs}
    s0 = st[hs[0]]; x, m, u, y = s0['x'], s0['m'], s0['u'], s0['y']
    env = Q.P4.env if J % 2 else Q.A.env
    out = {h: (np.abs(st[h]['p2'])*env(h, x, 2)/m**2, u*np.abs(st[h]['p3'])*env(h, x, 3)/m**2) for h in hs}
    if len(hs) == 2:
        eo2 = np.sqrt(env('++', x, 2)*env('+-', x, 2))
        out['off'] = (np.abs(s0['p2'])*eo2/m**2, u*np.abs(s0['p3'])*eo2/(2 - y)/m**2)
    return out


def interp_nodes(g, J, m, x, key_fn):
    """|Gamma_K2 - Gamma_K1| at the state's nodes x (linear in x between row nodes, held outside), 1/m-interpolated as the model."""
    Ds = sorted(g); ms = np.array([g[D]['m'] for D in Ds])
    vals = [np.interp(np.minimum(x, g[D]['xs'][-1]), g[D]['xs'], key_fn(g[D])) for D in Ds]
    if m <= ms[0]:
        return vals[0]
    if m >= ms[-1]:
        return vals[-1]
    i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
    t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
    return (1 - t)*vals[i0] + t*vals[i1]


def run_state(Dg, J, Down, mfl, pr, npn=256):
    Lm = E.lam_of(Dg); m = mfl*Lm
    rows_own = {float(Down)}
    out = dict(J=J, mfl=mfl, Down=Down)
    NC.PK = None
    g2, dsp2 = grid_for(J, 2, False)
    g1, dsp1 = grid_for(J, 1, False)
    g0m, _ = grid_for(J, 2, True)
    f2 = model_of(J, g2); f0m = model_of(J, g0m); f1 = model_of(J, g1)
    d2 = NC.dens(Dg, J, mfl, f2, npn); dm = NC.dens(Dg, J, mfl, f0m, npn)
    lam2, dF2 = NC.lam_of(d2); lamm, dFm = NC.lam_of(dm)
    lam0 = E.lam2(d2['base'], d2['Ns'], 0*dF2)
    lamK1 = NC.lam_of(NC.dens(Dg, J, mfl, f1, npn))[0]
    # chain's own 64-node value for the same model (consistency with q3's quadrature)
    d64 = NC.dens(Dg, J, mfl, f2, 64); lam64 = NC.lam_of(d64)[0]
    nu, x = d2['nu'], d2['x']
    own = g2.get(float(Down))
    xs = own['xs'] if own is not None else np.array([0.05, 0.8])
    x0, xmax = xs[0], xs[-1]
    hs = d2['hs']
    d2['abs'] = absdens(Dg, J, mfl, npn)
    # truncation: per node |Gamma_K2 - Gamma_K1| for each (kind, entry), common nodes
    tr = {}
    for kd in ('C2', 'C3H'):
        for h in (hs + (['off'] if len(hs) == 2 else [])):
            key = kd if J % 2 else (kd, h)
            def kf(r, key=key):
                return r[key]
            gd = {}
            for D in g2:
                if D not in g1:
                    continue
                x2r, x1r = np.round(g2[D]['xs'], 6), np.round(g1[D]['xs'], 6)
                xc = np.intersect1d(x2r, x1r)           # nodes kept at both orders (the Hampel filter may differ per order)
                if len(xc) < 2:
                    continue
                v2 = np.array([g2[D][key][int(np.argmin(np.abs(x2r - v)))] for v in xc])
                v1 = np.array([g1[D][key][int(np.argmin(np.abs(x1r - v)))] for v in xc])
                gd[D] = dict(m=g2[D]['m'], xs=xc, v=np.abs(v2 - v1))
            if isinstance(key, tuple) and key[1] == 'off':
                gd = {D: dict(v, v=v['v']*g2[D]['m']) for D, v in gd.items()}       # c_off = m Gamma_off (OFFSCALE), weight 1/m^2
            tr[(kd, h)] = interp_nodes(gd, J, m, x, lambda r: r['v'])
    # truncation error matrix (sign-free with the node-wise |dGamma|, x 3)
    n = len(hs); Mt = np.zeros((n, n))
    win = (x >= x0) & (x <= xmax)
    for i, h in enumerate(hs):
        a2, a3 = d2['abs'][h]
        Mt[i, i] = (a2*tr[('C2', h)] + a3*tr[('C3H', h)])[win].sum()/d2['Ns'][i]
    if n == 2:
        a2, a3 = d2['abs']['off']
        Mt[0, 1] = Mt[1, 0] = (a2*tr[('C2', 'off')] + a3*tr[('C3H', 'off')])[win].sum()/np.sqrt(d2['Ns'].prod())
    E_trunc = 3*np.sqrt((Mt*Mt).sum())
    # interpolation: linear vs spline, every-other-node subsample (own row only; other rows unchanged)
    fl = lin_model(J, g2); laml = NC.lam_of(NC.dens(Dg, J, mfl, fl, npn))[0]
    lsub = []
    for keep in (0, 1):
        gs = dict(g2)
        if own is not None:
            gs[float(Down)] = sub({float(Down): own}, keep)[float(Down)]
        lsub.append(NC.lam_of(NC.dens(Dg, J, mfl, model_of(J, gs), npn))[0])
    E_interp = 3*max(abs(lam2 - laml), *[abs(lam2 - v) for v in lsub])
    # small-nu window: nu < max(NU_LOW, x0 m) -- sign-free G_LOW; contribution of the model there is also removed from the central value
    lowsel = nu < max(NU_LOW, x0*m)
    E_low = errmat(d2, lowsel, G_LOW*gauge_factor(Dg, nu))
    # edge window x > xmax: sign-free 10 x the largest |Gamma| of the last four own-row nodes (at least 1e4 c2)
    gedge = 1e4
    if own is not None:
        for kd in ('C2', 'C3H'):
            for h in hs:
                key = kd if J % 2 else (kd, h)
                gedge = max(gedge, 10*np.max(np.abs(np.asarray(own[key])[-4:])))
    edsel = x > xmax
    E_edge = errmat(d2, edsel, gedge)
    # central value with the windows' model contributions removed (they are fully inside the sign-free budgets)
    lamc = NC.lam_of(d2, sel=~(lowsel | edsel))[0]
    # base: certified flat floor and Lemma C
    base_info = {}
    for i, h in enumerate(hs):
        st = E.hyb(Dg, J, mfl, h, npn)
        beta = 2*np.arcsinh((J + .5)/m)
        dP = (1 - 1/np.cosh(beta/2)) + np.minimum(2.0, st['pn']*Lm*(2*np.sinh(beta/2) - beta))
        cvb = (np.abs(st['p2'])*2*dP).sum() + (np.abs(st['p3'])*st['u']*dP).sum()
        Nf = abs(st['B2']) + abs(st['B3'])
        base_info[h] = dict(lam_hyb=float(st['base']/st['Ns']), flat=float((st['B2'] + st['B3'])/Nf), Nf_Ns=float(Nf/st['Ns']),
                            cv_Ns=float((st['cv2'] + st['cv3'])/st['Ns']), cvb_Ns=float(cvb/st['Ns']),
                            rig=float((Nf/5 - cvb)/st['Ns']))
    tot = E_trunc + E_interp + E_low + E_edge
    E_tl = 3*abs(lam2 - lamK1)                     # truncation at the lambda level: 3 x the signed effect of the last order kept
    tot_l = E_tl + E_interp + E_low + E_edge
    out.update(lam_hybrid=float(lam0), lam_model=float(lamm), lam_K2=float(lam2), lam_K2_64=float(lam64), lam_central=float(lamc),
               lam_lin=float(laml), lam_sub=[float(v) for v in lsub], E_trunc=float(E_trunc), E_interp=float(E_interp),
               E_low=float(E_low), E_edge=float(E_edge), G_edge=float(gedge), E_tot=float(tot), lower=float(lamc - tot),
               L=float(lamc/tot) if tot > 0 else np.inf, lam_K1=float(lamK1), E_trunc_lam=float(E_tl), E_tot_lam=float(tot_l),
               lower_lam=float(lamc - tot_l), L_lam=float(lamc/tot_l) if tot_l > 0 else np.inf, base=base_info, rows=sorted(float(D) for D in g2), rows_model=sorted(float(D) for D in g0m),
               despike=dsp2, dF_model=dFm.tolist(), dF_K2=dF2.tolist(), x0=float(x0), xmax=float(xmax))
    pr(f"  J={J:2d} m/Lam={mfl:.5f} rows {out['rows']} (model {out['rows_model']}): hybrid {lam0:+.5f} | model chain {lamm:+.5f} | "
       f"own row K2 {lam2:+.5f} (64 nodes {lam64:+.5f}; lin {laml:+.5f}; sub {lsub[0]:+.5f}/{lsub[1]:+.5f}) | central (windows cut) {lamc:+.5f}")
    pr(f"        errors: trunc {E_trunc:.2e}  interp {E_interp:.2e}  low(nu<{max(NU_LOW, x0*m):.0f}, G={G_LOW:.0e}) {E_low:.2e}  "
       f"edge(x>{xmax:.3f}, G={gedge:.2e}) {E_edge:.2e}  -> total {tot:.2e}; lower bound {lamc - tot:+.5f}; L = {out['L']:.1f}")
    pr(f"        lambda-level truncation: K1 chain {lamK1:+.5f} -> 3|K2 - K1| = {E_tl:.2e}; total {tot_l:.2e}; lower bound {lamc - tot_l:+.5f}; L = {out['L_lam']:.1f}")
    for h, b in base_info.items():
        pr(f"        base {h}: hybrid {b['lam_hyb']:+.5f}, flat {b['flat']:+.5f} (Nf/Ns {b['Nf_Ns']:.4f}), LR curvature {b['cv_Ns']:+.2e} N_s, "
           f"Lemma-C bound {b['cvb_Ns']:.2e} N_s, certified-floor route (Nf/5 - LemmaC)/Ns = {b['rig']:+.5f}")
    if dsp2:
        pr(f"        despiked: {dsp2}")
    return out


if __name__ == '__main__':
    Dg = float(sys.argv[1]); tag = sys.argv[2]
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    t0 = time.time()
    res = []
    pr(f"# ng_bind Dg={Dg:g}: per state: hybrid base | published-row chain | own-row chain (K2) | error budget (x3 safety on trunc/interp)")
    for spec in sys.argv[3:]:
        Js, Dw, mf = spec.split(':')
        for v in mf.split(','):
            res.append(run_state(Dg, int(Js), float(Dw), float(v), pr))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(NGC, 'ng_bind_%s.txt' % tag), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(NGC, 'ng_bind_%s.json' % tag), 'w'))
