"""
gs2_bind (GAP_SEAM2.md; copy of gapseam/gsm_bind.py: guard gs2_guard, rows also read from gapseam2/cache/rows, outputs
gapseam2/cache/gs2_bind_<TAG>.*, and GS2_XMAX (drop 2000-row nodes above XMAX, like GSM_XMIN below XMIN)).  Original docstring:
gsm_bind (GAP_SEAM.md): error-budgeted near-gap states of Omega^EE at Delta_g = 2000 on own-mass K2 rows with a PER-ROW heavy reference
Delta_ref (GAP_NEARGAP N1), plus the even-J parity-odd OO block through Gamma_OO = Gamma_++ (PARITY_ODD_CLOSE §1; checked at 2000 by
gsm_ooid.py).  Copy of gapnear/gn_bind.py (NEARGAP2's budget, neargap2/ng_bind helpers unchanged, in-process) with:
  * the row grids rebuilt by re-typed copies of p4_odd_analysis.grid / l9_analysis.grid in which the reference is looked up per row:
    Delta_ref(J, D) = REFMAP[(J, D)] if set, else 2e4 (the published value).  Every other step (loaders, basis, Hampel, de-shift,
    a-fit, Laplace pred) is the module's own function, called through the module (so all patches of w2_chain / q3 / ee1000_common /
    ng_chain apply exactly as before).  With REFMAP empty the copies reproduce gn_bind to every digit (checked: gsm_bind.py check).
  * the row loader also reads gapseam/cache/rows (and, as gn_bind, gapnear/cache/rows, gapscale/cache/rows);
  * HIGHSPIN's reference-zero rule (refbad) and ref_zero_x evaluated with the row's own reference;
  * run_state = ng_bind.run_state re-typed, plus the ++-diagonal (OO) 1x1 value and budget for even J >= 4:
    lambda_OO/N_s = (base_++ + dF_++)/N_s,++ (HYPOTHESES_CLOSE §3.1, with dF_OO = dF_++), each budget term's ++ entry, and the identity
    residual E_id = errmat(|Gamma_OO - Gamma_++| <= G_ID) on every node (G_ID from gsm_ooid, x10).
Usage: python3 gs2_bind.py TAG SPEC [SPEC ...]      SPEC = Dg:J:Down:mfl[,mfl...][:Dref]   (Dref applies to the row Down; 'tw', 'tw1000')
       python3 gsm_bind.py --hmono TAG J:D1/dref1,D2/dref2,...     (Gamma(x) tables, K2 [K2 - K1], per-row reference)
Output: cache/gs2_bind_<TAG>.txt / .json
"""
import os, sys, glob, pickle, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gs2_guard  # noqa: F401   (registers itself as gsm_guard, ng_guard, gn_guard, ...)
GSMR = os.path.join(HERE, '..', 'gapseam', 'cache', 'rows')      # read-only
GS2R = os.path.join(HERE, 'cache', 'rows')
GNR = os.path.join(HERE, '..', 'gapnear', 'cache', 'rows')
GSR = os.path.join(HERE, '..', 'gapscale', 'cache', 'rows')
GSW = os.path.join(HERE, '..', 'gapscale', 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'neargap2'))
import numpy as np
_WTS = {Dg: dict(np.load(os.path.join(GSW, 'weights_%d.npz' % Dg))) for Dg in (2000, 5000)}
import ng_chain as NC
import ng_bind as NB
W, Q, E = NC.E.W, NC.Q, NC.E
P4, A = Q.P4, Q.A
DREF0 = 20000.0
REFMAP = {}                         # (J, D) -> Delta_ref of that row
HIDE = [False]                      # True: only the rows NEARGAP2 / HIGHSPIN saw (Delta_g = 1000 reproduction)
REFBAD = os.environ.get('GN_REFBAD', '1') == '1'
REFBAD_TOL = 0.1
XMIN = float(os.environ.get('GSM_XMIN', '0'))       # diagnostic only (default 0: all nodes)
XMAX = float(os.environ.get('GS2_XMAX', '9'))       # GAP_SEAM2: drop 2000-row nodes above XMAX (default: none)
G_ID = float(os.environ.get('GSM_GID', '10'))     # c2; sign-free bound on |Gamma_OO - Gamma_++| (x10 the max measured, gsm_ooid)


def dref(J, D):
    return float(REFMAP.get((int(J), float(D)), DREF0))


_ngpk = W._pk


def _pk(pats):
    R = _ngpk(pats)
    if not NC.NOOWN and not HIDE[0]:
        for p in pats:
            for d in (GSR, GNR, GSMR, GS2R):
                for fn in sorted(glob.glob(os.path.join(d, p))):
                    R.update(pickle.load(open(fn, 'rb')))
    return R


W._pk = _pk
_gf_orig = NB.gauge_factor


def gauge_factor(Dg, nu):
    key = int(round(Dg))
    if key not in _WTS:
        return _gf_orig(Dg, nu)
    w = _WTS[key]; g = w['nu']; X = w['X'].real; phi = w['phi']
    with np.errstate(divide='ignore', invalid='ignore'):
        r = np.abs(X/phi)
    r = np.where(np.isfinite(r), r, 1.0)
    out = np.ones_like(nu)
    for i, v in enumerate(nu):
        if v <= 25.0:
            sel = np.abs(g - v) <= 0.5
            out[i] = max(1.0, float(r[sel].max()) if sel.any() else 1.0)
    return out


NB.gauge_factor = gauge_factor


# ------------------------------------------------------------------------------------------------ grids with a per-row reference
def p4_grid(J, K=2, Ds=None, kinds=('C2', 'C3H'), pr=None, k3src=None):
    """p4_odd_analysis.grid, re-typed; the reference of row D is dref(J, D)."""
    h = '+-'
    RC = {kd: P4.load(kd, J, k3=(K == 3 if k3src is None else k3src)) for kd in kinds}
    Ds = Ds or sorted(set(k[0] for k in RC['C2'] if k[0] < 1.0e4))
    out = {}
    for D in Ds:
        DR = dref(J, D)
        m = P4.mphys(D, J); mr = P4.mphys(DR, J)
        xs = sorted(set(round(k[1]/m, 6) for k in RC['C2'] if k[0] == D))
        xs = [x for x in xs if all(P4.getrec(RC[kd], D, x*m) is not None and P4.getrec(RC[kd], DR, x*m) is not None for kd in kinds)]
        xs = [x for x in xs if (XMIN <= x <= XMAX) or D < 1990]                  # diagnostic (GSM_XMIN): drop the 2000-row nodes below XMIN
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}
        c0 = None
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G = []
            for x in xs:
                nu = x*m; st = P4.getrec(RC[kd], D, nu); rf = P4.getrec(RC[kd], DR, nu); xr = nu/mr
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                ords = ('M0', 'M2', 'M4', 'M6') if K == 3 else ('M0', 'M2', 'M4')
                v = np.exp(dl)*P4.pred([st[o][-1, -1] for o in ords], [rf[o][-1, -1] for o in ords], K)*fk
                G.append((P4.SIG[h] if k == 3 else 1.0)*v*P4.env(h, xr, k)*P4.Dc(h, J, xr))
            rows[kd] = (k, np.array(G))
        sel = (xs <= P4.XFIT + 1e-9) & (xs*m >= P4.NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]

        def fit(g):
            A_ = np.array([P4.Dc(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A_, (g/P4.env(h, xs, 2))[sel], rcond=None)[0][0]
        a = fit(rows['C2'][1].real)
        a_raw = a
        if P4.DESHIFT:
            sh = {k: P4.shift(J, h, k, xs, m) for k in (2, 3)}
            for _ in range(3):
                a = fit(rows['C2'][1].real - a*sh[2])
            rows = {kd: (k, G - a*sh[k]) for kd, (k, G) in rows.items()}
        rec = dict(m=m, xs=xs, a=a, a_raw=a_raw, c0=c0, dref=DR)
        for kd, (k, G) in rows.items():
            sg = P4.SIG[h] if k == 3 else 1.0
            rec[kd] = m*m*(G.real/(a*P4.env(h, xs, k)) - sg*P4.Dc(h, J, xs))
            rec[(kd, 'im')] = m*G.imag/(a*P4.env(h, xs, k))
        out[D] = rec
    return out


def l9_grid(J, K, kinds=('C2', 'C3H'), XS=None, Ds=None, pr=None):
    """l9_analysis.grid, re-typed; the reference of row D is dref(J, D) (reference helicity basis per Delta_ref)."""
    RC = {k: A.Recs(k, J, K) for k in kinds}
    Rb = dict(A.L.records(J)); Rb.update({kk: v for kk, v in A._load('C2', J, '_k2').items()}); Rb.update(A._load('C2', J, ''))
    hels = ('++', '+-') if J >= 4 else ('++',)
    Ds = Ds or sorted(set(k[0] for k in RC['C2'].R if k[0] < 1.0e4))
    out = {}
    for D in Ds:
        DR = dref(J, D)
        Bref = A.L.basis(J, DR, Rb) if J >= 4 else dict(Ti=np.eye(1))
        m = A.mphys(D, J); mr = A.mphys(DR, J)
        B = A.L.basis(J, D, Rb) if J >= 4 else dict(Ti=np.eye(1))
        if B is None or Bref is None:
            continue
        xs = XS or sorted(set(round(k[1]/m, 6) for k in RC['C2'].R if k[0] == D))
        xs = [x for x in xs if all(RC[kd].get(D, x*m) is not None and RC[kd].get(DR, x*m) is not None for kd in kinds)]
        xs = [x for x in xs if (XMIN <= x <= XMAX) or D < 1990]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}
        c0 = None
        for kd in kinds:
            k = 2 if kd == 'C2' else 3
            G = {h: [] for h in hels}; Go = []
            for x in xs:
                nu = x*m; st = RC[kd].get(D, nu); rf = RC[kd].get(DR, nu); xr = nu/mr
                Ms, Mr = A.helmat(st, B, K), A.helmat(rf, Bref, K)
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                for j, h in enumerate(hels):
                    v = np.exp(dl)*A.pred([Mo[j, j] for Mo in Ms], [Mo[j, j] for Mo in Mr], K)*fk
                    G[h].append((A.SIG[h] if k == 3 else 1.0)*v*A.env(h, xr, k)*A.Dsh(h, J, xr))
                if len(hels) == 2:
                    r2 = Mr if kd == 'C2' else A.helmat(RC['C2'].get(DR, nu), Bref, K)
                    sq = np.sqrt(abs((r2[0][0, 0]/r2[0][1, 1]).real))
                    off = 0.5*np.array([Mo[0, 1] + Mo[1, 0] for Mo in Ms])
                    v = np.exp(dl)*A.pred(off, [Mo[0, 0] for Mo in Mr], K)*fk*sq
                    Go.append(v*np.sqrt(A.env('++', xr, k)*A.env('+-', xr, k)*A.Dsh('++', J, xr)*A.Dsh('+-', J, xr)))
            rows[kd] = (k, {h: np.array(G[h]) for h in hels}, np.array(Go))
        a = {}
        sel = (xs <= A.XFIT + 1e-9) & (xs*m >= A.NUFIT)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]

        def fit(h, g):
            A_ = np.array([A.Dsh(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A_, (g/A.env(h, xs, 2))[sel], rcond=None)[0][0]
        for h in hels:
            a[h] = fit(h, rows['C2'][1][h].real)
        if A.DESHIFT:
            sh = {(k, h): A.shift(J, h, k, xs, m) for k in (2, 3) for h in hels}
            for _ in range(3):
                for h in hels:
                    a[h] = fit(h, rows['C2'][1][h].real - a[h]*sh[(2, h)])
            for kd, (k, G, Go) in rows.items():
                for h in hels:
                    G[h] = G[h] - a[h]*sh[(k, h)]
        rec = dict(m=m, xs=xs, a=a, dref=DR)
        for kd, (k, G, Go) in rows.items():
            for h in hels:
                sg = A.SIG[h] if k == 3 else 1.0
                rec[(kd, h)] = m*m*(G[h].real/(a[h]*A.env(h, xs, k)) - sg*A.Dsh(h, J, xs))
                rec[(kd, h, 'im')] = m*G[h].imag/(a[h]*A.env(h, xs, k))
            if len(hels) == 2:
                eo = np.sqrt(A.env('++', xs, k)*A.env('+-', xs, k))
                rec[(kd, 'off')] = m*Go.real/np.sqrt(a['++']*a['+-'])/eo
                rec[(kd, 'off', 'im')] = m*Go.imag/np.sqrt(a['++']*a['+-'])/eo
        out[D] = rec
    return out


NC._P4grid = p4_grid            # ng_chain's Q.P4.grid lambda calls _P4grid(J, KK, ...)
W._A_grid = l9_grid             # w2_chain's Q.A.grid (_grid: XSUB / DROPX / Hampel) calls _A_grid(J, K, kinds, ...)


# ------------------------------------------------------------------------------------------------ reference-zero rule (per-row reference)
_RB = {}


def refdiag(J, D):
    pats = (['w1_p3_%s_J%d*.pkl', 'ee1000_p3_%s_J%d*.pkl', 'q2_p3_%s_J%d_*.pkl', 'p3_pts_%s_J%d.pkl'] if J % 2 else
            ['w1_pts_%s_J%d_k2*.pkl', 'ee1000_pts_%s_J%d_k2*.pkl', 'l9_pts_%s_J%d_k2.pkl', 'q2_pts_%s_J%d_k2*.pkl', 'u5_pts_%s_J%d_k2*.pkl'])
    m = E.mphys(D, J); ex = {}; DR = dref(J, D)
    for kind in ('EE', 'C3H'):
        R = {}
        for p in pats:
            for d in (NC.SFC, NC.NGC, GSR, GNR, GSMR, GS2R):
                for fn in sorted(glob.glob(os.path.join(d, p % (kind, J)))):
                    R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items()})
        out = []
        for (DD, nu), r in R.items():
            if DD != D or (DR, nu) not in R:
                continue
            rf = R[(DR, nu)]
            a = np.abs(np.asarray(rf['M0'])).max()
            out.append((nu/m, nu, float(np.abs(np.asarray(rf['M2'])).max()/a)))
        if not out:
            continue
        c = float(np.median([v*nu for _, nu, v in out]))
        for x, nu, v in out:
            ex[round(x, 6)] = max(ex.get(round(x, 6), -np.inf), v - c/nu)
    return sorted(ex.items())


def refbad(J, D):
    key = (J, D, dref(J, D))
    if key not in _RB:
        _RB[key] = {x for x, e in refdiag(J, D) if e > REFBAD_TOL}
    return _RB[key]


def ref_zero_x(J, D):
    from scipy.optimize import brentq
    h = '+-' if J % 2 else '++'
    if J < 1:
        return np.inf
    f = lambda t: E.Dc_fixed(h, J, t)
    ts = np.linspace(1e-5, 0.6, 12000); v = f(ts)
    i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if not len(i):
        return np.inf
    return brentq(f, ts[i[0]], ts[i[0] + 1])*E.mphys(dref(J, D), J)/E.mphys(D, J)


def _drop(r, xbad):
    n = len(r['xs']); sel = np.array([round(float(x), 6) not in xbad for x in r['xs']])
    if sel.all():
        return r
    return {k: (v[sel] if isinstance(v, np.ndarray) and v.shape[:1] == (n,) else v) for k, v in r.items()}


DROPLOG = {}


def _apply_refbad(J, g):
    if not REFBAD:
        return g
    out = {}
    for D, r in g.items():
        mine = D > 1990 or (J % 2 and J >= 17 and abs(D - (1000 + J)) < 1e-6)
        xb = set(refbad(J, D)) if mine else set()
        if mine:
            n = len(r['xs'])
            for k, v in r.items():
                if isinstance(v, np.ndarray) and v.shape[:1] == (n,) and np.issubdtype(v.dtype, np.number):
                    bad = ~np.isfinite(v.reshape(n, -1)).all(1)
                    xb |= {round(float(x), 6) for x in np.asarray(r['xs'])[bad]}
        if xb:
            DROPLOG[(J, D)] = sorted(x for x in xb if round(x, 6) in {round(float(v), 6) for v in r['xs']})
        out[D] = _drop(r, xb) if xb else r
    return out


def odd_grid(J):
    g = Q.P4.grid(J, K=2, kinds=('C2', 'C3H'))
    Ds = Q.row_filter(J, sorted(D for D in g if D <= 6000.0))
    return _apply_refbad(J, {D: g[D] for D in Ds})


def even_grid(J):
    g = Q.A.grid(J, 2, ('C2', 'C3H'))
    Ds = Q.row_filter(J, sorted(D for D in g if D <= 6000.0))
    return Q.align_rows(J, _apply_refbad(J, {D: g[D] for D in Ds}))


assert not Q.DIRECT400
Q.odd_grid = odd_grid
Q.even_grid = even_grid


def pr(s='', buf=None):
    print(s, flush=True)
    if buf is not None:
        buf.append(s)


# ------------------------------------------------------------------------------------------------ run_state (ng_bind's, + OO block)
def errdiag(d, sel, G, i=0):
    h = d['hs'][i]; a2, a3 = d['abs'][h]
    Gv = np.broadcast_to(G, d['nu'].shape)
    return ((a2 + a3)*Gv)[sel].sum()/d['Ns'][i]


def lam_pp(d, sel=None):
    """normalised ++ diagonal (base_++ + dF_++)/N_s,++ (= lambda_OO/N_s given Gamma_OO = Gamma_++)."""
    _, dF = NC.lam_of(d, sel=sel)
    return float((d['base'][0] + dF[0, 0])/d['Ns'][0])


def run_state(Dg, J, Down, mfl, out_pr, npn=256):
    Lm = E.lam_of(Dg); m = mfl*Lm
    out = dict(J=J, mfl=mfl, Down=Down)
    NC.PK = None
    g2, dsp2 = NB.grid_for(J, 2, False)
    g1, dsp1 = NB.grid_for(J, 1, False)
    g0m, _ = NB.grid_for(J, 2, True)
    f2 = NB.model_of(J, g2); f0m = NB.model_of(J, g0m); f1 = NB.model_of(J, g1)
    d2 = NC.dens(Dg, J, mfl, f2, npn); dm = NC.dens(Dg, J, mfl, f0m, npn)
    lam2, dF2 = NC.lam_of(d2); lamm, dFm = NC.lam_of(dm)
    lam0 = E.lam2(d2['base'], d2['Ns'], 0*dF2)
    d1 = NC.dens(Dg, J, mfl, f1, npn)
    lamK1 = NC.lam_of(d1)[0]
    d64 = NC.dens(Dg, J, mfl, f2, 64); lam64 = NC.lam_of(d64)[0]
    nu, x = d2['nu'], d2['x']
    own = g2.get(float(Down))
    xs = own['xs'] if own is not None else np.array([0.05, 0.8])
    x0, xmax = xs[0], xs[-1]
    hs = d2['hs']
    d2['abs'] = NB.absdens(Dg, J, mfl, npn)
    tr = {}
    for kd in ('C2', 'C3H'):
        for h in (hs + (['off'] if len(hs) == 2 else [])):
            key = kd if J % 2 else (kd, h)
            gd = {}
            for D in g2:
                if D not in g1:
                    continue
                x2r, x1r = np.round(g2[D]['xs'], 6), np.round(g1[D]['xs'], 6)
                xc = np.intersect1d(x2r, x1r)
                if len(xc) < 2:
                    continue
                v2 = np.array([g2[D][key][int(np.argmin(np.abs(x2r - v)))] for v in xc])
                v1 = np.array([g1[D][key][int(np.argmin(np.abs(x1r - v)))] for v in xc])
                gd[D] = dict(m=g2[D]['m'], xs=xc, v=np.abs(v2 - v1))
            if isinstance(key, tuple) and key[1] == 'off':
                gd = {D: dict(v, v=v['v']*g2[D]['m']) for D, v in gd.items()}
            tr[(kd, h)] = NB.interp_nodes(gd, J, m, x, lambda r: r['v'])
    n = len(hs); Mt = np.zeros((n, n))
    win = (x >= x0) & (x <= xmax)
    for i, h in enumerate(hs):
        a2, a3 = d2['abs'][h]
        Mt[i, i] = (a2*tr[('C2', h)] + a3*tr[('C3H', h)])[win].sum()/d2['Ns'][i]
    if n == 2:
        a2, a3 = d2['abs']['off']
        Mt[0, 1] = Mt[1, 0] = (a2*tr[('C2', 'off')] + a3*tr[('C3H', 'off')])[win].sum()/np.sqrt(d2['Ns'].prod())
    E_trunc = 3*np.sqrt((Mt*Mt).sum())
    fl = NB.lin_model(J, g2); dl_ = NC.dens(Dg, J, mfl, fl, npn); laml = NC.lam_of(dl_)[0]
    lsub, dsubs = [], []
    for keep in (0, 1):
        gs = dict(g2)
        if own is not None:
            gs[float(Down)] = NB.sub({float(Down): own}, keep)[float(Down)]
        ds_ = NC.dens(Dg, J, mfl, NB.model_of(J, gs), npn); dsubs.append(ds_)
        lsub.append(NC.lam_of(ds_)[0])
    E_interp = 3*max(abs(lam2 - laml), *[abs(lam2 - v) for v in lsub])
    lowsel = nu < max(NB.NU_LOW, x0*m)
    Glow = NB.G_LOW*gauge_factor(Dg, nu)
    E_low = NB.errmat(d2, lowsel, Glow)
    gedge = 1e4
    if own is not None:
        for kd in ('C2', 'C3H'):
            for h in hs:
                key = kd if J % 2 else (kd, h)
                gedge = max(gedge, 10*np.max(np.abs(np.asarray(own[key])[-4:])))
    edsel = x > xmax
    E_edge = NB.errmat(d2, edsel, gedge)
    csel = ~(lowsel | edsel)
    lamc = NC.lam_of(d2, sel=csel)[0]
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
    E_tl = 3*abs(lam2 - lamK1)
    tot_l = E_tl + E_interp + E_low + E_edge
    out.update(lam_hybrid=float(lam0), lam_model=float(lamm), lam_K2=float(lam2), lam_K2_64=float(lam64), lam_central=float(lamc),
               lam_lin=float(laml), lam_sub=[float(v) for v in lsub], E_trunc=float(E_trunc), E_interp=float(E_interp),
               E_low=float(E_low), E_edge=float(E_edge), G_edge=float(gedge), E_tot=float(tot), lower=float(lamc - tot),
               L=float(lamc/tot) if tot > 0 else np.inf, lam_K1=float(lamK1), E_trunc_lam=float(E_tl), E_tot_lam=float(tot_l),
               lower_lam=float(lamc - tot_l), L_lam=float(lamc/tot_l) if tot_l > 0 else np.inf, base=base_info,
               rows=sorted(float(D) for D in g2), rows_model=sorted(float(D) for D in g0m), drefs={f"{D:g}": float(g2[D].get('dref', DREF0)) for D in g2},
               despike=dsp2, dF_model=dFm.tolist(), dF_K2=dF2.tolist(), x0=float(x0), xmax=float(xmax),
               xzero=float(ref_zero_x(J, Down)))
    out_pr(f"  J={J:2d} m/Lam={mfl:.5f} rows {out['rows']} (Dref {out['drefs']}; model {out['rows_model']}): hybrid {lam0:+.5f} | "
           f"model chain {lamm:+.5f} | own row K2 {lam2:+.5f} (64 nodes {lam64:+.5f}; lin {laml:+.5f}; sub {lsub[0]:+.5f}/{lsub[1]:+.5f}) | "
           f"central (windows cut) {lamc:+.5f}")
    out_pr(f"        errors: trunc {E_trunc:.2e}  interp {E_interp:.2e}  low(nu<{max(NB.NU_LOW, x0*m):.0f}, G={NB.G_LOW:.0e}) {E_low:.2e}  "
           f"edge(x>{xmax:.3f}, G={gedge:.2e}) {E_edge:.2e}  -> total {tot:.2e}; lower bound {lamc - tot:+.5f}; L = {out['L']:.1f}")
    out_pr(f"        lambda-level truncation: K1 chain {lamK1:+.5f} -> 3|K2 - K1| = {E_tl:.2e}; total {tot_l:.2e}; lower bound "
           f"{lamc - tot_l:+.5f}; L = {out['L_lam']:.1f}   (ref zero x0 = {out['xzero']:.3f})")
    for h, b in base_info.items():
        out_pr(f"        base {h}: hybrid {b['lam_hyb']:+.5f}, flat {b['flat']:+.5f} (Nf/Ns {b['Nf_Ns']:.4f}), LR curvature {b['cv_Ns']:+.2e} N_s, "
               f"Lemma-C bound {b['cvb_Ns']:.2e} N_s, certified-floor route (Nf/5 - LemmaC)/Ns = {b['rig']:+.5f}")
    if dsp2:
        out_pr(f"        despiked: {dsp2}")
    if J % 2 == 0 and J >= 4:              # OO block (1x1) = normalised ++ diagonal, Gamma_OO = Gamma_++
        oc = lam_pp(d2, csel); o2 = lam_pp(d2); o1 = lam_pp(d1); ol = lam_pp(dl_); os_ = [lam_pp(v) for v in dsubs]
        o_hyb = float(d2['base'][0]/d2['Ns'][0])
        Et = 3*Mt[0, 0]; Etl = 3*abs(o2 - o1)
        Ei = 3*max(abs(o2 - ol), *[abs(o2 - v) for v in os_])
        El = errdiag(d2, lowsel, Glow); Ee = errdiag(d2, edsel, gedge)
        Eid = errdiag(d2, np.ones_like(nu, bool), G_ID)
        totO = Et + Ei + El + Ee + Eid; totOl = Etl + Ei + El + Ee + Eid
        out['OO'] = dict(lam_hyb=o_hyb, lam_K2=o2, lam_central=oc, lam_K1=o1, E_trunc=float(Et), E_trunc_lam=float(Etl), E_interp=float(Ei),
                         E_low=float(El), E_edge=float(Ee), E_id=float(Eid), E_tot=float(totO), E_tot_lam=float(totOl),
                         lower=float(oc - totO), lower_lam=float(oc - totOl), L=float(oc/totO), L_lam=float(oc/totOl))
        out_pr(f"        OO (= ++ diagonal, Gamma_OO = Gamma_++): hybrid {o_hyb:+.5f}  K2 {o2:+.5f}  central {oc:+.5f} | trunc {Et:.2e} "
               f"(lam {Etl:.2e}) interp {Ei:.2e} low {El:.2e} edge {Ee:.2e} id(G={G_ID:g}) {Eid:.2e} -> L = {oc/totO:.1f} (lam {oc/totOl:.1f})")
    return out


# ------------------------------------------------------------------------------------------------ H_mono-A tables
def hmono(J, specs, buf, xs=(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9)):
    """specs: list of (D, dref); each row built with its own reference."""
    HIDE[0] = False; NC.NOOWN = False
    res = {}
    gg = []; hdr = []
    for D, DR in specs:
        REFMAP.clear(); REFMAP[(J, float(D))] = DR
        g, _ = NB.grid_for(J, 2, False); g1, _ = NB.grid_for(J, 1, False)
        gg.append((D, DR, g.get(float(D)), g1.get(float(D))))
        if g.get(float(D)) is not None:
            hdr.append(f"{D:g}/ref {DR:g} (J^2/m={J*J/g[float(D)]['m']:.3f}, x0={ref_zero_x(J, D):.3f}, dropped {DROPLOG.get((J, float(D)), [])})")
        REFMAP.clear()
    pr(f"## H_mono-A, J = {J}: Gamma(x) in c2 at K = 2 [|K2 - K1|], rows " + ", ".join(hdr), buf)
    keys = ['C2', 'C3H'] if J % 2 else ([('C2', '++'), ('C3H', '++')] if J < 4 else
                                        [('C2', '++'), ('C2', '+-'), ('C3H', '++'), ('C3H', '+-'), ('C2', 'off'), ('C3H', 'off')])
    for key in keys:
        pr(f"   {str(key):16s} x: " + "  ".join(f"{D:>9g}/{DR:<6g}" for D, DR, _, _ in gg), buf)
        for x in xs:
            row = []
            for D, DR, r, r1 in gg:
                if r is None or key not in r or not (r['xs'][0] - 1e-9 <= x <= r['xs'][-1] + 1e-9):
                    row.append(f"{'—':>16s}"); continue
                fac = r['m'] if (isinstance(key, tuple) and key[1] == 'off') else 1.0
                v = np.interp(x, r['xs'], r[key])*fac
                s = (abs(np.interp(x, r1['xs'], r1[key])*fac - v) if (r1 is not None and key in r1) else np.nan)
                res.setdefault(str(key), {}).setdefault(f"{D:g}/{DR:g}", []).append((x, float(v), float(s)))
                row.append(f"{v:+9.1f} [{s:5.0f}]")
            pr(f"   {'':16s} {x:.2f} " + "  ".join(row), buf)
    return res


def mfl_of(tok, Dg, J):
    if tok == 'tw':
        return float(E.mgap(Dg, J))
    if tok == 'tw1000':
        return float(E.mgap(1000, J))
    return float(tok)


def run(specs, buf):
    out = []
    for spec in specs:
        parts = spec.split(':')
        Dg, J, Down, mf = parts[:4]
        Dg, J, Down = int(Dg), int(J), float(Down)
        REFMAP.clear()
        if len(parts) > 4:
            REFMAP[(J, Down)] = float(parts[4])
        HIDE[0] = (Dg == 1000)
        pr(f"\n# Delta_g = {Dg}, J = {J}: own row Delta = {Down:g} (Dref {dref(J, Down):g}); twist edge m/Lam = {E.mgap(Dg, J):.7f}; "
           f"rows hidden: {HIDE[0]}; Lam = {E.lam_of(Dg):.3f}", buf)
        for tok in mf.split(','):
            mfl = mfl_of(tok, Dg, J)
            r = run_state(Dg, J, Down, mfl, lambda s: pr(s, buf))
            r.update(Dg=Dg, tok=tok, dropped={f"{k[1]:g}": v for k, v in DROPLOG.items() if k[0] == J})
            out.append(r)
    HIDE[0] = False; REFMAP.clear()
    return out


if __name__ == '__main__':
    t0 = time.time()
    buf = []
    if sys.argv[1] == '--hmono':
        tag = 'hmono_' + sys.argv[2]
        res = {}
        for spec in sys.argv[3:]:
            Js, rows = spec.split(':')
            J = int(Js)
            res[J] = hmono(J, [(float(a.split('/')[0]), float(a.split('/')[1])) for a in rows.split(',')], buf)
        json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'cache', 'gs2_%s.json' % tag), 'w'))
        open(os.path.join(HERE, 'cache', 'gs2_%s.txt' % tag), 'w').write("\n".join(buf) + "\n")
    else:
        tag = sys.argv[1]
        res = run(sys.argv[2:], buf)
        json.dump(res, open(os.path.join(HERE, 'cache', 'gs2_bind_%s.json' % tag), 'w'), default=float)
        open(os.path.join(HERE, 'cache', 'gs2_bind_%s.txt' % tag), 'w').write("\n".join(buf) + "\n")
    pr(f"[{time.time() - t0:.0f}s]", buf)
