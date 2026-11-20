"""
ee1000_chain (EE_DG1000_LARGEJ.md): per-state positivity of Omega^EE at Delta_g = 1000 for the large-J near-gap states (even J >= 32,
odd J >= 17), with the K2 crystal rows at the states' own masses (ee1000_rows) -- the w2/q3 chain (EE_NEARGAP.md) restricted to hybrid
states (flat NF02 + leading-Regge curvature, scalar remainder 0, as q3 Q3_EXTRA_M), with three changes, each switchable:
  (1) EE1000_JACFIX=1 (default): stable complex Jacobi polynomial in the +- de-shift (ee1000_common; the EE_NEARGAP "K2 +- breakdown");
  (2) EE1000_MODEL = spline (q3: natural cubic spline in x, held beyond the last node; 1/m-linear between rows, held outside)
                   | shape (default): per row and key, least-squares fit
                         Gamma(x) = c(x) D_h(x) + s(x) x D_h'(x)/(2J+1) + r(x),   c, s, r polynomials of degree EE1000_K in x
                     (D_h the LO crystal shape of that helicity).  The oscillation of Gamma is carried by the exact D_h, D_h', so the
                     fit is free of the aliasing that a spline of an oscillating Gamma suffers (EE_NEARGAP §5: P_J period vs node
                     spacing); c, s, r are held at the last node beyond it (the oscillation continues with D_h);
  (3) EE1000_NPN (default 256): Gauss-Legendre nodes in p for the curvature and correction sums (q3: 64).
States: for each J, m/Lam in {m_gap(J), 1.05, 1.075, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2, 2.5, 3, 4, 5} (twist-allowed only; m_gap is the
lowest twist-allowed mass, Delta = Delta_g + J); even J: 2x2 (++, +-), odd J: 1x1 (+-).  Off-diagonal: spline of c_off = m Gamma_off
(q3 OFFSCALE), EE-own helicity basis, S1_S4 row alignment.  Hampel despike as w2.
Output: cache/ee1000_chain_<tag>.txt, cache/ee1000_recs_<tag>.pkl (per state: base, N_s, dF, lam).
Usage: EE1000_JEVEN=32,36,40 EE1000_JODD=17,31 [EE1000_MODEL=shape|spline] [W2_ROWS=1040,...] EE1000_TAG=x python3 ee1000_chain.py
"""
import os, sys, io, glob, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import ee1000_common as E

Q, W, C = E.Q, E.W, E.C
JEVEN = [int(v) for v in os.environ.get('EE1000_JEVEN', '').split(',') if v]
JODD = [int(v) for v in os.environ.get('EE1000_JODD', '').split(',') if v]
MODEL = os.environ.get('EE1000_MODEL', 'shape')
KDEG = int(os.environ.get('EE1000_K', '4'))
NPN = int(os.environ.get('EE1000_NPN', '256'))
TAG = os.environ.get('EE1000_TAG', 'x')
DG = int(os.environ.get('EE1000_DG', '1000'))
MGRID = [1.05, 1.075, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.75, 2.0, 2.5, 3.0, 4.0, 5.0]
GSCALE = float(os.environ.get('EE1000_GSCALE', '1'))          # common factor on all corrections (t* scans)
MPOW = float(os.environ.get('EE1000_MPOW', '0'))              # Gamma(x; m) = Gamma_row(x) (m_row/m)^MPOW above the largest row (0 = held)
SHAPEMIN = int(os.environ.get('EE1000_SHAPEMIN', '30'))       # rows with fewer nodes: spline
XMAX = float(os.environ.get('EE1000_XMAX', '0.8'))            # row nodes used: x <= XMAX (K2 validated to 0.8; beyond: extrapolated)
XMIN = float(os.environ.get('EE1000_XMIN', '0'))              # and x >= XMIN


def trim(g):
    out = {}
    for D, r in g.items():
        if D >= 1e4:
            continue
        sel = (r['xs'] <= XMAX + 1e-9) & (r['xs'] >= XMIN - 1e-9)
        out[D] = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == sel.shape else v) for k, v in r.items()}
    return out


# ------------------------------------------------------------------------------------------------ loaders: + ee1000 caches
def ee_even(J):
    return W._pk(['u5_pts_EE_J%d_k2*.pkl' % J, 'w1_pts_EE_J%d_k2*.pkl' % J, 'ee1000_pts_EE_J%d_k2*.pkl' % J])


W.ee_even = ee_even
_q2e = Q._q2_even


def _q2_even(kind, J):
    R = _q2e(kind, J)
    R.update(W._pk(['ee1000_pts_%s_J%d_k2*.pkl' % (kind, J)]))
    return R


Q._q2_even = _q2_even


def _odd(kind, J, anti=False, k3=False):
    if kind == 'C2':
        R = W._pk(['w1_p3_EE_J%d*.pkl' % J, 'ee1000_p3_EE_J%d.pkl' % J])
    else:
        R = Q._orig_p4load_raw(kind, J, k3) if J < 11 else {}
        R.update(W._pk(['q2_p3_%s_J%d_*.pkl' % (kind, J), 'w1_p3_%s_J%d*.pkl' % (kind, J), 'ee1000_p3_%s_J%d.pkl' % (kind, J)]))
    R = {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}
    return Q._strip(R, J)


Q.P4.load = _odd


# ------------------------------------------------------------------------------------------------ direct state values (EE1000_DSTATE=1)
# state record := direct all-orders matrix (ee1000_direct, TPS normalisation), M2 = M4 = 0; reference K2 (w2 W2_DSTATE convention).
DSTATE = os.environ.get('EE1000_DSTATE', '0') == '1'
DS_USED = []
if DSTATE:
    import glob as _g
    import n1_analysis as N1
    import ee1000_direct as ED
    for _fn in sorted(_g.glob(os.path.join(C, '*shadow_k2*.pkl'))):
        N1.SHK2.update(pickle.load(open(_fn, 'rb')))

    def _dstate(R, J, kind):
        db = ED.dload(J, kind)
        R = dict(R)
        I2 = dict(Ti=np.eye(2 if J >= 4 else 1))
        for (D, nu), r in db.items():
            k = (float(D), round(float(nu), 9))
            if D >= 1e4 or k not in R:
                continue
            H, l = N1.Hrep(r, I2, D, J)
            st = R[k]
            M = np.exp(l - (st['logC'] - st['logK']))*H.mean(0)
            R[k] = dict(st, M0=M, M2=0*M, M4=0*M)
            DS_USED.append((kind, J, D, round(nu/st['m'], 4)))
        return R
    W._dstate = _dstate
    W.DSTATE = True
    _q2o = W._q2_even_orig
    W._q2_even_orig = lambda kind, J: {**_q2o(kind, J), **W._pk(['ee1000_pts_%s_J%d_k2*.pkl' % (kind, J)])}


# ------------------------------------------------------------------------------------------------ shape model
def _dshape(h, J, x):
    D = E.Dc_fixed(h, J, x)
    hh = 1e-6
    Dp = (E.Dc_fixed(h, J, x + hh) - E.Dc_fixed(h, J, x - hh))/(2*hh)
    return D, x*Dp/(2*J + 1)


def _basis(h, J, x, xm, K):
    D, Dp = _dshape(h, J, x)
    t = np.minimum(x, xm)/xm
    P = [t**k for k in range(K + 1)]
    return np.array([p*D for p in P] + [p*Dp for p in P] + [p for p in P]).T


class ShapeFit:
    """Gamma(x) = c(x) D_h + s(x) x D_h'/(2J+1) + r(x) on one row; keys (kind, h) (even) or kind (odd, h = '+-')."""
    def __init__(self, J, row, keys, K=KDEG):
        self.J, self.K, self.xm = J, K, row['xs'][-1]
        self.coef, self.res = {}, {}
        x = row['xs']
        for key in keys:
            h = key[1] if isinstance(key, tuple) else '+-'
            A = _basis(h, J, x, self.xm, K)
            y = np.asarray(row[key], float)
            sc = np.abs(A).max(0); sc[sc == 0] = 1
            c, *_ = np.linalg.lstsq(A/sc, y, rcond=None)
            c = c/sc
            self.coef[key] = (h, c)
            r = y - A @ c
            self.res[key] = (np.sqrt(np.mean(r**2)), np.sqrt(np.mean(y**2)), np.abs(r).max())

    def __call__(self, key, x):
        h, c = self.coef[key]
        return _basis(h, self.J, np.asarray(x, float), self.xm, self.K) @ c


def row_models(J, g, keys):
    """per row D: callable key -> Gamma(x) (shape fit for the diagonal keys; spline of c_off = m Gamma_off for 'off')."""
    from scipy.interpolate import CubicSpline
    mods = {}
    for D, r in g.items():
        # shape fit only on dense rows (3(K+1) parameters per key; the 16-19-node rows of EE_NEARGAP are resolved by the spline at J <= 28)
        sf = (ShapeFit(J, r, [k for k in keys if not (isinstance(k, tuple) and k[1] == 'off')])
              if MODEL == 'shape' and len(r['xs']) >= SHAPEMIN else None)
        spl = {k: CubicSpline(r['xs'], r[k], bc_type='natural') for k in keys}
        xm = r['xs'][-1]

        def f(key, x, sf=sf, spl=spl, xm=xm, m=r['m']):
            x = np.asarray(x, float)
            if isinstance(key, tuple) and key[1] == 'off':
                return spl[key](np.minimum(x, xm))*m          # c_off = m Gamma_off (held beyond the last node)
            if sf is not None:
                return sf(key, x)
            return spl[key](np.minimum(x, xm))
        mods[D] = (r['m'], f, sf)
    return mods


def model_eval(mods, key, x, m):
    """1/m-linear between rows, held outside (MPOW: power-law in m above the largest row)."""
    Ds = sorted(mods, key=lambda D: mods[D][0]); ms = np.array([mods[D][0] for D in Ds])
    if m <= ms[0]:
        return mods[Ds[0]][1](key, x)
    if m >= ms[-1]:
        return mods[Ds[-1]][1](key, x)*(ms[-1]/m)**MPOW
    i1 = int(np.searchsorted(ms, m)); i0 = i1 - 1
    t = (1/ms[i0] - 1/m)/(1/ms[i0] - 1/ms[i1])
    return (1 - t)*mods[Ds[i0]][1](key, x) + t*mods[Ds[i1]][1](key, x)


# ------------------------------------------------------------------------------------------------ states
def mlist(J):
    mg = E.mgap(DG, J)
    return [mg] + [v for v in MGRID if v > mg + 1e-6]


def even_states(J, mods, pr):
    recs = []
    for mfl in mlist(J):
        st = {h: E.hyb(DG, J, mfl, h, NPN) for h in ('++', '+-')}
        m = st['++']['m']; x = st['++']['x']
        G = {}
        for kd in ('C2', 'C3H'):
            for e in ('++', '+-', 'off'):
                G[(kd, e)] = model_eval(mods, (kd, e), x, m) if mods else 0*x
        G[('C2', 'off')] = G[('C2', 'off')]/m; G[('C3H', 'off')] = G[('C3H', 'off')]/m     # c_off -> Gamma_off (q3 OFFSCALE)
        base = np.array([st[h]['base'] for h in ('++', '+-')]); Ns = np.array([st[h]['Ns'] for h in ('++', '+-')])
        dF = np.zeros((2, 2))
        for i, h in enumerate(('++', '+-')):
            dF[i, i] = E.dF_state(st[h], G, h)
        dF[0, 1] = dF[1, 0] = E.dF_off(st['++'], st['+-'], G)
        dF = GSCALE*dF
        recs.append(dict(J=J, mfl=mfl, hs=['++', '+-'], base=base, Ns=Ns, dF=dF, lam=E.lam2(base, Ns, dF), lam0=E.lam2(base, Ns, 0*dF),
                         diag=[dF[i, i]/Ns[i] for i in range(2)], off=abs(dF[0, 1])/np.sqrt(Ns.prod()),
                         curv=[(st[h]['cv2'] + st[h]['cv3'])/st[h]['Ns'] for h in ('++', '+-')]))
    return recs


def odd_states(J, mods, pr):
    recs = []
    for mfl in mlist(J):
        st = E.hyb(DG, J, mfl, '+-', NPN)
        m, x, u = st['m'], st['x'], st['u']
        g2 = model_eval(mods, 'C2', x, m) if mods else 0*x
        g3 = model_eval(mods, 'C3H', x, m) if mods else 0*x
        d = (st['p2']*Q.P4.env('+-', x, 2)*g2/m**2).sum() + (st['p3']*u*Q.P4.env('+-', x, 3)*g3/m**2).sum()
        d = GSCALE*d
        base, Ns = np.array([st['base']]), np.array([st['Ns']])
        recs.append(dict(J=J, mfl=mfl, hs=['+-'], base=base, Ns=Ns, dF=np.array([[d]]), lam=(st['base'] + d)/st['Ns'],
                         lam0=st['base']/st['Ns'], diag=[d/st['Ns']], off=0.0, curv=[(st['cv2'] + st['cv3'])/st['Ns']]))
    return recs


def run(pr):
    allrecs, info = [], {}
    for J in JEVEN:
        g = trim(Q.even_grid(J))
        keys = [(kd, e) for kd in ('C2', 'C3H') for e in ('++', '+-', 'off')]
        mods = row_models(J, g, keys) if g else {}
        info[J] = {D: dict(m=r['m'], n=len(r['xs']), xs=r['xs'], res={k: v for k, v in (mods[D][2].res.items() if mods[D][2] else [])})
                   for D, r in g.items()}
        recs = even_states(J, mods, pr)
        allrecs += recs
        report(J, g, mods, recs, pr)
    for J in JODD:
        g = trim(Q.odd_grid(J))
        mods = row_models(J, g, ['C2', 'C3H']) if g else {}
        info[J] = {D: dict(m=r['m'], n=len(r['xs']), xs=r['xs'], res={k: v for k, v in (mods[D][2].res.items() if mods[D][2] else [])})
                   for D, r in g.items()}
        recs = odd_states(J, mods, pr)
        allrecs += recs
        report(J, g, mods, recs, pr)
    return allrecs, info


def report(J, g, mods, recs, pr):
    rows = ", ".join(f"{D:.0f} (m={r['m']:.1f}, {len(r['xs'])} nodes, x <= {r['xs'][-1]:.3f}, J^2/m={J*J/r['m']:.2f})" for D, r in sorted(g.items()))
    pr(f"  J={J:2d} rows: {rows or 'NONE (uncorrected)'}")
    for D in sorted(mods):
        sf = mods[D][2]
        if sf is not None:
            pr(f"        shape-fit residual rms/rms(Gamma) [max|res|] at {D:.0f}: " +
               "; ".join(f"{k if isinstance(k, str) else k[0] + ' ' + k[1]}: {a/b:.3f} [{c:.0f}]" for k, (a, b, c) in sf.res.items()))
    w = min(recs, key=lambda r: r['lam'])
    dg = max(abs(d) for r in recs for d in r['diag'])
    pr(f"      J={J:2d}: min lam/N_s {w['lam']:+.4f} (uncorr. {min(r['lam0'] for r in recs):+.4f}) at m/Lam={w['mfl']:.4f}; "
       f"max |dF_hh|/N_h {dg:.4f}; max |dF_+-|/sqrt(NN) {max(r['off'] for r in recs):.4f}")
    pr("         per m/Lam: " + "  ".join(f"{r['mfl']:.3f}: {r['lam']:+.3f} [" + "/".join(f"{d:+.3f}" for d in r['diag']) + "]" for r in recs))


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    pr(f"# ee1000_chain: Delta_g = {DG}; model {MODEL} (K = {KDEG}); NPN = {NPN}; JACFIX {E.JACFIX}; rows filter W2_ROWS={os.environ.get('W2_ROWS', '')}; "
       f"despike {os.environ.get('W2_DESPIKE', '1')}; GSCALE {GSCALE}; MPOW {MPOW}; x-range [{XMIN}, {XMAX}]")
    recs, info = run(pr)
    w = min(recs, key=lambda r: r['lam'])
    pr(f"# overall (these J): min lam/N_s = {w['lam']:+.4f} at J={w['J']}, m/Lam={w['mfl']:.4f}; states {len(recs)}; all positive: "
       f"{all(r['lam'] > 0 for r in recs)}")
    pr(f"# despiked nodes (J, Delta, x): {W.DESPIKE_LOG}")
    pr(f"# direct state values used (EE1000_DSTATE): {sorted(set(DS_USED))}")
    pr("# row alignment: " + "; ".join(f"{J},{D:.0f}: {ov:+.2f}{'' if sD > 0 else ' FLIP'}" for (J, D), (ov, sD) in sorted(Q.ALIGN_LOG.items())))
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'ee1000_chain_%s.txt' % TAG), 'w').write(buf.getvalue())
    pickle.dump(dict(recs=recs, info=info, env={k: v for k, v in os.environ.items() if k.startswith(('EE1000', 'W2_', 'Q3_'))}),
                open(os.path.join(C, 'ee1000_recs_%s.pkl' % TAG), 'wb'))
