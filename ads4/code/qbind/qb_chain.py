"""
qb_chain (QBIND_COROLLARY.md; copy of neargap2/ng_chain.py, row loaders also read neargap2/cache (read-only) and qbind/cache): the EE_NEARGAP main chain (w2_chain in EE mode, through ee1000_common: Jacobi fix, EE helicity basis,
Hampel despike, local de-shift, hybrid states) run UNCHANGED in-process, with two redirections only:
  * every row loader (w2_chain._pk) also reads neargap2/cache (the own-mass K2 rows of ng_rows.py);
  * outputs go to neargap2/cache/ng_chain_<tag>.{txt,pkl}.
Options (env): Q3_JEVEN / Q3_JODD / Q3_EXTRA_M / W2_ROWS as for w2_chain; NG_K = Laplace order of the K2 rows (default 2; 1 and 0 for
the truncation estimate).
Helpers for the analyses: models(grids) (= q3_analysis.run's model construction), split(...) (per-x-window dF of a hybrid state).
Usage: [env] python3 ng_chain.py Dg tag
"""
import os, sys, io, glob, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import qb_guard  # noqa: F401
os.chdir(SF)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import ee1000_common as E
W, Q = E.W, E.Q
NGC = os.path.join(HERE, 'cache')
NG2 = os.path.join(HERE, '..', 'neargap2', 'cache')       # read-only: NEARGAP2's own-mass rows
SFC = os.path.join(SF, 'cache')
KK = int(os.environ.get('NG_K', '2'))
NOOWN = os.environ.get('NG_NOOWN', '0') == '1'       # 1: ignore the neargap2 rows (reproduces the published chain)


def _pk(pats):
    R = {}
    for p in pats:
        fns = sorted(glob.glob(os.path.join(SFC, p)))
        if p.startswith('w1_p3_'):              # the odd ee1000 rows (Delta = 1000 + J, EE_DG1000_LARGEJ) carry the 1000 near-gap states
            fns += sorted(glob.glob(os.path.join(SFC, 'ee1000_p3_' + p[len('w1_p3_'):])))
        if not NOOWN:
            fns += sorted(glob.glob(os.path.join(NG2, p))) + sorted(glob.glob(os.path.join(NGC, p)))
        for fn in fns:
            R.update(pickle.load(open(fn, 'rb')))
    return R


W._pk = _pk
_Agrid = Q.A.grid
Q.A.grid = lambda J, K, *a, **k: _Agrid(J, KK, *a, **k)
PK = None                                             # even J: Laplace truncation inside pred (Recs loading needs K = 2)
_Apred = Q.A.pred
Q.A.pred = lambda a, r, K: _Apred(a, r, K if PK is None else PK)
_P4grid = Q.P4.grid
Q.P4.grid = lambda J, K=2, *a, **k: _P4grid(J, KK, *a, **k)


def models(grids):
    """the model functions exactly as q3_analysis.run builds them."""
    emods, omods = {}, {}
    for J, g in grids['even'].items():
        if Q.OFFSCALE and J >= 4:
            g = {D: {**r, **{(kd, 'off'): r[(kd, 'off')]*r['m'] for kd in ('C2', 'C3H')}} for D, r in g.items()}
        if len(g) == 1:
            D = list(g)[0]; g = {D: g[D], D + 1e-6: dict(g[D], m=g[D]['m']*(1 + 1e-9))}
        emods[J] = Q.NG.model(g, J, ('C2', 'C3H'))
    for J, g in grids['odd'].items():
        if len(g) == 1:
            D = list(g)[0]; g = {D: g[D], D + 1e-6: dict(g[D], m=g[D]['m']*(1 + 1e-9))}
        omods[J] = Q.odd_model(g)
    return emods, omods


def dens(Dg, J, mfl, f, npn=256):
    """per-node correction densities of a hybrid state: dict h -> (nu, d2, d3) and 'off' -> (nu, o2, o3); plus base, Ns."""
    Lm = E.lam_of(Dg)
    hs = ['+-'] if J % 2 else (['++'] if J < 4 else ['++', '+-'])
    st = {h: E.hyb(Dg, J, mfl, h, npn) for h in hs}
    s0 = st[hs[0]]; x, m, u, y = s0['x'], s0['m'], s0['u'], s0['y']
    nu = s0['pn']*Lm
    G = f(x, m)
    out = dict(hs=hs, nu=nu, x=x, base=np.array([st[h]['base'] for h in hs]), Ns=np.array([st[h]['Ns'] for h in hs]))
    if J % 2:
        env = Q.P4.env
        out['+-'] = (st['+-']['p2']*env('+-', x, 2)*G['C2']/m**2, st['+-']['p3']*u*env('+-', x, 3)*G['C3H']/m**2)
    else:
        G = dict(G)
        if len(hs) == 2:
            G[('C2', 'off')] = G[('C2', 'off')]/m; G[('C3H', 'off')] = G[('C3H', 'off')]/m
            eo2 = np.sqrt(Q.A.env('++', x, 2)*Q.A.env('+-', x, 2))
            out['off'] = (s0['p2']*eo2*G[('C2', 'off')]/m, s0['p3']*u*eo2/(2 - y)*G[('C3H', 'off')]/m)
        for h in hs:
            out[h] = (st[h]['p2']*Q.A.env(h, x, 2)*G[('C2', h)]/m**2, st[h]['p3']*u*Q.A.env(h, x, 3)*G[('C3H', h)]/m**2)
    return out


def lam_of(d, sel=None, scale=1.0):
    hs = d['hs']; n = len(hs)
    sel = np.ones_like(d['nu'], bool) if sel is None else sel
    dF = np.zeros((n, n))
    for i, h in enumerate(hs):
        dF[i, i] = (d[h][0][sel].sum() + d[h][1][sel].sum())*scale
    if n == 2:
        dF[0, 1] = dF[1, 0] = (d['off'][0][sel].sum() + d['off'][1][sel].sum())*scale
    return E.lam2(d['base'], d['Ns'], dF), dF


if __name__ == '__main__':
    Dg = int(sys.argv[1]); tag = sys.argv[2]
    buf = io.StringIO()

    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    pr(f"# ng_chain Dg={Dg} tag={tag} K={KK} JEVEN={Q.JEVEN} JODD={Q.JODD} extra={os.environ.get('Q3_EXTRA_M', '')} "
       f"rows={os.environ.get('W2_ROWS', 'all')} noown={NOOWN}")
    grids, out = Q.run([Dg], pr)
    pr(f"# despiked: {W.DESPIKE_LOG}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(NGC, 'qb_chain_%s.txt' % tag), 'w').write(buf.getvalue())
    pickle.dump(dict(grids=grids, recs=W.perstate(out[Dg]['recs']), summary={k: v for k, v in out[Dg].items() if k != 'recs'}),
                open(os.path.join(NGC, 'qb_chain_%s.pkl' % tag), 'wb'))
