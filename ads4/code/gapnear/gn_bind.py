"""
gn_bind (GAP_NEARGAP.md): the error-budgeted near-gap states of Omega^EE at Delta_g = 2000 / 5000 on their own-mass K2 rows
(gn_rows.py), with NEARGAP2's budget run UNCHANGED (neargap2/ng_bind.run_state, read-only import), side by side with the same
protocol at Delta_g = 1000 (NEARGAP2's own rows, the published / ee1000 rows; the gapnear rows hidden), and the H_mono-A test
(Gamma(x) in c2 at fixed (J, x) on the rows 1000 + J, 2000 + J, 5000 + J).
In-process redirections only (nothing in neargap2/, spinfunc/, gapscale/ is written; gn_guard refuses it):
  * the row loader w2_chain._pk also reads gapnear/cache/rows and gapscale/cache/rows (unless HIDE[0]: then only the rows NEARGAP2 saw);
  * q3_analysis.odd_grid / even_grid: identical copies with the row cut D <= 3200 raised to D <= 6000 (5000 rows), and with the
    reference-zero rule of HIGHSPIN §2.1 (hs_core.refbad: drop a K2 node whose Delta_ref = 2e4 reference record has excess
    |r2/r0| > 0.1 -- a property of the reference, decided before Gamma is looked at) applied to the gapnear rows (GN_REFBAD=1, default)
    and, for the like-for-like comparison, to the ee1000 rows of odd J >= 17 (the published shape chain dropped that node too);
  * ng_bind.gauge_factor at Delta_g = 2000 / 5000 uses the min-norm weight of gapscale/gs_weights.py (cache/weights_<Dg>.npz, read
    only), exactly as gapscale/gs_bind.py.
Usage: python3 gn_bind.py TAG SPEC [SPEC ...]       SPEC = Dg:J:Down:mfl[,mfl...]   (mfl 'tw' = twist edge of Dg, 'tw1000' = that of 1000)
       python3 gn_bind.py --hmono J[,J...]          (Gamma tables)
Output: cache/gn_bind_<TAG>.txt / .json
"""
import os, sys, glob, pickle, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gn_guard  # noqa: F401   (registers itself as ng_guard)
GNR = os.path.join(HERE, 'cache', 'rows')
GSR = os.path.join(HERE, '..', 'gapscale', 'cache', 'rows')
GSW = os.path.join(HERE, '..', 'gapscale', 'cache')
sys.path.insert(0, os.path.join(HERE, '..', 'neargap2'))
import numpy as np
_WTS = {Dg: dict(np.load(os.path.join(GSW, 'weights_%d.npz' % Dg))) for Dg in (2000, 5000)}
import ng_chain as NC
import ng_bind as NB
W, Q, E = NC.E.W, NC.Q, NC.E
DREF = 20000.0
HIDE = [False]                      # True: only the rows NEARGAP2 / HIGHSPIN saw (Delta_g = 1000 reproduction)
REFBAD = os.environ.get('GN_REFBAD', '1') == '1'
REFBAD_TOL = 0.1
_ngpk = W._pk


def _pk(pats):
    R = _ngpk(pats)
    if not NC.NOOWN and not HIDE[0]:
        for p in pats:
            for d in (GSR, GNR):
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


# ------------------------------------------------------------------------------------------------ reference-zero rule
_RB = {}


def refdiag(J, D):
    """(x, excess) per K2 node of row D (max over the EE/C2 and C3H records): HIGHSPIN hs_core.refdiag, re-typed (reads every row
    file the chain reads)."""
    pats = (['w1_p3_%s_J%d*.pkl', 'ee1000_p3_%s_J%d*.pkl', 'q2_p3_%s_J%d_*.pkl', 'p3_pts_%s_J%d.pkl'] if J % 2 else
            ['w1_pts_%s_J%d_k2*.pkl', 'ee1000_pts_%s_J%d_k2*.pkl', 'l9_pts_%s_J%d_k2.pkl', 'q2_pts_%s_J%d_k2*.pkl', 'u5_pts_%s_J%d_k2*.pkl'])
    m = E.mphys(D, J); ex = {}
    for kind in ('EE', 'C3H'):
        R = {}
        for p in pats:
            for d in (NC.SFC, NC.NGC, GSR, GNR):
                for fn in sorted(glob.glob(os.path.join(d, p % (kind, J)))):
                    R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items()})
        out = []
        for (DD, nu), r in R.items():
            if DD != D or (DREF, nu) not in R:
                continue
            rf = R[(DREF, nu)]
            a = np.abs(np.asarray(rf['M0'])).max()
            out.append((nu/m, nu, float(np.abs(np.asarray(rf['M2'])).max()/a)))
        if not out:
            continue
        c = float(np.median([v*nu for _, nu, v in out]))
        for x, nu, v in out:
            ex[round(x, 6)] = max(ex.get(round(x, 6), -np.inf), v - c/nu)
    return sorted(ex.items())


def refbad(J, D):
    if (J, D) not in _RB:
        _RB[(J, D)] = {x for x, e in refdiag(J, D) if e > REFBAD_TOL}
    return _RB[(J, D)]


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
    return brentq(f, ts[i[0]], ts[i[0] + 1])*E.mphys(DREF, J)/E.mphys(D, J)


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
        if mine:            # a node with a non-finite entry at this Laplace order (K = 1 at a reference zero, row 5012) is dropped too
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


# ------------------------------------------------------------------------------------------------ H_mono-A tables
def hmono(J, buf, xs=(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9)):
    HIDE[0] = False; NC.NOOWN = False
    g, _ = NB.grid_for(J, 2, False)
    g1, _ = NB.grid_for(J, 1, False)
    Ds = [D for D in sorted(g) if D >= 990]
    pr(f"## H_mono-A, J = {J}: Gamma(x) in c2 at K = 2 [K2 - K1 step], rows " +
       ", ".join(f"{D:g} (m={g[D]['m']:.1f}, J^2/m={J*J/g[D]['m']:.3f}, x_refzero={ref_zero_x(J, D):.3f}, dropped {DROPLOG.get((J, D), [])})"
                 for D in Ds), buf)
    keys = ['C2', 'C3H'] if J % 2 else ([('C2', '++'), ('C3H', '++')] if J < 4 else
                                        [('C2', '++'), ('C2', '+-'), ('C3H', '++'), ('C3H', '+-'), ('C2', 'off'), ('C3H', 'off')])
    res = {}
    for key in keys:
        pr(f"   {str(key):16s} x: " + "  ".join(f"{D:>16g}" for D in Ds), buf)
        for x in xs:
            row = []
            for D in Ds:
                r, r1 = g[D], g1.get(D)
                if key not in r or not (r['xs'][0] - 1e-9 <= x <= r['xs'][-1] + 1e-9):
                    row.append(f"{'—':>16s}"); continue
                fac = r['m'] if (isinstance(key, tuple) and key[1] == 'off') else 1.0
                v = np.interp(x, r['xs'], r[key])*fac
                s = (abs(np.interp(x, r1['xs'], r1[key])*fac - v) if (r1 is not None and key in r1) else np.nan)
                res.setdefault(str(key), {}).setdefault(f"{D:g}", []).append((x, float(v), float(s)))
                row.append(f"{v:+9.1f} [{s:5.0f}]")
            pr(f"   {'':16s} {x:.2f} " + "  ".join(row), buf)
    return res


# ------------------------------------------------------------------------------------------------ states
def mfl_of(tok, Dg, J):
    if tok == 'tw':
        return float(E.mgap(Dg, J))
    if tok == 'tw1000':
        return float(E.mgap(1000, J))
    return float(tok)


def run(specs, buf):
    out = []
    for spec in specs:
        Dg, J, Down, mf = spec.split(':')
        Dg, J, Down = int(Dg), int(J), float(Down)
        HIDE[0] = (Dg == 1000)
        pr(f"\n# Delta_g = {Dg}, J = {J}: own row Delta = {Down:g}; twist edge m/Lam = {E.mgap(Dg, J):.7f}; "
           f"rows hidden: {HIDE[0]}; Lam = {E.lam_of(Dg):.3f}", buf)
        for tok in mf.split(','):
            mfl = mfl_of(tok, Dg, J)
            r = NB.run_state(Dg, J, Down, mfl, lambda s: pr(s, buf))
            r.update(Dg=Dg, tok=tok, dropped={f"{k[1]:g}": v for k, v in DROPLOG.items() if k[0] == J})
            out.append(r)
    HIDE[0] = False
    return out


if __name__ == '__main__':
    t0 = time.time()
    buf = []
    if sys.argv[1] == '--hmono':
        res = {}
        for J in [int(v) for v in sys.argv[2].split(',')]:
            res[J] = hmono(J, buf)
        tag = 'hmono_' + sys.argv[2].replace(',', '_')
        json.dump({str(k): v for k, v in res.items()}, open(os.path.join(HERE, 'cache', 'gn_%s.json' % tag), 'w'))
    else:
        tag = sys.argv[1]
        res = run(sys.argv[2:], buf)
        json.dump(res, open(os.path.join(HERE, 'cache', 'gn_bind_%s.json' % tag), 'w'), default=float)
    pr(f"[{time.time() - t0:.0f}s]", buf)
    open(os.path.join(HERE, 'cache', 'gn_%s.txt' % (tag if tag.startswith('hmono') else 'bind_' + tag)), 'w').write("\n".join(buf) + "\n")
