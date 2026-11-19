"""
w2 (EE_NEARGAP.md): the S1_S4 near-gap chain (q3_analysis: K2 linear chain, aligned off-diagonals, c_off = m Gamma_off model, combined
near-gap states, C^_3 without de-shift) with the C^_2 crystal rows of EE (u5 + w1 caches) in place of q^2's, for ALL even and odd J that
have EE rows; C^_3 rows (C3H, o_H q) unchanged (q2/l9/p3 caches + w1 rows for the new J).  Generalises u6_ee_neargap:
  * even J: any J (l8 q^2 records are NOT mixed into the EE C2 rows: Recs patched); odd J: p3_odd_action rows of EE (w1_p3_EE_*);
  * helicity basis per row (W2_BASIS): 'q2' (u6 convention: rank-1 basis of the q^2 order-0 matrix, where q^2 records exist on the row,
    else EE's own) or 'ee' (EE's own order-0 matrix on the row and at the reference, w1 'refb' points);
  * local k = 2 de-shift (1 +- u)/2 with u = u_EE (exact rational, u2) in EE mode, u_T in q^2 mode (both sectors);
  * every record: logK -> logK - i pi (Delta + J) (q3._strip; asserted);
  * W2_MODES = 'ee' (default) or 'q2,ee' (both chains, side by side where both exist).
Per-state records (base, N_s, dF) are pickled for the join with the band model (w4).
Usage: Q3_JEVEN=... Q3_JODD=... [Q3_EXTRA_M=1.4,1.5,...] W2_BASIS=q2 W2_TAG=x python3 w2_chain.py 300 1000
       -> cache/w2_chain_<tag>.txt, cache/w2_recs_<tag>.pkl
"""
import os, sys, io, glob, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import q3_analysis as Q
import u2_hclosed_struct as U2

C = os.path.join(HERE, 'cache')
BASIS = os.environ.get('W2_BASIS', 'q2')
MODES = os.environ.get('W2_MODES', 'ee').split(',')
TAG = os.environ.get('W2_TAG', '')
MODE = {'ee': False}
_F = U2.numfuncs()
_uq = Q.u_spin
_uE = lambda v: (lambda a, b: (a - b)/(a + b))(_F['EE'][0](np.asarray(v, complex)), _F['EE'][1](np.asarray(v, complex)))


def _pk(pats):
    R = {}
    for p in pats:
        for fn in sorted(glob.glob(os.path.join(C, p))):
            R.update(pickle.load(open(fn, 'rb')))
    return R


def ee_even(J):
    return _pk(['u5_pts_EE_J%d_k2*.pkl' % J, 'w1_pts_EE_J%d_k2*.pkl' % J])


# ---------------------------------------------------------------- direct (RQMC) STATE values with K2 references (W2_DSTATE=1)
# The n1 contour method fails at the heavy reference for nu <~ 35 at large J (masking; e.g. J = 24, nu = 25/33: direct/K2 = 0.15/0.83),
# while the state points agree with K2 to ~1 %.  So: state record M0 := direct all-orders matrix (TPS normalisation), M2 = M4 = 0;
# reference = K2 (M2/M0 ~ 1e-3 at Delta_ref = 2e4).  Applied to the raw records before the phase strip; the basis keeps the K2 M0.
DSTATE = os.environ.get('W2_DSTATE', '0') == '1'
DS_NUMIN = float(os.environ.get('W2_DS_NUMIN', '30'))
DS_LOG = []


def _dstate(R, J, kind):
    if not DSTATE:
        return R
    import n1_analysis as N1
    for fn in sorted(glob.glob(os.path.join(C, '*shadow_k2*.pkl'))):
        N1.SHK2.update(pickle.load(open(fn, 'rb')))
    fn = os.path.join(C, 'w3_n1_%s_J%d.pkl' % (kind, J))
    if not os.path.exists(fn):
        return R
    R = dict(R)
    I2 = dict(Ti=np.eye(2 if J >= 4 else 1))
    for (D, nu), r in pickle.load(open(fn, 'rb')).items():
        if D >= 1e4 or nu < DS_NUMIN:
            continue
        k = (float(D), round(float(nu), 9))
        if k not in R:
            continue
        H, l = N1.Hrep(r, I2, D, J)
        st = R[k]
        M = np.exp(l - (st['logC'] - st['logK']))*H.mean(0)
        R[k] = dict(st, M0=M, M2=0*M, M4=0*M)
        DS_LOG.append((kind, J, D, round(nu, 2)))
    return R


# ---------------------------------------------------------------- even loaders
_q3_load = Q._load                       # q3's loader (J >= 11: q2 caches; J < 11: l9 caches), with phase strip
_q2_even_orig = Q._q2_even


def _q2_even(kind, J):
    R = _q2_even_orig(kind, J)
    R.update(_pk(['w1_pts_%s_J%d_k2*.pkl' % (kind, J)]))
    return R


Q._q2_even = _q2_even


def _load(kind, J, sfx):
    if MODE['ee'] and kind == 'C2':
        return Q._strip(_dstate(ee_even(J), J, 'EE'), J) if sfx == '_k2' else {}
    if MODE['ee'] and kind == 'C3H' and DSTATE and sfx == '_k2':
        R = _q3_load(kind, J, sfx) if J < 11 else {}
        R = dict(R); R.update(Q._strip(_dstate(_q2_even_orig(kind, J) | _pk(['w1_pts_%s_J%d_k2*.pkl' % (kind, J)]), J, 'C3H'), J))
        return R
    if J < 11 and sfx == '_k2':
        R = _q3_load(kind, J, sfx)
        R.update(Q._strip(_pk(['w1_pts_%s_J%d_k2*.pkl' % (kind, J)]), J))
        return R
    return _q3_load(kind, J, sfx)


Q.A._load = _load
_Recs_init = Q.A.Recs.__init__


def _recs_init(self, kind, J, K):
    _Recs_init(self, kind, J, K)
    if MODE['ee'] and kind == 'C2':                  # drop the l8 q^2 records mixed in by Recs (K = 2, kind C2)
        self.R = dict(_load('C2', J, '_k2'))


Q.A.Recs.__init__ = _recs_init
_orig_basis = Q.A.L.basis
BASIS_LOG = {}


def _basis(J, D, R, *a, **k):
    Rq = dict(Q.A.L.records(J)); Rq.update(_q3_load('C2', J, '_k2')); Rq.update(Q._strip(_pk(['w1_pts_C2_J%d_k2*.pkl' % J]), J))
    Re = Q._strip(ee_even(J), J)
    src = 'q2'
    if BASIS == 'ee' or (MODE['ee'] and _orig_basis(J, D, Rq, *a, **k) is None):
        src = 'ee'
    B = _orig_basis(J, D, Rq if src == 'q2' else Re, *a, **k)
    if B is None and src == 'q2':
        src = 'ee'; B = _orig_basis(J, D, Re, *a, **k)
    BASIS_LOG[(MODE['ee'], J, D)] = src
    return B


Q.A.L.basis = _basis


# ---------------------------------------------------------------- odd loaders
def _odd(kind, J, anti=False, k3=False):
    if MODE['ee'] and kind == 'C2':
        R = _pk(['w1_p3_EE_J%d*.pkl' % J])
    else:
        R = Q._orig_p4load_raw(kind, J, k3) if J < 11 else {}
        R.update(_pk(['q2_p3_%s_J%d_*.pkl' % (kind, J), 'w1_p3_%s_J%d*.pkl' % (kind, J)]))
    R = {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}
    return Q._strip(R, J)


Q.P4.load = _odd


# ---------------------------------------------------------------- hybrid states for every J with grids (W2_HYB_ALL=1)
if os.environ.get('W2_HYB_ALL', '1') == '1':
    import inspect, textwrap
    _src = textwrap.dedent(inspect.getsource(Q.neargap))
    _old = "            if J < 11 or (mfl <= 1.3 + 1e-9 and J <= 20):          # J <= 20: exact-skeleton states exist up to 1.3 Lam\n"
    assert _old in _src
    _src = _src.replace(_old, "            if (round(mfl, 9), J) in _exist:\n")
    _src = _src.replace("    states = list(load_states(Dg))\n",
                        "    states = list(load_states(Dg))\n    _exist = set((round(s_[2]['mfl'], 9), s_[2]['J']) for s_ in states)\n")
    # gauge of the Ĉ'_EE weight at small nu (W2_GAUGE): the chain's local form uses psi2 node weights; 'min' / 'can' rescale the p2
    # weights of every beyond-leading-Regge term (graviton correction dF and scalar remainder rem; NOT the flat B2 nor the leading-Regge
    # curvature cv2, which are gauge-invariant by Lemma N) by Re X_gauge/psi2 at nu = p Lam <= W2_GNU (u3 X_min, x2 X_can).
    _o1 = "        rem = (p2*envk(2, x)*fac*r0['rho'][2]).sum()"
    _o2 = "base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2, p3=p3,"
    assert _o1 in _src and _o2 in _src
    _src = _src.replace(_o1, "        rem = (p2*_GW(Dg, Lm*pn)*envk(2, x)*fac*r0['rho'][2]).sum()")
    _src = _src.replace(_o2, "base=B2 + B3 + cv2 + cv3 + rem, Ns=Ns, r0=r0, x=x, y=y, p2=p2*_GW(Dg, Lm*pn), p3=p3,")
    exec(compile(_src, 'w2_neargap', 'exec'), Q.__dict__)

GAUGE = os.environ.get('W2_GAUGE', 'none')
GNU = float(os.environ.get('W2_GNU', '25'))
_GWC = {}


def _gw(Dg, nu):
    nu = np.asarray(nu, float)
    if GAUGE == 'none':
        return np.ones_like(nu)
    if Dg not in _GWC:
        a = pickle.load(open(os.path.join(C, 'x2_chi_EE_can_%d.pkl' % Dg), 'rb'))
        b = pickle.load(open(os.path.join(C, 'u3_chi_EE_%d.pkl' % Dg), 'rb'))
        phi = np.interp(a['nu'], b['nu'], b['phi'].real)
        X = a['X'] if GAUGE == 'can' else a['X_minnorm']
        _GWC[Dg] = (a['nu'], X.real/phi)
    g, r = _GWC[Dg]
    return np.where(nu <= GNU, np.interp(nu, g, r), 1.0)


Q._GW = _gw


_rf_orig = Q.row_filter
ONLYROW = [float(v) for v in os.environ.get('W2_ROWS', '').split(',') if v]      # keep only these Delta rows (diagnostic)
DMIN = float(os.environ.get('W2_DMIN', '0'))


def _row_filter(J, Ds):
    Ds = [D for D in _rf_orig(J, Ds) if D >= DMIN]
    return [D for D in Ds if D in ONLYROW] if ONLYROW else Ds


Q.row_filter = _row_filter


if os.environ.get('W2_INTERP', 'spline') == 'linear':      # x-interpolation of the row grids: natural cubic spline (q3) or linear
    class _Lin:
        def __init__(self, x, y, bc_type=None):
            self.x, self.y = np.asarray(x), np.asarray(y)
        def __call__(self, t):
            t = np.asarray(t, float)
            s = (self.y[1] - self.y[0])/(self.x[1] - self.x[0])
            return np.where(t < self.x[0], self.y[0] + s*(t - self.x[0]), np.interp(t, self.x, self.y))
    Q.NG.CubicSpline = _Lin; Q.CubicSpline = _Lin

XLOW = os.environ.get('W2_XLOW', 'spline')         # crystal corrections below a row's first x node: spline (q3) | hold | zero
_ng_model = Q.NG.model


def _model(rows, J, kinds):
    f, Ds = _ng_model(rows, J, kinds)
    if XLOW == 'spline':
        return f, Ds
    x0 = max(rows[D]['xs'][0] for D in rows)
    def g(x, m):
        x = np.asarray(x, float)
        out = f(np.maximum(x, x0), m)
        if XLOW == 'zero':
            out = {k: np.where(x < x0, 0.0, v) for k, v in out.items()}
        return out
    return g, Ds


Q.NG.model = _model


DIRECT = os.environ.get('W2_DIRECT', '0') == '1'     # EE mode: replace (J, Delta) rows by w3's direct-RQMC rows where computed
_A_grid = Q.A.grid


def _grid(J, K, kinds=('C2', 'C3H'), *a, **k):
    g = _A_grid(J, K, kinds, *a, **k)
    if os.environ.get('W2_XSUB', '') == '19':          # x-resolution test: keep only the 19 standard nodes of each row
        X19 = np.array([0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8])
        for D, r in g.items():
            sel = np.array([np.min(np.abs(X19 - x)) < 1e-4 for x in r['xs']])
            for k, v in list(r.items()):
                if isinstance(v, np.ndarray) and v.shape == sel.shape:
                    r[k] = v[sel]
    for spec in [v for v in os.environ.get('W2_DROPX', '').split(',') if v]:      # "J:Delta:x1/x2/..." drop nodes (outlier tests)
        Js, Dsx, xl = spec.split(':')
        if int(Js) == J and float(Dsx) in g:
            r = g[float(Dsx)]
            sel = np.array([all(abs(x - float(v)) > 1e-4 for v in xl.split('/')) for x in r['xs']])
            for k, v in list(r.items()):
                if isinstance(v, np.ndarray) and v.shape == sel.shape:
                    r[k] = v[sel]
    if os.environ.get('W2_DESPIKE', '1') == '1':    # Hampel filter: drop isolated K2 outliers (|G - local median| > 6 x robust scale)
        for D, r in g.items():
            n = len(r['xs']); bad = np.zeros(n, bool)
            for k in [kk for kk in r if isinstance(kk, tuple) and len(kk) == 2 and kk[1] in ('++', '+-')]:
                v = np.asarray(r[k], float)
                scale = np.median(np.abs(v)) + 1e-12
                for i in range(n):
                    w = np.r_[v[max(0, i - 3):i], v[i + 1:i + 4]]
                    med = np.median(w); mad = np.median(np.abs(w - med))
                    if abs(v[i] - med) > 6*max(mad, scale):
                        bad[i] = True
            bad[0] = False                                      # never drop the first node (small-x end, no left neighbours)
            if bad.any():
                DESPIKE_LOG.append((J, D, [round(float(x), 4) for x in r['xs'][bad]]))
                for k, v in list(r.items()):
                    if isinstance(v, np.ndarray) and v.shape == bad.shape:
                        r[k] = v[~bad]
    if DIRECT and MODE['ee'] and K == 2:
        fn = os.path.join(C, 'w3_direct_grids.pkl')
        dg = pickle.load(open(fn, 'rb')) if os.path.exists(fn) else {}
        for (Jr, Dr), v in dg.items():
            if Jr == J:
                g[Dr] = v['direct']; DIRECT_LOG.append((J, Dr))
    return g


DIRECT_LOG = []
DESPIKE_LOG = []
Q.A.grid = _grid


def set_mode(ee):
    MODE['ee'] = ee
    Q.u_spin = _uE if ee else _uq                    # q3._shift_generic reads the module-level name u_spin at call time
    Q.ALIGN_LOG.clear()


def perstate(recs):
    """per state: lambda/N_s (corrected, uncorrected), diag/N_h (signed, per helicity), off/sqrt(NN)."""
    out = []
    for r in recs:
        dg = [r['dF'][k, k]/r['Ns'][k] for k in range(len(r['hs']))]
        of = abs(r['dF'][0, 1])/np.sqrt(r['Ns'].prod()) if len(r['hs']) == 2 else 0.0
        out.append(dict(mfl=r['mfl'], J=r['J'], DF=r['DF'], hs=r['hs'], lam=Q.lam(r, 1.0), lam0=Q.lam(r, 0.0), diag=dg, off=of,
                        cov=r['cov'], base=r['base'], Ns=r['Ns'], dF=r['dF']))
    return out


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    Dgs = [int(a) for a in sys.argv[1:]] or [300, 1000]
    t0 = time.time()
    res, store = {}, {}
    for md in MODES:
        ee = md == 'ee'
        set_mode(ee)
        pr(f"\n######## C^_2 rows: {'EE' if ee else 'q^2'};  basis mode {BASIS}; JEVEN {Q.JEVEN}; JODD {Q.JODD}; extra m {os.environ.get('Q3_EXTRA_M', '')}")
        grids, out = Q.run(Dgs, pr)
        pr("# basis source per (J, Delta): " + "; ".join(f"{J},{D:.0f}:{s}" for (e, J, D), s in sorted(BASIS_LOG.items()) if e == ee))
        res[md] = out
        store[md] = dict(grids=grids, recs={Dg: perstate(out[Dg]['recs']) for Dg in Dgs},
                         summary={Dg: {k: v for k, v in out[Dg].items() if k != 'recs'} for Dg in Dgs})
    if len(MODES) == 2:
        pr("\n# side by side (per J: min lam/N_s | max |dF_+-|/sqrt(NN) | max |dF_hh|/N_h):")
        for Dg in Dgs:
            for J in sorted(res['ee'][Dg]['perJ']):
                a = res['q2'][Dg]['perJ'].get(J); b = res['ee'][Dg]['perJ'].get(J)
                if a and b:
                    pr(f"   Delta_g = {Dg}, J = {J:2d}: q^2 {a['min']:+.4f} | {a['off']:.4f} | {a['diag']:.4f}   EE {b['min']:+.4f} | "
                       f"{b['off']:.4f} | {b['diag']:.4f}   (uncorrected {a['unc']:+.4f})")
    pr(f"# despiked nodes (J, Delta, x): {DESPIKE_LOG}")
    pr(f"# direct state points used (W2_DSTATE): {len(DS_LOG)}: {sorted(set(DS_LOG))[:40]}")
    pr(f"# phase strip: {len(Q.PHASE_LOG)} records, max dev {max(Q.PHASE_LOG) if Q.PHASE_LOG else 0:.1e}; non-finite dropped: {len(Q.NONFINITE)}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'w2_chain_%s.txt' % TAG), 'w').write(buf.getvalue())
    pickle.dump(store, open(os.path.join(C, 'w2_recs_%s.pkl' % TAG), 'wb'))
