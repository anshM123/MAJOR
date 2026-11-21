"""
hc_chain (HYPOTHESES_CLOSE.md): the project chains (sm_chainfix modes: curvature double count removed, SM_CURV=row; ee1000 Jacobi fix)
re-run with the hc_points DIRECT nodes injected into the K2 rows, plus the odd-row Hampel despike of st_odd_despike.py (same rule as
w2_chain._grid: window +-3, |G - median| > 6 max(MAD, median|G|), first node kept), all in-process.  Nothing on disk is modified;
outputs only cache/hc_chain_<mode>_<HC_TAG>.txt/.pkl.
Modes:
  ee1000        : ee1000_chain.run() (hybrid large-J states, shape model; env EE1000_*, W2_*)
  w2 Dg [Dg]    : the EE_NEARGAP main chain (q3/w2 over exact-skeleton + hybrid states; env Q3_*, W2_*), ee1000 loaders + Jacobi fix
  z4 300        : z4_chain (EE_DG300_SMALLNU: direct small-nu nodes + non-local min-norm X below nu = 25; env Z4_*, Q3_*, W2_ROWS),
                  Jacobi fix applied
Injection (record loaders, below every chain model; the helicity basis is computed WITHOUT the injected records):
  HC_INJ   = "J:D[,D];J:D..."  rows whose hc points are injected as row nodes (state Delta = D and its own reference), both kinds;
  HC_IXMIN, HC_IXMAX (default 0, 9): x-range of the injected nodes; K2 nodes within 1e-6 (relative) of an injected nu are replaced;
  HC_Z3    = "J:D[,D];..." (z4 mode): rows whose hc points are served to z4 through z3_nodes.zpts (EE/C3H), in addition to z2's own;
             z3's exact-scalar cache is redirected to cache/hc_exact.pkl;
  HC_EDGE_MODEL = piece (default; ee1000): protocol model on x <= HC_XB (0.8), separate shape fit (degree HC_EK = 0) on the nodes
             x >= HC_XB used beyond HC_XB | global: one shape fit over all nodes;
  HC_JK = i: replica i of the injected hc nodes left out (jackknife);  HC_ODDDESPIKE = 1 (default) | 0;  HC_DESPIKE_INJ = 1 (default: injected nodes are subject to the Hampel rules) | 0 (exempt).
Usage: HC_TAG=x HC_INJ=40:1040 EE1000_JEVEN=40 EE1000_XMAX=0.9001 SM_CURV=row python3 hc_chain.py ee1000
"""
import os, sys, io, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

mode = sys.argv[1]
TAG = os.environ.get('HC_TAG', 'x')
C = os.path.join(HERE, 'cache')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


def _parse(s):
    out = {}
    for item in [v for v in s.split(';') if v]:
        J, Ds = item.split(':')
        out[int(J)] = [float(d) for d in Ds.split(',')]
    return out


INJ = _parse(os.environ.get('HC_INJ', ''))
Z3ROWS = _parse(os.environ.get('HC_Z3', ''))
IXMIN = float(os.environ.get('HC_IXMIN', '0'))
IXMAX = float(os.environ.get('HC_IXMAX', '9'))
ODDDESP = os.environ.get('HC_ODDDESPIKE', '1') == '1'
DESP_INJ = os.environ.get('HC_DESPIKE_INJ', '1') == '1'
HCJK = None if os.environ.get('HC_JK', '') == '' else int(os.environ['HC_JK'])     # replica left out (injected hc nodes)

t0 = time.time()
import w2_chain as W2
import hc_nodes as N
Q = W2.Q
INJLOG = []
NOINJ = [False]
INJ_NU = {}          # (J, D) -> injected state nus (for despike exemption / reporting)


def _inject_dict(R, kind_file, J, strip):
    """R: loader output {(D, nu): rec} (already phase-stripped if strip); returns R with the hc nodes of INJ[J] added/overriding."""
    if NOINJ[0] or J not in INJ:
        return R
    R = dict(R)
    N.set_resid(R, J, kind_file)
    for D in INJ[J]:
        add = N.records(kind_file, J, D, jk=HCJK, xmin=IXMIN, xmax=IXMAX)
        if strip:
            add = Q._strip(add, J)
        for k, v in add.items():
            near = [kk for kk in R if kk[0] == k[0] and abs(kk[1] - k[1]) < 1e-6*k[1]]
            for kk in near:
                R.pop(kk)
            R[k] = v
            if k[0] == D:
                INJ_NU.setdefault((J, D), set()).add(k[1])
                INJLOG.append((kind_file, J, D, round(k[1]/N.mphys(D, J), 4), 'replaced K2' if near else 'new'))
    return R


# ------------------------------------------------------------------------------------------------ odd-row despike (st_odd_despike rule)
DESPLOG = []
_p4grid = Q.P4.grid


def _grid_despiked(J, *a, **k):
    g = _p4grid(J, *a, **k)
    if not ODDDESP:
        return g
    for D, r in g.items():
        n = len(r['xs']); bad = np.zeros(n, bool)
        for kd in ('C2', 'C3H'):
            if kd not in r:
                continue
            v = np.asarray(r[kd], float)
            scale = np.median(np.abs(v)) + 1e-12
            for i in range(n):
                w = np.r_[v[max(0, i - 3):i], v[i + 1:i + 4]]
                med = np.median(w); mad = np.median(np.abs(w - med))
                if abs(v[i] - med) > 6*max(mad, scale):
                    bad[i] = True
        bad[0] = False
        if not DESP_INJ and (J, D) in INJ_NU:
            inj = np.array([any(abs(x*r['m'] - v) < 1e-6*v for v in INJ_NU[(J, D)]) for x in r['xs']])
            bad &= ~inj
        if bad.any():
            DESPLOG.append((J, float(D), [round(float(x), 4) for x in r['xs'][bad]]))
            for kk, vv in list(r.items()):
                if isinstance(vv, np.ndarray) and vv.shape == bad.shape:
                    r[kk] = vv[~bad]
    return g


Q.P4.grid = _grid_despiked

# even-row despike exemption of injected nodes: w2's Hampel acts in W2._grid (= Q.A.grid as installed by w2_chain, AFTER the inner
# l9 grid W2._A_grid); wrap Q.A.grid there: w2's filter off for rows with injected nodes, the same rule applied to the non-injected nodes
if not DESP_INJ:
    _Wgrid = Q.A.grid

    def _grid_ex(J, K, kinds=('C2', 'C3H'), *a, **k):
        if J not in INJ:
            return _Wgrid(J, K, kinds, *a, **k)
        old = os.environ.get('W2_DESPIKE', '1')
        os.environ['W2_DESPIKE'] = '0'
        try:
            g = _Wgrid(J, K, kinds, *a, **k)
        finally:
            os.environ['W2_DESPIKE'] = old
        if old != '1':
            return g
        for D, r in g.items():
            n = len(r['xs']); bad = np.zeros(n, bool)
            inj = np.array([any(abs(x*r['m'] - v) < 1e-5*v for v in INJ_NU.get((J, D), ())) for x in r['xs']])
            for kk in [kk for kk in r if isinstance(kk, tuple) and len(kk) == 2 and kk[1] in ('++', '+-')]:
                v = np.asarray(r[kk], float); scale = np.median(np.abs(v)) + 1e-12
                for i in range(n):
                    w = np.r_[v[max(0, i - 3):i], v[i + 1:i + 4]]
                    med = np.median(w); mad = np.median(np.abs(w - med))
                    if abs(v[i] - med) > 6*max(mad, scale):
                        bad[i] = True
            bad[0] = False; bad &= ~inj
            if bad.any():
                W2.DESPIKE_LOG.append((J, D, [round(float(x), 4) for x in r['xs'][bad]]))
                for kk, vv in list(r.items()):
                    if isinstance(vv, np.ndarray) and vv.shape == bad.shape:
                        r[kk] = vv[~bad]
        return g
    Q.A.grid = _grid_ex

import sm_curvfix as CF          # noqa: E402  (wraps Q.A.grid / Q.P4.grid: curvature image subtracted after despike)

EDGE_MODEL = os.environ.get('HC_EDGE_MODEL', 'piece')     # ee1000 mode: 'piece' | 'global'
XB = float(os.environ.get('HC_XB', '0.8'))
EK = int(os.environ.get('HC_EK', '0'))
EREP = os.environ.get('HC_EREP', 'shape')                # edge representation: shape (degree HC_EK) | linear | spline
EDGEFIT = {}

if mode in ('ee1000', 'w2'):
    import ee1000_chain as CH
    W, E = CH.W, CH.E
    if mode == 'ee1000' and EDGE_MODEL == 'piece':
        # piecewise row model: the protocol's model (shape fit K = EE1000_K, or spline) built on the nodes x <= HC_XB only, used for
        # x <= HC_XB; beyond, a separate shape fit Gamma = c D_h + s x D_h'/(2J+1) + r (coefficients polynomials of degree HC_EK in x,
        # held beyond the last node) on the nodes x >= HC_XB (the direct edge nodes).  Without nodes beyond HC_XB this is the protocol's
        # held extrapolation.  Off-diagonal keys: the protocol's spline (held) throughout.
        _row_models = CH.row_models

        def row_models(J, g, keys):
            gb = {}; ge = {}
            for D, r in g.items():
                sb = r['xs'] <= XB + 1e-9; se = r['xs'] >= XB - 1e-9
                gb[D] = {k: (v[sb] if isinstance(v, np.ndarray) and v.shape == sb.shape else v) for k, v in r.items()}
                ge[D] = {k: (v[se] if isinstance(v, np.ndarray) and v.shape == se.shape else v) for k, v in r.items()}
            mb = _row_models(J, gb, keys)
            out = {}
            for D in mb:
                m_, fb, sfb = mb[D]
                re_ = ge[D]
                dk = [k for k in keys if not (isinstance(k, tuple) and k[1] == 'off')]
                sfe = None
                if re_['xs'].size >= 2 and re_['xs'][-1] > XB + 1e-6:
                    if EREP == 'shape' and len(re_['xs']) >= 3*(EK + 1) + 1:
                        sfe = CH.ShapeFit(J, re_, dk, K=EK)
                    elif EREP in ('linear', 'spline'):
                        from scipy.interpolate import CubicSpline, interp1d
                        xe_ = re_['xs']
                        ip = {k: (interp1d(xe_, re_[k]) if EREP == 'linear' else CubicSpline(xe_, re_[k], bc_type='natural')) for k in dk}
                        sfe = lambda key, x, ip=ip, xe_=xe_: ip[key](np.clip(x, xe_[0], xe_[-1]))
                EDGEFIT[(J, D)] = (sfe, re_['xs'])

                def f(key, x, fb=fb, sfe=sfe):
                    x = np.asarray(x, float)
                    v = fb(key, x)
                    if sfe is None or (isinstance(key, tuple) and key[1] == 'off'):
                        return v
                    return np.where(x > XB, sfe(key, x), v)
                out[D] = (m_, f, sfb)
            return out
        CH.row_models = row_models
elif mode == 'z4':
    import ee1000_common as E    # noqa: F401  (Jacobi fix)
    import z3_nodes as Z3
    Z3.EXPK = os.path.join(C, 'hc_exact.pkl')
    _zpts = Z3.zpts

    def zpts(kind, J):
        P = _zpts(kind, J)
        seen = {}
        for D in Z3ROWS.get(J, []):
            for k, v in N.pts(kind, J, D).items():
                if k[0] == N.DREF:
                    assert k not in seen and k not in P, ('reference key collision', J, D, k)
                    seen[k] = D
                P[k] = v
        return P
    Z3.zpts = zpts
    _z3shadow = Z3.shadow
    CURKIND = ['EE']

    def _shadow(D, J):             # z3 converts with z3.shadow; use the row-consistent shadow (resid match, per kind) for the hc rows
        if J in Z3ROWS:
            kd = CURKIND[0]
            if (float(D), J, kd) not in N.RESID:
                N.set_resid(Z4._ee_orig(J) if kd == 'EE' else Z4._q2e_orig(kd, J), J, kd)
            (b, logK, res), _ = N.shadow(D, J, N.RESID.get((float(D), J, kd)))
            return (b, logK, res)
        return _z3shadow(D, J)
    Z3.shadow = _shadow
    _z3rec = Z3.records

    def _records(kind, J, *a, **k):
        CURKIND[0] = kind
        return _z3rec(kind, J, *a, **k)
    Z3.records = _records
    import z4_chain as Z4
    # z4 injects C3H through q3._q2_even, which the w2 loader uses only for J >= 11; for J < 11 (C3H from the l9 caches) inject here
    _q3l = W2._q3_load

    def _q3_load_inj(kind, J, sfx):
        R = _q3l(kind, J, sfx)
        if kind == 'C3H' and J < 11 and sfx == '_k2' and Z4.DIRECT and J in Z4.JS:
            R = dict(R)
            R.update(Q._strip(Z3.records('C3H', J, Z4.EST, Z4.JK, numin=Z4.NUMIN, numax=Z4.NUMAX), J))
        return R
    W2._q3_load = _q3_load_inj
else:
    raise SystemExit('mode?')

# ---- loader wrappers (after the chain modules installed theirs)
_ee = W2.ee_even


def ee_even(J):
    R = _ee(J)
    N.set_resid(R, J, 'EE')
    return _inject_dict(R, 'EE', J, strip=False)


W2.ee_even = ee_even
_q2e = Q._q2_even


def q2_even(kind, J):
    R = _q2e(kind, J)
    if kind == 'C3H':
        N.set_resid(R, J, 'C3H')
        R = _inject_dict(R, 'C3H', J, strip=False)
    return R


Q._q2_even = q2_even
_p4load = Q.P4.load


def p4load(kind, J, anti=False, k3=False):
    R = _p4load(kind, J, anti, k3)
    return _inject_dict(R, 'EE' if kind == 'C2' else kind, J, strip=True)


Q.P4.load = p4load
_basis = Q.A.L.basis


def basis_noinj(*a, **k):
    NOINJ[0] = True
    try:
        return _basis(*a, **k)
    finally:
        NOINJ[0] = False


Q.A.L.basis = basis_noinj

pr(f"# hc_chain {mode} {sys.argv[2:]}: HC_INJ={os.environ.get('HC_INJ', '')} x in [{IXMIN}, {IXMAX}]; HC_Z3={os.environ.get('HC_Z3', '')}; "
   f"odd despike {ODDDESP}; despike of injected nodes {DESP_INJ}; SM_CURV={CF.MODE}; env " +
   " ".join(f"{k}={v}" for k, v in sorted(os.environ.items()) if k.startswith(('Q3_', 'W2_', 'EE1000', 'Z4_', 'SM_'))))
out = {}
if mode == 'w2':
    Dgs = [int(a) for a in sys.argv[2:]] or [1000]
    grids, res = Q.run(Dgs, pr)
    out = {Dg: W2.perstate(res[Dg]['recs']) for Dg in Dgs}
    out['grids'] = grids
elif mode == 'ee1000':
    recs, info = CH.run(pr)
    w = min(recs, key=lambda r: r['lam'])
    pr(f"# overall (these J): min lam/N_s = {w['lam']:+.4f} at J={w['J']}, m/Lam={w['mfl']:.4f}; states {len(recs)}")
    out = dict(recs=recs, info=info)
elif mode == 'z4':
    W2.set_mode(True)
    Dgs = [int(a) for a in sys.argv[2:]] or [300]
    grids, res = Q.run(Dgs, pr)
    Z4.node_table(pr)
    for s in Z4.LOG:
        pr("# " + s)
    if Z3.BEST_LOG:
        pr("# estimator per node (kind, J, nu, chosen, jk rel err dir, cv): " + "; ".join(str(t) for t in sorted(set(Z3.BEST_LOG))))
    out = {Dg: W2.perstate(res[Dg]['recs']) for Dg in Dgs}
    out['grids'] = grids


def node_report(J):
    """per injected row: Gamma at the injected nodes (direct; curvature image subtracted, as the chain uses it) against the baseline
    model (rows WITHOUT injection, nodes x <= HC_BXMAX (0.8), shape fit K = EE1000_K (ee1000 mode) / natural spline (q3), held beyond
    the last node) evaluated at the same x, and the K2 value where a K2 node existed at that x."""
    from scipy.interpolate import CubicSpline
    odd = J % 2 == 1
    gi = Q.odd_grid(J) if odd else Q.even_grid(J)
    NOINJ[0] = True
    try:
        gb = Q.odd_grid(J) if odd else Q.even_grid(J)
    finally:
        NOINJ[0] = False
    bx = float(os.environ.get('HC_BXMAX', '0.8'))
    keys = ['C2', 'C3H'] if odd else [(kd, h) for kd in ('C2', 'C3H') for h in ('++', '+-')]
    for D in INJ.get(J, []):
        if D not in gi or D not in gb:
            pr(f"   node report J={J} D={D:g}: row missing (inj {D in gi}, base {D in gb})"); continue
        ri, rb = gi[D], gb[D]
        sel = rb['xs'] <= bx + 1e-9
        rbt = {k: (v[sel] if isinstance(v, np.ndarray) and v.shape == sel.shape else v) for k, v in rb.items()}
        if mode == 'ee1000' and len(rbt['xs']) >= CH.SHAPEMIN and CH.MODEL == 'shape':
            sf = CH.ShapeFit(J, rbt, keys)
            model = lambda k, x: float(sf(k, np.array([x]))[0])
            lab = 'shape K=%d' % CH.KDEG
        else:
            sp = {k: CubicSpline(rbt['xs'], rbt[k], bc_type='natural') for k in keys}
            model = lambda k, x: float(sp[k](min(x, rbt['xs'][-1])))
            lab = 'spline'
        injx = sorted(round(v/ri['m'], 6) for v in INJ_NU.get((J, D), ()))
        pr(f"   node report J={J} Delta={D:g} (m={ri['m']:.1f}): Gamma (c2; curvature image subtracted) at injected nodes -- direct | "
           f"baseline model [{lab}, nodes x <= {bx}, held beyond] | K2 (if a K2 node existed there)")
        for x in injx:
            i = int(np.argmin(np.abs(ri['xs'] - x)))
            if abs(ri['xs'][i] - x) > 1e-5:
                pr(f"     x={x:.4f}: dropped from the row (despike / incomplete kinds)"); continue
            jb = np.where(np.abs(rb['xs'] - x) < 1e-5)[0]
            s_ = f"     x={x:.4f}:"
            for k in keys:
                kk = k if isinstance(k, str) else k[0] + k[1]
                s_ += f"  {kk}: {ri[k][i]:+8.0f} | {model(k, x):+8.0f} | " + (f"{rb[k][jb[0]]:+8.0f}" if jb.size else "    --  ")
            pr(s_)


if os.environ.get('HC_NODES', '1') == '1' and mode in ('ee1000', 'w2', 'z4'):
    for J in sorted(INJ):
        try:
            node_report(J)
        except Exception as ex:
            pr(f"   node report J={J}: failed ({ex})")
for (J, D), (sfe, xe) in sorted(EDGEFIT.items()):
    if sfe is not None and hasattr(sfe, 'res'):
        pr(f"# edge fit J={J} D={D:g} (K={EK}) on x = {np.round(xe, 4).tolist()}: residual rms/rms [max] " +
           "; ".join(f"{k if isinstance(k, str) else k[0] + ' ' + k[1]}: {a/b:.3f} [{c:.0f}]" for k, (a, b, c) in sfe.res.items()))
pr(f"# injected nodes (kind, J, D, x, action): {sorted(set(INJLOG))}")
pr(f"# hc_nodes log: {sorted(set(map(str, N.LOG)))[:20]}")
pr(f"# even despiked (J, Delta, x): {W2.DESPIKE_LOG}")
pr(f"# odd despiked (J, Delta, x): {DESPLOG}")
pr(f"# phase strip: {len(Q.PHASE_LOG)} records, max dev {max(Q.PHASE_LOG) if Q.PHASE_LOG else 0:.1e}; non-finite dropped: {len(Q.NONFINITE)}")
pr("# " + CF.summary()[:600])
pr(f"[{time.time() - t0:.0f}s]")
open(os.path.join(C, 'hc_chain_%s_%s.txt' % (mode, TAG)), 'w').write(buf.getvalue())
pickle.dump(out, open(os.path.join(C, 'hc_chain_%s_%s.pkl' % (mode, TAG)), 'wb'))
