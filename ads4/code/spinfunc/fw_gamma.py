"""
fw_gamma (FOLD_WINDOW_LARGEJ.md): graviton Gamma_{k,h} (k = 2: C^'_EE, k = 3: C3H; h = ++, +-) on a large-J near-gap row from the fw_points
DIRECT nodes (state + Delta_ref = 2e4 on common random numbers), in exactly the l9/q3/w2 EE-mode conventions (EE's own helicity basis
from order-0 K2 records -- fw_rows, fixed engine -- at the state (x >= 0.15) and at the reference (refb), logK - i pi (Delta + J), spinning
local de-shift with u_EE for k = 2 (none for k = 3), Dc_fixed (ee1000_common Jacobi fix), per-row a_h fitted on x <= 0.25, nu >= 15),
obtained by running q3's grid on the injected records (w2 Hampel despike OFF).  Then the smooth-part extraction of SMOOTH_BOUND §1.3 /
sm_foldfit:  y = Gamma - kappa-image (sm_curvfix.curv_image, leading-Regge curvature in partial-wave phase), least squares on the fold
window NUMIN <= nu <= NUMAX with the oscillating class {D_h, x D_h'/(2J+1)} (+ x^2 D_h) plus a constant r (and variant r0 + r1 nu);
errors: replica jackknife (records recomputed leaving one common-seed replica out) and the residual-scaled lstsq error.
Estimator FW_EST = dir (plain CRN ratio) | cv (scalar control variate with the exact scalar ratio from sm_exact, fw_scalar.install).
Usage: python3 fw_gamma.py J:Delta[,J:Delta...] [NUMIN NUMAX]    -> cache/fw_gamma.txt (appended), cache/fw_gamma.pkl
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ['W2_DESPIKE'] = '0'
os.environ['SM_CURV'] = 'none'
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import ee1000_common as E          # EE mode, Jacobi fix, w2/q3 loaders
import sm_curvfix as CF
import z3_nodes as Z3
import fw_fix as FX

Q, W, C = E.Q, E.W, E.C
DREF = 20000.0
EST = os.environ.get('FW_EST', 'dir')
CUR = {}                            # (J) -> dict(D=..., jk=None)


def fwpts(kind, J, D=None):
    """direct points; legacy file fw_pts_<kind>_J<J>.pkl (state and reference keyed (Delta, nu): a reference at a nu shared by two
    states of the same J is the LAST one computed -- §2.6), then the per-state file fw_pts_<kind>_J<J>_D<D>.pkl (overrides)."""
    R = {}
    f = os.path.join(C, 'fw_pts_%s_J%d.pkl' % (kind, J))
    if os.path.exists(f):
        R.update(pickle.load(open(f, 'rb')))
    D = D if D is not None else CUR.get(J, {}).get('D')
    if D is not None:
        g = os.path.join(C, 'fw_pts_%s_J%d_D%g.pkl' % (kind, J, D))
        if os.path.exists(g):
            R.update(pickle.load(open(g, 'rb')))
    return R


FIXLOG = []


def fwk2(kind, J):
    f = os.path.join(C, 'fw_k2_%s_J%d.pkl' % (kind, J))
    return FX.fix_k2(pickle.load(open(f, 'rb')), J, FIXLOG) if os.path.exists(f) else {}


_z3shadow = Z3.shadow


def _shadow(D, J):                   # the fw (patched-basis) shadows take precedence for the fw rows
    try:
        return FX.bcur(D, J)
    except KeyError:
        return _z3shadow(D, J)


Z3.shadow = _shadow


def direct_records(kind, J, jk=None):
    """z3.to_record of every fw node that has state and reference (dir estimator; jk = replica left out)."""
    P = fwpts(kind, J)
    by = {}
    for (D, nu), r in P.items():
        by.setdefault(nu, {})[D] = FX.fix_direct(r, J, D, nu)
    out = {}
    for nu, dd in by.items():
        if DREF not in dd or len(dd) < 2:
            continue
        R = len(dd[DREF]['reps'])
        sel = None if jk is None else (np.arange(R) != jk)
        for D, r in dd.items():
            fac = 1.0
            if EST == 'cv' and D != DREF:      # scalar control variate: MC scalar ratio replaced by the exact one (fw_scalar objects)
                import fw_scalar as FS
                if (round(float(D), 9), int(J)) not in getattr(Z3, '_exact_override', {}):
                    FS.install(J, D)
                o = Z3.scalar_check(J, nu, dd, D)
                fac = o['cv'] if jk is None else o['cv_jk'][jk]
            out[(float(D), round(nu, 9))] = Z3.to_record(r, D, J, fac, sel)
    return out


_ee0 = W.ee_even
_q20 = Q._q2_even


def ee_even(J):
    if J not in CUR:
        return _ee0(J)
    R = dict(fwk2('EE', J))
    if os.environ.get('FW_NODIRECT', '0') != '1':
        R.update(direct_records('EE', J, CUR[J].get('jk')))
    return R


def q2_even(kind, J):
    if J not in CUR or kind != 'C3H':
        return _q20(kind, J)
    R = dict(fwk2('C3H', J))
    if os.environ.get('FW_NODIRECT', '0') != '1':
        R.update(direct_records('C3H', J, CUR[J].get('jk')))
    return R


W.ee_even = ee_even
Q._q2_even = q2_even


# ------------------------------------------------------------------------------------------------ generalised helicity basis
# At the Delta_g = 300 near-gap states (J^2/Delta >~ 10) the K2-order shadow matrix b is not converged (orders 1 -> 2 change b by 3 % at
# (400, 100)), so the order-0 matrices M0 = S b^-T are NOT symmetric (asym 0.7) and l8's symmetric rank-1 basis fails.  A wrong b is
# a fixed linear map on one index, so the x-independent two-helicity structure survives as M0(x) = A diag(d(x)) B^T with A != B:
# the 4-component span of {M0(x)} is 2-dimensional (sv3/sv1 = 4e-16 on x in [0.15, 0.4] at (400, 100)).  gbasis finds the two rank-1
# matrices a_h b_h^T in that span (det(A1 + t A2) = 0), Ti = A^-1 (left), TiR = B^-1 (right); d_h made real on the basis records by a
# constant phase per h (absorbed, like the modulus, in the per-row a_h); ++ identified by flatness as l8.basis.
GB_X = dict(state=(float(os.environ.get('FW_XB0', '0.15')), float(os.environ.get('FW_XB1', '0.4'))),
            ref=(float(os.environ.get('FW_XR0', '0.15')), float(os.environ.get('FW_XR1', '0.45'))))
GB_LOG = {}


def gbasis(J, D, R):
    import l7_dc2_analysis as An
    m = Z3.mphys(D, J)
    lo, hi = GB_X['ref' if D >= 1e4 else 'state']
    pts = sorted((k[1]/m, R[k]) for k in R if k[0] == D and lo - 1e-9 <= k[1]/m <= hi + 1e-9 and 'slog' in R[k])
    if len(pts) < 3:
        return None
    V = np.array([p[1]['M0'].ravel()/np.abs(p[1]['M0']).max() for p in pts])
    _, sv, Vh = np.linalg.svd(V)
    A1, A2 = Vh[0].reshape(2, 2), Vh[1].reshape(2, 2)
    a = np.linalg.det(A2); c = np.linalg.det(A1)
    b = A1[0, 0]*A2[1, 1] + A1[1, 1]*A2[0, 0] - A1[0, 1]*A2[1, 0] - A1[1, 0]*A2[0, 1]
    cols, rows_ = [], []
    for t in np.roots([a, b, c]):
        Rm = A1 + t*A2
        j = int(np.argmax(np.linalg.norm(Rm, axis=0))); av = Rm[:, j]
        i = int(np.argmax(np.abs(av))); bv = Rm[i, :]/av[i]
        cols.append(av); rows_.append(bv)
    A = np.array(cols).T; Bm = np.array(rows_).T
    Ti, TiR = np.linalg.inv(A), np.linalg.inv(Bm)
    dd = np.array([np.diag(Ti @ p[1]['M0'] @ TiR.T) for p in pts])          # (n, 2)
    ph = np.exp(-1j*np.angle(np.median(dd.real, axis=0) + 1j*np.median(dd.imag, axis=0)))
    Ti = np.diag(ph) @ Ti
    dd = dd*ph[None, :]
    offr = max(np.abs((Ti @ p[1]['M0'] @ TiR.T)[0, 1]) + np.abs((Ti @ p[1]['M0'] @ TiR.T)[1, 0]) for p in pts)/np.abs(dd).max()
    imag = float(np.max(np.abs(dd.imag)/np.abs(dd)))
    import l7_dc2_analysis as An
    xs = np.array([p[0] for p in pts])
    flat = []
    for j in range(2):
        q = np.array([dd[k, j]*np.exp(p[1]['logC'] - p[1]['logK'] - p[1]['slog'])/p[1]['sS'][0] for k, p in enumerate(pts)])/An.W('++', J, xs)
        flat.append(np.std(np.abs(q/q[0])))
    if flat[0] > flat[1]:
        P_ = np.array([[0, 1], [1, 0]]); Ti = P_ @ Ti; TiR = P_ @ TiR
    out = dict(T=np.linalg.inv(Ti), Ti=Ti, TiR=TiR, sv3=sv[2]/sv[0], sv2=sv[1]/sv[0], imag=imag, offres=offr, n=len(pts), xr=(lo, hi))
    GB_LOG[(J, D)] = {k: v for k, v in out.items() if k in ('sv2', 'sv3', 'imag', 'offres', 'n', 'xr')}
    return out


if (os.environ.get('FW_MINNODES') or os.environ.get('FW_AFIT', 'full') != 'full') and not getattr(W._A_grid, '_fw', False):
    # FW_MINNODES: preliminary runs with fewer than 5 nodes on a row (l9 grid's minimum);
    # FW_AFIT=D: the per-row a_h fitted with D_h alone (l9: D_h, x^2, x^4 on x <= 0.25, nu >= 15 -- ill-conditioned when every node has
    # x <~ 0.03, e.g. (1100, 100), where it returned a_+- off by 10^2).  a_h is a per-row constant absorbed in c0 D_h; r scales by 1/(1+eps).
    import inspect, textwrap
    _gs = textwrap.dedent(inspect.getsource(W._A_grid))
    _a_old = "A = np.array([Dsh(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T"
    assert 'if len(xs) < 5 or' in _gs and _a_old in _gs
    Q.A.__dict__['_MINN'] = int(os.environ.get('FW_MINNODES', '5'))
    Q.A.__dict__['_AFULL'] = os.environ.get('FW_AFIT', 'full') == 'full'
    _gs = _gs.replace('if len(xs) < 5 or', 'if len(xs) < _MINN or').replace(
        _a_old, "A = np.array([Dsh(h, J, xs[sel])] + ([xs[sel]**2, xs[sel]**4] if _AFULL else [])).T")
    exec(compile(_gs, 'l9_grid_fw', 'exec'), Q.A.__dict__)
    Q.A.__dict__['grid']._fw = True
    W._A_grid = Q.A.__dict__['grid']
GBASIS = os.environ.get('FW_BASIS', 'gen') == 'gen'
_basis0 = Q.A.L.basis
_helmat0 = Q.A.helmat


def _basis(J, D, R, *a, **k):
    if J in CUR and GBASIS:
        return gbasis(J, D, Q._strip(ee_even(J), J))
    return _basis0(J, D, R, *a, **k)


def _helmat(r, B, K):
    TR = B.get('TiR', B['Ti'])
    return [B['Ti'] @ r[o] @ TR.T for o in ('M0', 'M2', 'M4', 'M6')[:K + 1]] + [0*r['M0']]*(3 - K)


Q.A.L.basis = _basis
Q.A.helmat = _helmat


def row(J, D, jk=None):
    CUR.clear(); CUR[J] = dict(jk=jk, D=float(D))
    Z3._SH.clear()
    m = Z3.mphys(D, J)
    nus = sorted(set(k[1] for k in fwpts('EE', J) if k[0] == float(D)) & set(k[1] for k in fwpts('C3H', J) if k[0] == float(D)))
    XS = [nu/m for nu in nus if (DREF, nu) in fwpts('EE', J) and (DREF, nu) in fwpts('C3H', J)]
    g = Q.A.grid(J, 2, ('C2', 'C3H'), XS=XS, Ds=[float(D)])     # exact x = nu/m (grid's default rounds x to 1e-6)
    CUR.clear()
    return g.get(float(D))


def dshape(h, J, x):
    D = E.Dc_fixed(h, J, x); hh = 1e-6
    return D, (E.Dc_fixed(h, J, x + hh) - E.Dc_fixed(h, J, x - hh))/(2*hh)*x/(2*J + 1)


def fit(x, y, h, J, basis='cs', w=None):
    """least squares (weights w = 1/sigma^2 if given: then the error is the propagated one, max'ed with the residual-scaled one)."""
    D, Dp = dshape(h, J, x)
    cols = [D, Dp] + ([x*x*D] if 'x2' in basis else []) + [np.ones_like(x)] + ([x - x.mean()] if 'lin' in basis else [])
    A = np.array(cols).T
    sw = np.ones_like(y) if w is None else np.sqrt(w)
    c, *_ = np.linalg.lstsq(A*sw[:, None], y*sw, rcond=None)
    res = y - A @ c
    dof = max(len(y) - A.shape[1], 1)
    Aw = A*sw[:, None]
    chi2 = float(((res*sw)**2).sum())/dof
    cov = np.linalg.inv(Aw.T @ Aw)*(max(chi2, 1.0) if w is not None else (res @ res)/dof)
    ir = 2 + ('x2' in basis)
    return c, np.sqrt(np.diag(cov)), ir, np.sqrt(np.mean(res**2)), np.sqrt(np.mean(y**2))


def analyse(J, D, numin=8.0, numax=34.0, pr=print, bases=('cs', 'cs+x2', 'cs+lin')):
    g0 = row(J, D)
    if g0 is None:
        pr(f"J={J} Delta={D:g}: no row (records missing)"); return None
    R = len(next(iter(fwpts('EE', J).values()))['reps'])
    gj = [row(J, D, jk=i) for i in range(R)]
    m = g0['m']; xs = g0['xs']; nus = xs*m
    sel = (nus >= numin - 1e-9) & (nus <= numax + 1e-9)
    out = dict(J=J, D=D, m=m, nus=nus, a=g0['a'], fits={})
    pr(f"J={J} Delta={D:g} m={m:.2f}: nodes nu = {np.round(nus, 2).tolist()};  a = {g0['a']}")
    pr(f"   basis ({'generalised' if GBASIS else 'l8 symmetric'}): " + "; ".join(f"Delta={k[1]:g}: " + ", ".join(
        f"{kk}={vv:.1e}" if isinstance(vv, float) else f"{kk}={vv}" for kk, vv in v.items()) for k, v in GB_LOG.items() if k[0] == J))
    for kd in ('C2', 'C3H'):
        for h in ('++', '+-'):
            key = (kd, h)
            if key not in g0:
                continue
            img = CF.curv_image(kd, h, J, m, xs, Q.A.env, Q.A.SIG)
            y0 = np.asarray(g0[key], float) - img
            yj = np.array([np.asarray(gg[key], float) - img for gg in gj])
            ej = np.sqrt((R - 1)/R*np.sum((yj - yj.mean(0))**2, axis=0))
            out[(kd, h)] = dict(G=np.asarray(g0[key], float), img=img, y=y0, err=ej, Gim=g0.get((kd, h, 'im')))
            pr(f"   {kd:3s} {h}: Gamma      = {np.round(g0[key]).tolist()}")
            pr(f"           kappa img  = {np.round(img).tolist()}")
            pr(f"           G - img    = {np.round(y0).tolist()}  (+- jk {np.round(ej).tolist()})")
            WF = os.environ.get('FW_WFIT', '0') == '1'
            wts = 1.0/np.maximum(ej[sel], 1.0)**2 if WF else None
            for b in bases:
                if sel.sum() < 3 + ('x2' in b) + ('lin' in b) + 1:
                    continue
                c, e, ir, res, rms = fit(xs[sel], y0[sel], h, J, b, wts)
                cj = np.array([fit(xs[sel], yy[sel], h, J, b, wts)[0][ir] for yy in yj])
                ejk = np.sqrt((R - 1)/R*np.sum((cj - cj.mean())**2))
                out['fits'][(kd, h, b)] = (c[ir], e[ir], ejk, c[0], res, rms)
                pr(f"           fit{'w' if WF else ' '}[{b:7s}] on nu in [{numin}, {numax}] ({sel.sum()} nodes): r = {c[ir]:+7.1f} +- {e[ir]:5.1f} (lsq) "
                   f"+- {ejk:5.1f} (jk);  c0 = {c[0]:+8.1f};  resid {res:.0f} of {rms:.0f}")
    return out


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    lo = float(sys.argv[2]) if len(sys.argv) > 2 else 8.0
    hi = float(sys.argv[3]) if len(sys.argv) > 3 else 34.0
    store = {}
    pk = os.path.join(C, 'fw_gamma.pkl')
    if os.path.exists(pk):
        store = pickle.load(open(pk, 'rb'))
    for spec in sys.argv[1].split(','):
        J, D = int(spec.split(':')[0]), float(spec.split(':')[1])
        pr(f"######## fw_gamma J={J} Delta={D:g} estimator {EST}")
        o = analyse(J, D, lo, hi, pr)
        if o is not None:
            store[(J, D, lo, hi, EST)] = o
    pickle.dump(store, open(pk, 'wb'))
    with open(os.path.join(C, 'fw_gamma.txt'), 'a') as f:
        f.write("\n".join(lines) + "\n")
