"""
oo_chain (HYPOTHESES_CLOSE2.md, H_OO at J >= 12): the even-J parity-odd (OO, X^O_same) block of Omega^EE at Delta_g = 300 in the
VALIDATED small-nu treatment (z4 protocol of EE_DG300_SMALLNU §6 / HYPOTHESES_CLOSE §4), i.e. the hc_oo construction moved from the w2
chain (convention A: K2 + spline below x = 0.05, local de-shift) into the z4 chain:
  (1) OO rows = hc_oo.grid_oo (copied below, unchanged conventions: 3x3 p3 action, OO = entry [-1, -1], heavy reference 2e4 order by order,
      ++ shapes (X^O_same couples like ++ at LO), per-row a_O on x <= 0.25, nu >= 15, local de-shift (1 +- u_EE)/2 for k = 2, none for k = 3,
      curvature image subtracted ('row')), built from the K2 OO records (cache/hc_p3_<kind>_J<J>.pkl of the previous agent + this agent's
      per-state cache/oo_p3_<kind>_J<J>_D<D>.pkl) PLUS the direct small-nu OO records (cache/oo_pts_<kind>_J<J>_D<D>.pkl, oo_points.py:
      n3 contour, 3x3 p8 integrand, state + reference on CRN): M0 = E b_odd^-T (b_odd = the p2_shadow_odd matrix whose stored fit residual
      equals the shadow_resid of the row's K2 OO records), M2 = M4 = 0, logC = logN; a direct node replaces a K2 node at the same nu.
      Estimators (OO_EST): dir | cv (scalar control variate with the exact scalar ratio, z3) | best (per node and kind, the smaller
      jackknife error of the OO element);  OO_JK = i: replica i left out.
  (2) the z4 NON-LOCAL k = 2 term for the OO entry: on every OO row with direct nodes, Gamma_c,OO (complex forward-normalised correction,
      curvature-subtracted, built from the OO row exactly as z4._A_grid builds it from the ++ row) and Gamma_NL = a_chi Gamma_c are put in
      z4_chain.NLROWS[(J, D)] under the '++' key ('+-' := 0), so that the chain's ++ density below Z4_NUNL = 25 is Re[(Y/phi) Gamma_NL,OO]
      (min-norm X); the local OO values at direct nodes nu < Z4_NUFLAT are flattened as z4 does; Z4_EXTRAP / Z4_CONST act on the OO rows.
      OO_DUPNL = 1 (default): a single-row J's duplicate row (q3 holds a single row everywhere) carries the same NL row (z4 leaves it local).
  (3) w2's Hampel despike on the OO rows (OO_DESPIKE = 1, default; same rule, first node kept; as the EE rows in z4) or none (= hc_oo).
Two chain runs in one process: run0 = the unmodified z4 chain (EE block; direct EE hc nodes of OO_Z3EE injected exactly as hc_chain z4
HC_Z3), run1 = the even rows of J in OO_J REPLACED by the OO rows (++ keys = Gamma_OO; +- and off-diagonal keys 0).  Per state:
  lambda_OO/N_s = (base_++ + dF_OO)/N_s,++   (1x1; N_same = tan^2(pi tau/2) >= 0 drops out of the sign), vs the EE ++ diagonal and the EE
  2x2 minimum of run0.
Usage: OO_J=18,20,22 OO_TAG=x [OO_EST=best OO_DIRECT=1 OO_JK= OO_DESPIKE=1 Z4_EXTRAP=hold Z4_NL=1 ...] python3 oo_chain.py
       -> cache/oo_chain_<tag>.txt / .pkl
"""
import os, sys, io, glob, time, pickle
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
os.environ.setdefault('SM_CURV', 'row')
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
OOJ = [int(v) for v in os.environ.get('OO_J', '18').split(',') if v]
TAG = os.environ.get('OO_TAG', 'x')
EST = os.environ.get('OO_EST', 'best')
OJK = None if os.environ.get('OO_JK', '') == '' else int(os.environ['OO_JK'])
ODIRECT = os.environ.get('OO_DIRECT', '1') == '1'
ODESP = os.environ.get('OO_DESPIKE', '1') == '1'
DUPNL = os.environ.get('OO_DUPNL', '1') == '1'
ONUMIN = float(os.environ.get('OO_NUMIN', '0'))
ODROWS = [float(v) for v in os.environ.get('OO_DROWS', '').split(',') if v]     # restrict direct OO rows (default: all found)
DGS = [int(v) for v in os.environ.get('OO_DG', '300').split(',')]
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


def _parse(s):
    out = {}
    for item in [v for v in s.split(';') if v]:
        J, Ds = item.split(':')
        out[int(J)] = [float(d) for d in Ds.split(',')]
    return out


# ------------------------------------------------------------------------------------------------ z4 setup (as hc_chain.py z4 mode)
import w2_chain as W2            # noqa: E402
import hc_nodes as N             # noqa: E402
Q = W2.Q
import sm_curvfix as CF          # noqa: E402
import ee1000_common as E        # noqa: E402,F401  (Jacobi fix)
import z3_nodes as Z3            # noqa: E402
Z3.EXPK = os.path.join(C, 'hc_exact.pkl')          # read: the hc exact scalar records (EE z4 rows); OO uses its own file (below)
Z3EE = _parse(os.environ.get('OO_Z3EE', '18:336,400;20:400;22:336'))
_zpts = Z3.zpts


def zpts(kind, J):
    P = _zpts(kind, J)
    for D in Z3EE.get(J, []):
        for k, v in N.pts(kind, J, D).items():
            if k[0] == N.DREF:
                assert k not in P, ('reference key collision', J, D, k)
            P[k] = v
    return P


Z3.zpts = zpts
_z3shadow = Z3.shadow
CURKIND = ['EE']


def _shadow(D, J):
    if J in Z3EE:
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
os.environ.setdefault('Z4_JS', ','.join(str(J) for J in sorted(Z3EE)))
import z4_chain as Z4            # noqa: E402
_q3l = W2._q3_load


def _q3_load_inj(kind, J, sfx):
    R = _q3l(kind, J, sfx)
    if kind == 'C3H' and J < 11 and sfx == '_k2' and Z4.DIRECT and J in Z4.JS and J in Z3EE:
        R = dict(R)
        R.update(Q._strip(Z3.records('C3H', J, Z4.EST, Z4.JK, numin=Z4.NUMIN, numax=Z4.NUMAX), J))
    return R


W2._q3_load = _q3_load_inj
P4 = Q.P4


# ------------------------------------------------------------------------------------------------ OO records
def _merge_files(fs):
    R = {}; owner = {}
    for f in fs:
        if not os.path.exists(f):
            continue
        d = pickle.load(open(f, 'rb'))
        for k, v in d.items():
            if k[0] == DREF:
                assert k not in owner or owner[k] == f, ('reference key shared between files', k, owner.get(k), f)
                owner[k] = f
            R[k] = v
    return R


def load_oo_k2(kind, J):
    fs = [os.path.join(C, 'hc_p3_%s_J%d.pkl' % (kind, J))] + sorted(glob.glob(os.path.join(C, 'oo_p3_%s_J%d_D*.pkl' % (kind, J))))
    R = _merge_files(fs)
    return {(k[0], k[1]): v for k, v in R.items() if k[2] == 1.0}


_OSC = {}


def odd_shadow(D, J, resid):
    if not _OSC:
        fs = ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
              sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl'))) +
              sorted(glob.glob(os.path.join(C, 'hc_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'oo_p2shadow*.pkl'))))
        for f in fs:
            if os.path.exists(f):
                for k, v in pickle.load(open(f, 'rb')).items():
                    _OSC.setdefault(k, []).append((os.path.basename(f), v))
    c = _OSC.get((float(D), int(J), 6), [])
    rr = np.asarray(resid, float).ravel()
    for f, v in c:
        if np.allclose(np.asarray(v[2], float).ravel(), rr, rtol=1e-10, atol=0):
            return v, f
    raise KeyError(('no odd shadow with matching resid', D, J, [f for f, _ in c]))


OOEXPK = os.path.join(C, 'oo_exact.pkl')
SELLOG = []


def _cvfac(J, nu, dd, D, jk):
    """scalar control variate factor (z3.scalar_check with the OO-own exact cache)."""
    old = Z3.EXPK
    hc = os.path.join(C, 'hc_exact.pkl')
    key = (round(float(D), 9), int(J)); keyr = (round(DREF, 9), int(J))
    db = pickle.load(open(hc, 'rb')) if os.path.exists(hc) else {}
    Z3.EXPK = hc if (key in db and keyr in db) else OOEXPK
    try:
        o = Z3.scalar_check(J, nu, dd, D)
    finally:
        Z3.EXPK = old
    return o['cv'] if jk is None else o['cv_jk'][jk]


def oo_pts(kind, J, D):
    f = os.path.join(C, 'oo_pts_%s_J%d_D%g.pkl' % (kind, J, D))
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def oo_direct(kind, J, D, Rk2, est=EST, jk=OJK):
    """direct OO records {(D, nu): rec, (DREF, nu): rec} of one state row (M0 = E b_odd^-T, 3x3; M2 = M4 = 0)."""
    P = oo_pts(kind, J, D)
    by = {}
    for (DD, nu), r in P.items():
        by.setdefault(nu, {})[DD] = r
    res = {}
    for (DD, nu), r in Rk2.items():
        if 'shadow_resid' in r:
            res.setdefault(DD, r['shadow_resid'])
    out = {}
    for nu, dd in sorted(by.items()):
        if DREF not in dd or D not in dd or nu < ONUMIN:
            continue
        R = len(dd[DREF]['reps'])
        sel = np.ones(R, bool) if jk is None else (np.arange(R) != jk)
        st, rf = dd[D], dd[DREF]
        e_ = est
        if est == 'best':
            f1 = lambda s_: st['reps'][s_].mean(0)[-1, -1]/rf['reps'][s_].mean(0)[-1, -1]
            f2 = lambda s_: ((st['reps'][s_].mean(0)[-1, -1]/st['repsS'][s_].mean())/(rf['reps'][s_].mean(0)[-1, -1]/rf['repsS'][s_].mean()))
            e1 = float(Z3._jk(f1, R)[2]/abs(Z3._jk(f1, R)[0])); e2 = float(Z3._jk(f2, R)[2]/abs(Z3._jk(f2, R)[0]))
            e_ = 'cv' if e2 < e1 else 'dir'
            SELLOG.append((kind, J, D, round(nu, 3), e_, round(e1, 4), round(e2, 4)))
        fac = _cvfac(J, nu, dd, D, jk) if e_ == 'cv' else 1.0
        for DDD, rr, fc in ((D, st, fac), (DREF, rf, 1.0)):
            (b, logK, rs), src = odd_shadow(DDD, J, res[DDD])
            binv = np.linalg.inv(np.asarray(b))
            reps = rr['reps'][sel]
            H = np.array([rep @ binv.T for rep in reps]).mean(0)*fc
            Z = 0*H
            out[(float(DDD), round(float(nu), 9))] = dict(M0=H, M2=Z, M4=Z, M6=Z, logC=complex(rr['logN']), logK=complex(logK),
                                                          direct=True, shsrc=src, m=Z3.mphys(DDD, J), shadow_resid=rs)
    return out


def getrec(R, D, nu):
    k = (float(D), round(float(nu), 9))
    if k in R:
        return R[k]
    c = [kk for kk in R if kk[0] == float(D) and abs(kk[1] - nu) < 1e-7*nu]
    return R[c[0]] if c else None


DIRNU = {}         # (J, D) -> direct state nus


def oo_records(kind, J):
    R = load_oo_k2(kind, J)
    if ODIRECT:
        Ds = sorted(set(k[0] for k in R if k[0] < DREF))
        for D in Ds:
            if ODROWS and D not in ODROWS:
                continue
            add = oo_direct(kind, J, D, R)
            for k, v in add.items():
                for kk in [kk for kk in R if kk[0] == k[0] and abs(kk[1] - k[1]) < 1e-6*k[1]]:
                    R.pop(kk)
                R[k] = v
                if k[0] == D:
                    DIRNU.setdefault((J, D), set()).add(k[1])
    return Q._strip(R, J)


# ------------------------------------------------------------------------------------------------ OO grid (hc_oo.grid_oo, records generalised)
def grid_oo(J, K=2, curv=True):
    h = '++'
    RC = {'C2': oo_records('EE', J), 'C3H': oo_records('C3H', J)}
    Ds = sorted(set(k[0] for k in RC['C2'] if k[0] < DREF))
    out = {}
    for D in Ds:
        m = P4.mphys(D, J); mr = P4.mphys(DREF, J)
        xs = sorted(set(round(k[1]/m, 6) for k in RC['C2'] if k[0] == D))
        xs = [x for x in xs if all(getrec(RC[kd], D, x*m) is not None and getrec(RC[kd], DREF, x*m) is not None for kd in RC)]
        if len(xs) < 5 or min(xs) > 0.1:
            continue
        xs = np.array(xs)
        rows = {}; c0 = None; lo = {}
        for kd in ('C2', 'C3H'):
            k = 2 if kd == 'C2' else 3
            G = []; G0 = []
            for x in xs:
                nu = x*m; st = getrec(RC[kd], D, nu); rf = getrec(RC[kd], DREF, nu); xr = nu/mr
                dl = st['logC'] - st['logK'] - rf['logC'] + rf['logK']
                if c0 is None:
                    c0 = dl.real
                dl = dl - c0
                fk = (m*m/(mr*mr)) if k == 3 else 1.0
                ords = ('M0', 'M2', 'M4')
                v = np.exp(dl)*P4.pred([st[o][-1, -1] for o in ords], [rf[o][-1, -1] for o in ords], K)*fk
                v0 = np.exp(dl)*st['M0'][-1, -1]/rf['M0'][-1, -1]*fk
                G.append(v*P4.env(h, xr, k)*P4.Dc(h, J, xr)); G0.append(v0)
            rows[kd] = (k, np.array(G))
            lo[kd] = np.array(G0)
        sel = (xs <= 0.25 + 1e-9) & (xs*m >= 15.0)
        if sel.sum() < 3:
            sel = np.argsort(xs)[:3]

        def fit(g):
            A_ = np.array([P4.Dc(h, J, xs[sel]), xs[sel]**2, xs[sel]**4]).T
            return np.linalg.lstsq(A_, (g/P4.env(h, xs, 2))[sel], rcond=None)[0][0]
        a = fit(rows['C2'][1].real); a_raw = a
        sh = {k: P4.shift(J, h, k, xs, m) for k in (2, 3)}
        for _ in range(3):
            a = fit(rows['C2'][1].real - a*sh[2])
        rows = {kd: (k, G - a*sh[k]) for kd, (k, G) in rows.items()}
        rec = dict(m=m, xs=xs, a=a, a_raw=a_raw)
        for kd, (k, G) in rows.items():
            g = m*m*(G.real/(a*P4.env(h, xs, k)) - P4.Dc(h, J, xs))
            rec[(kd, 'raw')] = g.copy()
            if curv:
                g = g - CF.curv_image(kd, h, J, m, xs, Q.A.env, Q.A.SIG)
            rec[kd] = g
            rec[(kd, 'im')] = m*G.imag/(a*P4.env(h, xs, k))
        rec['lo_c3c2'] = lo['C3H']/lo['C2']
        dn = DIRNU.get((J, D), set())
        rec['direct'] = np.array([any(abs(x*m - v) < 1e-6*v for v in dn) for x in xs])
        out[D] = rec
    return out


NODELOG = []


def oo_even_rows(J, g):
    """OO grid -> even-row dict (++ keys = Gamma_OO, +- / off = 0); NL rows for z4 (++ = OO), flattening, Hampel."""
    out = {}
    h = '++'
    for D, r in g.items():
        m, xs = r['m'], r['xs']; nus = xs*m
        dn = r['direct']
        gre = np.array(r['C2'], float); g3 = np.array(r['C3H'], float)
        if dn.any() and Z4.DIRECT:
            sel = nus < 40.0
            gim = m*r[('C2', 'im')][sel] - m*m*Z4._ahat_im(J, h, xs[sel], m)/Q.A.env(h, xs[sel], 2)
            gc = gre[sel] + 1j*gim
            hm, hp = Z4.hhat(nus[sel])
            nl = dict(m=m, nus=nus[sel], direct=dn[sel], _curvfixed=True)
            nl['++'] = gc*(hm + hp); nl[('++', 'c')] = gc; nl[('++', 'C3H')] = g3[sel].copy()
            nl['+-'] = 0*gc; nl[('+-', 'c')] = 0*gc; nl[('+-', 'C3H')] = 0*g3[sel]
            Z4.NLROWS[(J, float(D))] = nl
            for i, v in enumerate(nus[sel]):
                NODELOG.append((J, float(D), float(v), bool(dn[sel][i]), complex(gc[i]), complex(nl['++'][i]),
                                float((Z4.yphi(v)*nl['++'][i]).real), float(g3[sel][i]), float(Z4.xphi(v).real)))
            lo = (nus < Z4.NUFLAT) & dn
            if lo.any():
                i1 = int(np.argmax(~lo & (nus >= Z4.NUFLAT)))
                gre = np.where(lo, gre[i1], gre)
        else:
            Z4.NLROWS.pop((J, float(D)), None)
        z = 0*xs
        row = {'m': m, 'xs': xs, 'a': {'++': r['a'], '+-': r['a']},
               ('C2', '++'): gre, ('C3H', '++'): g3, ('C2', '+-'): z.copy(), ('C3H', '+-'): z.copy(),
               ('C2', '++', 'im'): r[('C2', 'im')], ('C3H', '++', 'im'): r[('C3H', 'im')], ('C2', '+-', 'im'): z.copy(),
               ('C3H', '+-', 'im'): z.copy(), ('C2', 'off'): z.copy(), ('C3H', 'off'): z.copy(), ('C2', 'off', 'im'): z + 1e-30,
               ('C3H', 'off', 'im'): z.copy()}
        if ODESP:                    # w2's Hampel rule (window +-3, 6 x max(MAD, median|G|), first node kept), on the ++ (= OO) keys
            n = len(xs); bad = np.zeros(n, bool)
            for kk in (('C2', '++'), ('C3H', '++')):
                v = np.asarray(row[kk], float); scale = np.median(np.abs(v)) + 1e-12
                for i in range(n):
                    w = np.r_[v[max(0, i - 3):i], v[i + 1:i + 4]]
                    med = np.median(w); mad = np.median(np.abs(w - med))
                    if abs(v[i] - med) > 6*max(mad, scale):
                        bad[i] = True
            bad[0] = False
            if bad.any():
                W2.DESPIKE_LOG.append(('OO', J, D, [round(float(x), 4) for x in xs[bad]]))
                for kk, vv in list(row.items()):
                    if isinstance(vv, np.ndarray) and vv.shape == bad.shape:
                        row[kk] = vv[~bad]
        out[D] = row
    if DUPNL and len(out) == 1:
        D = list(out)[0]
        if (J, float(D)) in Z4.NLROWS:
            Z4.NLROWS[(J, float(D) + 1e-6)] = Z4.NLROWS[(J, float(D))]
    return out


# ------------------------------------------------------------------------------------------------ main
if __name__ == '__main__':
    t0 = time.time()
    W2.set_mode(True)
    for J in OOJ:
        if J not in Z4.JS:
            Z4.JS.append(J)
    Q.JEVEN[:] = sorted(OOJ)                  # as hc_oo: the chain runs over the OO J only (both runs)
    Q.JODD[:] = []
    pr(f"# oo_chain: OO block (even-J parity-odd X^O_same) in the z4 chain, Delta_g = {DGS}; OO_J = {OOJ}; direct OO nodes {ODIRECT} "
       f"(est {EST}, jk {OJK}, nu >= {ONUMIN}); Z4: NL {Z4.NL}, gauge {Z4.GAUGE}, extrap {Z4.EXTRAP} {os.environ.get('Z4_CONST', '')}, "
       f"NUNL {Z4.NUNL}, EE direct rows {os.environ.get('OO_Z3EE', '18:336,400;20:400;22:336')} (Z4 est {Z4.EST}); Hampel on OO rows {ODESP}; "
       f"dup NL {DUPNL}; SM_CURV {CF.MODE}; JEVEN {Q.JEVEN} JODD {Q.JODD}; rows {os.environ.get('W2_ROWS', 'all')}; "
       f"extra m {os.environ.get('Q3_EXTRA_M', '')}")
    G = {}
    for J in OOJ:
        g = grid_oo(J)
        G[J] = g
        for D, r in sorted(g.items()):
            pr(f"  J={J} Delta={D:g} m={r['m']:.2f} a_O={r['a']:.6g} (no de-shift {r['a_raw']:.6g}); nodes {len(r['xs'])} "
               f"(direct {int(r['direct'].sum())}); x = " + " ".join(f"{x:.4f}{'*' if d else ''}" for x, d in zip(r['xs'], r['direct'])))
            for kd in ('C2', 'C3H'):
                pr(f"    {kd:3s} Gamma_OO (curv-subtracted): " + " ".join(f"{v:+7.0f}" for v in r[kd]))
            c = r['lo_c3c2']
            pr(f"    LO order-0 C3H/C2 x (2-y) [++ coupling: x-independent, ~2]: " +
               " ".join(f"{(cc*(2 - x*x)).real:+.3f}" for x, cc in zip(r['xs'], c)))
    # run0: unmodified z4 chain (EE block)
    if os.environ.get('OO_RUN0', '1') == '1':
        res0 = Q.run(DGS, lambda s: None)[1]
    else:
        res0 = None
    for k in [k for k in Z4.NLROWS if k[0] in OOJ]:
        Z4.NLROWS.pop(k)
    _grid = Q.A.grid
    Q.A.grid = lambda J, K, kinds=('C2', 'C3H'), *a, **k: (oo_even_rows(J, G[J]) if J in G and G[J] else _grid(J, K, kinds, *a, **k))
    pr("# OO rows: " + "; ".join(f"J={J}: " + ",".join(f"{D:g}" for D in sorted(G[J])) for J in OOJ))
    res1 = Q.run(DGS, lambda s: None)[1]
    Q.A.grid = _grid
    store = {}
    for Dg in DGS:
        pr(f"## Delta_g = {Dg}")
        r1 = {(round(r['mfl'], 9), r['J'], r['DF']): r for r in res1[Dg]['recs']}
        r0 = {(round(r['mfl'], 9), r['J'], r['DF']): r for r in res0[Dg]['recs']} if res0 else {}
        for J in OOJ:
            ks = sorted(k for k in r1 if k[1] == J and len(r1[k]['hs']) == 2)
            rows = []
            for k in ks:
                b = r1[k]
                oo = (b['base'][0] + b['dF'][0, 0])/b['Ns'][0]
                oo0 = b['base'][0]/b['Ns'][0]
                if k in r0:
                    a = r0[k]
                    lam_ee = Q.lam(a, 1.0); pp_ee = (a['base'][0] + a['dF'][0, 0])/a['Ns'][0]; dfe = a['dF'][0, 0]/a['Ns'][0]
                else:
                    lam_ee = pp_ee = dfe = np.nan
                rows.append(dict(mfl=k[0], DF=k[2], lam_ee=lam_ee, pp_ee=pp_ee, oo=oo, oo0=oo0, dfoo=b['dF'][0, 0]/b['Ns'][0], dfpp=dfe))
            store[(Dg, J)] = rows
            if not rows:
                pr(f"  J={J}: no states"); continue
            w = min(rows, key=lambda t: t['oo'])
            pr(f"  J={J:2d}: min lambda_OO/N_s = {w['oo']:+.4f} at m/Lam={w['mfl']:.4f} ({w['DF']}) [LO {w['oo0']:+.4f}]; max |dF_OO|/N_++ "
               f"{max(abs(t['dfoo']) for t in rows):.4f} (EE ++: {np.nanmax([abs(t['dfpp']) for t in rows]):.4f}); EE-block min lam "
               f"{np.nanmin([t['lam_ee'] for t in rows]):+.4f}")
            near = [t for t in rows if t['mfl'] <= 1.3 + 1e-9]
            pr("      m/Lam [DF]: OO (dF_OO/N) | EE ++ diag (dF_++/N) | EE 2x2 min:")
            for t in near:
                pr(f"        {t['mfl']:.4f} [{t['DF']}]: {t['oo']:+.4f} ({t['dfoo']:+.4f}) | {t['pp_ee']:+.4f} ({t['dfpp']:+.4f}) | {t['lam_ee']:+.4f}")
            far = [t for t in rows if t['mfl'] > 1.3 + 1e-9]
            if far:
                pr("      beyond 1.3: " + "  ".join(f"{t['mfl']:.2f}: {t['oo']:+.4f}" for t in far))
    pr("# NL OO nodes (J, D, nu, direct, Gamma_c, Gamma_NL, NL density, C3H, Re X/phi):")
    for t in NODELOG:
        pr(f"    J={t[0]} D={t[1]:g} nu={t[2]:7.3f} {'dir' if t[3] else 'K2 '} Gc {t[4].real:+8.0f}{t[4].imag:+8.0f}i  NL {t[5].real:+8.0f}"
           f"{t[5].imag:+8.0f}i  dens {t[6]:+8.0f}  C3H {t[7]:+7.0f}  ReX/phi {t[8]:+.3f}")
    if SELLOG:
        pr("# OO estimator per node (kind, J, D, nu, chosen, jk rel err dir, cv): " + "; ".join(str(t) for t in sorted(set(SELLOG))))
    pr(f"# despiked: {W2.DESPIKE_LOG}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(C, 'oo_chain_%s.txt' % TAG), 'w').write(buf.getvalue())
    pickle.dump(dict(grids=G, states=store, nodes=NODELOG, sel=SELLOG), open(os.path.join(C, 'oo_chain_%s.pkl' % TAG), 'wb'))
