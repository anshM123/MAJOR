"""
hc_nodes (HYPOTHESES_CLOSE.md): the hc_points direct points turned into chain records (K2-record format of the l9/q3/w2/p4 chains:
M0 = direct all-orders matrix E b^-T (even; b the <TTO> shadow matrix) or E/b_OO (odd, 1x1), M2 = M4 = M6 = 0, logC = logN,
logK = the shadow's logK), for the state AND its own reference (both direct; one file per state, so no reference key is ever shared
between two states).  The shadow b used for a state/reference is the one whose stored fit residual matches the 'shadow_resid' of the
existing K2 record of the same (Delta, J) in the row the point is injected into (the K2 rows store M = S b^-T with whatever b was
current; fw_fix's rule) -- so the direct records are in exactly the structure/shadow basis of the row.
Estimators: 'dir' (plain ratio), jackknife over replicas (jk = replica left out).
Functions: records(kind, J, D, jk=None, xmin=0, xmax=9) -> {(Delta, nu): rec};  inject(R, kind, J, rows) (adds/overrides keys);
           shadow_even(D, J, resid), shadow_odd(D, J, resid).
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
LOG = []


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def pfile(kind, J, D):
    return os.path.join(C, 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D))


def pts(kind, J, D):
    f = pfile(kind, J, D)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


_SC = {}


def _cands(odd):
    key = 'odd' if odd else 'even'
    if key not in _SC:
        if odd:
            fs = ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                  sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl'))) +
                  sorted(glob.glob(os.path.join(C, 'hc_p2shadow*.pkl'))))
        else:
            fs = ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
                  [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
                  sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_shadow_k2*.pkl'))) +
                  sorted(glob.glob(os.path.join(C, 'y1_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'hc_shadow_k2*.pkl'))))
        d = {}
        for f in fs:
            if os.path.exists(f):
                for k, v in pickle.load(open(f, 'rb')).items():
                    d.setdefault(k, []).append((os.path.basename(f), v))
        _SC[key] = d
    return _SC[key]


def shadow(D, J, resid=None):
    odd = J % 2 == 1
    k = (float(D), int(J), 6) if odd else (float(D), int(J))
    c = _cands(odd).get(k, [])
    if not c:
        raise KeyError(('no shadow', D, J))
    if resid is not None:
        rr = np.asarray(resid, float).ravel()
        for f, v in c:
            if np.allclose(np.asarray(v[2], float).ravel(), rr, rtol=1e-10, atol=0):
                return v, f
        LOG.append(('shadow: no resid match, first candidate used', D, J, c[0][0]))
    return c[0][1], c[0][0]


RESID = {}      # (D, J, kind) -> shadow_resid of the K2 records of the row (set by set_resid from the loaded K2 records)


def set_resid(R, J, kind):
    for (D, nu), r in R.items():
        if 'shadow_resid' in r and (float(D), J, kind) not in RESID:
            RESID[(float(D), J, kind)] = r['shadow_resid']


def to_record(r, D, J, kind, sel=None):
    (b, logK, res), src = shadow(D, J, RESID.get((float(D), J, kind)))
    reps = r['reps'] if sel is None else r['reps'][sel]
    if J % 2 == 1:
        bi = 1.0/np.asarray(b)[-1, -1]
        H = np.array([[np.mean([rep[0, 0] for rep in reps])*bi]])
    else:
        binv = np.linalg.inv(b)
        na = r['E'].shape[0]
        H = np.array([rep[:na] @ binv.T for rep in reps]).mean(0)
    Z = 0*H
    return dict(M0=H, M2=Z, M4=Z, M6=Z, logC=complex(r['logN']), logK=complex(logK), direct=True, hc=True, shsrc=src,
                m=mphys(D, J), shadow_resid=res)


def records(kind, J, D, jk=None, xmin=0.0, xmax=9.0, numin=0.0, numax=1e9):
    P = pts(kind, J, D)
    m = mphys(D, J)
    out = {}
    nus = sorted(set(k[1] for k in P))
    for nu in nus:
        if (D, nu) not in P or (DREF, nu) not in P or not (xmin - 1e-9 <= nu/m <= xmax + 1e-9) or not (numin <= nu <= numax):
            continue
        R = len(P[(DREF, nu)]['reps'])
        sel = None if jk is None else (np.arange(R) != jk)
        for DD in (D, DREF):
            out[(float(DD), round(float(nu), 9))] = to_record(P[(DD, nu)], DD, J, kind, sel)
    return out


def inject(R, kind, J, rows, jk=None, **kw):
    """R: existing K2 records {(D, nu): rec} of (kind, J) (keys may carry a third n0 entry for odd p3 caches -- caller strips);
    rows: list of state Deltas whose hc points are injected.  Reference keys must not collide between rows."""
    set_resid(R, J, kind)
    R = dict(R)
    seen = {}
    for D in rows:
        add = records(kind, J, D, jk, **kw)
        for k, v in add.items():
            if k[0] == DREF:
                assert seen.get(k, D) == D, ('reference key shared between rows', k, seen.get(k), D)
                seen[k] = D
            if k in R:
                LOG.append(('override', kind, J, k[0], round(k[1], 4)))
            R[k] = v
    return R


def rel_errors(kind, J, D):
    """per node: x, max relative replica error of the state and of the reference matrix."""
    P = pts(kind, J, D); m = mphys(D, J); out = []
    for (DD, nu), r in sorted(P.items()):
        if DD == D and (DREF, nu) in P:
            e1 = float(np.max(np.abs(r['err'])/np.maximum(np.abs(r['E']), 1e-300)))
            rr = P[(DREF, nu)]
            e2 = float(np.max(np.abs(rr['err'])/np.maximum(np.abs(rr['E']), 1e-300)))
            out.append((nu/m, nu, e1, e2, r.get('method', '')))
    return out
