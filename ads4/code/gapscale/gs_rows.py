"""
gapscale/gs_rows.py (GAP_SCALING.md §4): own-mass K2 crystal rows of an ODD-J near-gap state (parity-odd +-), e.g. the binding state
J = 5 at Delta_g = 2000 (Delta = 2005), with the heavy reference Delta_ref = 2e4 at the same nu, on the dense x grid.
This is neargap2/ng_rows.py's odd branch re-typed with every write redirected to gapscale/cache/rows/ (ng_rows itself cannot be used:
its guard and caches live in neargap2/cache, which this note must not touch).  Engines (spinfunc/w1_rows, p3_odd_action, p2_shadow_odd,
w0_ee_common) are imported read-only and run unchanged.  Shadows are looked up read-only in spinfunc/cache and neargap2/cache.
A process-local guard refuses any write outside gapscale/cache.
Usage: OMP_NUM_THREADS=1 python3 gs_rows.py KINDS J Delta [x1,x2,...|dense]
"""
import os, sys, glob, pickle, time, builtins
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
NG = os.path.join(HERE, '..', 'neargap2', 'cache')
GSC = os.path.realpath(os.path.join(HERE, 'cache'))
ROWS = os.path.join(GSC, 'rows')
os.makedirs(ROWS, exist_ok=True)

# ---- guard: only gapscale/cache may be written
_open, _replace = builtins.open, os.replace


def _ok(p):
    try:
        p = os.path.realpath(os.fspath(p))
    except TypeError:
        return True
    return p == '/dev/null' or p.startswith(GSC + os.sep)


def _gopen(f, mode='r', *a, **k):
    if any(c in mode for c in 'wax+') and not _ok(f):
        raise PermissionError('gs_rows guard: refusing to write %s' % f)
    return _open(f, mode, *a, **k)


def _grep(s, d, *a, **k):
    if not _ok(d):
        raise PermissionError('gs_rows guard: refusing to move onto %s' % d)
    return _replace(s, d, *a, **k)


builtins.open = _gopen; os.replace = _grep
sys.path.insert(0, SF)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
os.chdir(SF)
import w1_rows as W1

SFC = os.path.join(SF, 'cache')
DREF = 20000.0
DENSE = [round(0.025*i, 4) for i in range(1, 37)]


def setup_odd():
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    P._fn = lambda kind, J: os.path.join(ROWS, 'w1_p3_%s_J%d_gs.pkl' % (kind, J))
    P2._DB = os.path.join(ROWS, 'gs_p2shadow_db.pkl')
    sfn = os.path.join(ROWS, 'gs_p2shadow.pkl')
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(SFC, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(SFC, 'w1_p2shadow*.pkl'))) + [os.path.join(NG, 'ng_p2shadow.pkl'), sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = orig_all(float(Delta), int(J))
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P2.shadow_odd = shadow_odd
    return P


if __name__ == '__main__':
    kinds, J, D = sys.argv[1].split(','), int(sys.argv[2]), float(sys.argv[3])
    xs = DENSE if len(sys.argv) < 5 or sys.argv[4] == 'dense' else [float(v) for v in sys.argv[4].split(',')]
    assert J % 2 == 1
    P = setup_odd()
    m = W1.mphys(D, J)
    t0 = time.time(); n = 0
    for kind in kinds:
        for x in xs:
            for DD in (D, DREF):
                fn = P._fn(kind, J)
                R = {(k[0], k[1]): v for k, v in (pickle.load(open(fn, 'rb')).items() if os.path.exists(fn) else [])}
                if not W1.have(R, DD, x*m):
                    t1 = time.time()
                    P.get(kind, J, DD, [x*m], verbose=False); n += 1
                    print(f"  {kind} J={J} D={DD:g} x={x:.4f} nu={x*m:.3f}: {time.time() - t1:.1f}s", flush=True)
    print(f"done odd {','.join(kinds)} J={J} Delta={D:.6g} m={m:.3f} J^2/m={J*J/m:.3f}: {n} new points in {time.time() - t0:.0f}s",
          flush=True)
