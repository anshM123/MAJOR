"""
hs_points (HIGHSPIN.md): DIRECT crystal points (spinfunc/hc_points.py, unchanged, imported read-only; n1 adaptive re-centred contour,
RQMC, 8 replicas, state and reference Delta_ref = 2e4 on common random numbers; odd J -> parity-odd 1x1 integrand, even J -> w0 2x2 EE
integrand + C3H) -- exactly neargap2/ng_points.py, with the cache moved to highspin/cache (hc_points' own HC_CACHE switch).
Existing direct points of the same state in neargap2/cache are copied into highspin/cache first (read there, written here), so they are
not recomputed.
x-list: comma list of x = nu/m(Delta, J), or 'gk15:a:b' (the 15 Gauss-Kronrod nodes on [a, b], Gauss-7 nodes first), or a mix joined
by '+'; 'gk15:a:b:xcap' skips the GK nodes above xcap.
Usage: python3 hs_points.py J Delta x-list      ('seed' as x-list: only copy the existing points)
"""
import os, sys, time, pickle, shutil
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import hs_guard  # noqa: F401  (first)
SF = os.path.join(HERE, '..', 'spinfunc')
NG = os.path.join(HERE, '..', 'neargap2', 'cache')
HSC = os.path.join(HERE, 'cache')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import numpy as np

GK_X = [0.949107912342759, 0.741531185599394, 0.405845151377397, 0.0,
        0.991455371120813, 0.864864423359769, 0.586087235467691, 0.207784955007898]


def gk15(a, b):
    c, h = 0.5*(a + b), 0.5*(b - a)
    out = []
    for t in GK_X:
        out += [c - h*t] if t == 0 else [c - h*t, c + h*t]
    return out


def parse_x(spec):
    xs = []
    for part in spec.split('+'):
        if part.startswith('gk15:'):                     # gk15:a:b[:xcap]  (nodes above xcap skipped: negligible weight)
            p = part.split(':'); cap = float(p[3]) if len(p) > 3 else 9.0
            xs += [x for x in gk15(float(p[1]), float(p[2])) if x <= cap]
        else:
            xs += [float(v) for v in part.split(',') if v]
    return xs


def seed_from_neargap2(J, D):
    """copy existing direct points of this state (neargap2/cache first, else spinfunc/cache: HC's hc_pts) into highspin/cache."""
    for kind in ('EE', 'C3H'):
        f = 'hc_pts_%s_J%d_D%g.pkl' % (kind, J, D)
        dst = os.path.join(HSC, f)
        for d in (NG, os.path.join(SF, 'cache')):
            src = os.path.join(d, f)
            if os.path.exists(src) and not os.path.exists(dst):
                pickle.dump(pickle.load(open(src, 'rb')), open(dst, 'wb'))
                print(f"# copied {f} from {os.path.relpath(d, HERE)}", flush=True)


if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2])
    seed_from_neargap2(J, D)
    if sys.argv[3] == 'seed':
        sys.exit(0)
    os.chdir(SF)
    os.environ['HC_CACHE'] = HSC
    for k, v in (('HC_METHOD', 'n1'), ('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')):
        os.environ.setdefault(k, v)
    import warnings
    warnings.filterwarnings('ignore')
    import hc_points as HP
    HP.STATE[0] = D
    m = HP.mphys(D, J)
    xs = parse_x(sys.argv[3])
    nus = [x*m for x in xs]
    HP.setup(odd=(J % 2 == 1))
    t0 = time.time()
    print(f"# hs_points [n1] J={J} Delta={D:g} m={m:.3f} x={[round(x, 6) for x in xs]} -> {HP.fn('EE', J)}", flush=True)
    HP.run_n1(J, nus, D, HP.Z2.MQ, HP.Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
