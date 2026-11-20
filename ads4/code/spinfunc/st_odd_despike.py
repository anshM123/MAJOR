"""Run sm_chainfix ee1000 with w2's Hampel despike rule extended to the ODD rows (keys 'C2', 'C3H'; same rule as w2_chain._grid:
window +-3, |G - median| > 6 max(MAD, median|G|), first node never dropped). In-process patch only; inserted below sm_curvfix
(stage 'post' semantics: despike first, then the curvature image is subtracted on the surviving nodes)."""
import os, sys, runpy
sys.dont_write_bytecode = True
H = '/home/user/MAJOR/ads4/code/spinfunc'; sys.path.insert(0, H); os.chdir(H)
os.environ.setdefault('L1_DEG', '6'); os.environ.setdefault('W2_BASIS', 'ee'); os.environ.setdefault('W2_MODES', 'ee')
import numpy as np
import w2_chain as W2
Q = W2.Q
LOGD = []
_orig = Q.P4.grid


def grid_despiked(J, *a, **k):
    g = _orig(J, *a, **k)
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
        if bad.any():
            LOGD.append((J, float(D), [round(float(x), 4) for x in r['xs'][bad]]))
            for kk, vv in list(r.items()):
                if isinstance(vv, np.ndarray) and vv.shape == bad.shape:
                    r[kk] = vv[~bad]
    return g


Q.P4.grid = grid_despiked
import atexit
atexit.register(lambda: print("# odd despike log (J, Delta, x dropped):", LOGD, flush=True))
sys.argv = ['sm_chainfix.py'] + (sys.argv[1:] or ['ee1000'])
runpy.run_path(os.path.join(H, 'sm_chainfix.py'), run_name='__main__')
