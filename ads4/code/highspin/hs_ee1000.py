"""
hs_ee1000 (HIGHSPIN.md §2): the published ee1000 shape chain (spinfunc/ee1000_chain.run, in-process, unchanged; sm_curvfix row/kappa
imported first, as sm_chainfix does) for the listed J, once as published and once with the K2 nodes at the reference zero
(hs_core.refbad) removed from every row.  Nothing is written outside highspin/cache.
Usage: [EE1000_JEVEN=.. EE1000_JODD=.. EE1000_K=.. EE1000_XMIN=0.05 W2_DMIN=1000] python3 hs_ee1000.py tag
"""
import os, sys, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import hs_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF); os.chdir(SF)
for k, v in (('L1_DEG', '6'), ('W2_BASIS', 'ee'), ('W2_MODES', 'ee'), ('SM_CURV', 'row'), ('W2_DMIN', '1000'), ('EE1000_XMIN', '0.05')):
    os.environ.setdefault(k, v)
import numpy as np
import sm_curvfix as CF  # noqa: F401
import ee1000_chain as CH
import hs_core as H      # (loaders of ng_chain added; refbad)
tag = sys.argv[1]
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(s + "\n")


recs0, _ = CH.run(lambda s: None)
_eg, _og = CH.Q.even_grid, CH.Q.odd_grid
DROPPED = []


def _drop(g, J):
    out = {}
    for D, r in g.items():
        bad = H.refbad(J, D)
        keep = np.array([round(x, 6) not in bad for x in r['xs']])
        if (~keep).any():
            DROPPED.append((J, D, [float(x) for x in r['xs'][~keep]]))
        n = len(r['xs'])
        out[D] = {k: (v[keep] if isinstance(v, np.ndarray) and v.shape == (n,) else v) for k, v in r.items()}
    return out


CH.Q.even_grid = lambda J: _drop(_eg(J), J)
CH.Q.odd_grid = lambda J: _drop(_og(J), J)
recs1, _ = CH.run(lambda s: None)
pr(f"# hs_ee1000 {tag}: env " + " ".join(f"{k}={v}" for k, v in sorted(os.environ.items()) if k.startswith(('EE1000', 'W2_', 'SM_'))))
pr(f"# reference-zero nodes removed: {DROPPED}")
for a, b in zip(recs0, recs1):
    pr(f"  J={a['J']:2d} m/Lam={a['mfl']:.4f}: published {a['lam']:+.4f} -> reference-zero node removed {b['lam']:+.4f}")
open(os.path.join(HERE, 'cache', 'hs_ee1000_%s.txt' % tag), 'w').write(buf.getvalue())
