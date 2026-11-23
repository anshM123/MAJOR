"""
ef_shadow (EXT_FOLD.md §2): the parity-odd <TTO> shadow coefficient for a state with no K2 row (J = 7, Delta = 1007), needed to turn
its direct points into chain records (hc_nodes.to_record: M0 = mean replica / b_OO, logK = the shadow's logK).  Computed with
spinfunc/p2_shadow_odd.shadow_all (read-only import; defaults xi = 0.5, 8 triples, seed 7 -- exactly what every cached odd shadow used:
p2_shadow_odd.shadow_odd / w1_rows.setup_odd call it with defaults), after reproducing two cached shadows bit for bit.
Usage: python3 ef_shadow.py -> cache/ef_p2shadow.pkl, cache/ef_shadow.txt
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc'); sys.path.insert(0, SF)
os.chdir(SF)
import glob
import numpy as np
import p2_shadow_odd as P2
EFC = os.path.join(HERE, 'cache')
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def cached(key):
    for fn in ([os.path.join(SF, 'cache', 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(SF, 'cache', 'q2_p2shadow*.pkl'))) +
               sorted(glob.glob(os.path.join(SF, 'cache', 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(SF, 'cache', 'ee1000_p2shadow*.pkl')))):
        db = pickle.load(open(fn, 'rb'))
        if key in db:
            return db[key], os.path.basename(fn)
    return None, None


if __name__ == '__main__':
    out = {}
    for D, J in ((800.0, 7), (1023.0, 23), (1007.0, 7)):
        t0 = time.time()
        r = P2.shadow_all(D, J)
        key = (D, J, P2.T.DEG)
        c, src = cached(key)
        if c is not None:
            same = np.array_equal(np.asarray(c[0]), np.asarray(r['b'])) and complex(c[1]) == complex(r['logK'])
            dev = float(np.abs(np.asarray(c[0]) - np.asarray(r['b'])).max()/np.abs(np.asarray(c[0])).max())
            pr(f"  check {key}: cached ({src}) vs recomputed: bitwise {same}; rel dev b {dev:.1e}; dlogK {abs(complex(c[1]) - complex(r['logK'])):.1e}")
        else:
            out[key] = (r['b'], r['logK'], r['resid'])
            pr(f"  new {key}: b_OO = {complex(np.asarray(r['b'])[-1, -1]):.6e}, logK = {complex(r['logK']):.6f}, resid {r['resid']}")
        pr(f"   [{time.time() - t0:.0f}s]")
    pickle.dump(out, open(os.path.join(EFC, 'ef_p2shadow.pkl'), 'wb'))
    open(os.path.join(EFC, 'ef_shadow.txt'), 'w').write("\n".join(lines) + "\n")
