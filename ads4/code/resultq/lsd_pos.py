"""
lsd_pos (LIGHTSIDE_DIST.md): reference whole-function values Omega_pos = E[piece x class](a, b) and (b, a) from EIN's engine (ein_run.Etable,
read-only import), for k = 2, 4, 5, at the angle pairs / polarisations used by lsd_run.  Output: cache/lsd_pos_<tag>.pkl.

    (ulimit -v 3000000; OMP_NUM_THREADS=1 python3 lsd_pos.py <tag> [n])      tags: P1, P2 (see CONFIGS)
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


def rnull(rng):
    v = rng.normal(size=2); v /= np.linalg.norm(v); s = rng.uniform(0.5, 2.0)
    return np.array([s, s*v[0], s*v[1]])


def config(tag):
    """(phi1, phi3, p2, p4) for a tag."""
    seeds = {'P1': (7, 0.7, -0.5), 'P2': (8, -2.3, 1.8), 'P3': (9, 0.25, -0.05)}
    s, a, b = seeds[tag]
    rng = np.random.default_rng(s)
    p2, p4 = rnull(rng), rnull(rng)
    return a, b, p2, p4


if __name__ == '__main__':
    import ein_run as R
    tag = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    a, b, p2, p4 = config(tag)
    t0 = time.time(); out = {'a': a, 'b': b, 'p2': p2, 'p4': p4, 'n': n}
    ks = tuple(int(x) for x in sys.argv[3].split(",")) if len(sys.argv) > 3 else (4, 2, 5)
    for k in ks:
        out[(k, 'ab')] = R.Etable(k, a, b, p2, p4, n=n)
        out[(k, 'ba')] = R.Etable(k, b, a, p2, p4, n=n)
        print(f'k = {k} done [{time.time() - t0:.0f}s]', flush=True)
    with open(os.path.join(HERE, 'cache', f'lsd_pos_{tag}_n{n}' + ('' if len(sys.argv) <= 3 else '_k' + sys.argv[3].replace(',', '')) + '.pkl'), 'wb') as fh:
        pickle.dump(out, fh)
    print('saved', flush=True)
