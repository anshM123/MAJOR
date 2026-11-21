"""
bj_eeq (BELOW_FOLD_SPINJ.md §2): production of closed-form spin-J EE images (bj_ee) for a state (Delta, J) and its reference
(Delta_ref = 2e4, same J), on a nu grid, kind EE2 (k = 2: C^'_EE, the Result D kind) [and EE4 on request].
Stores, per (Delta, J): images for several (N, Kr) truncations (convergence checks), keyed res[(N, Kr)] = bj_ee.images output.
Usage: python3 bj_eeq.py Delta J depth N[,N2,...] [nus-spec] [kinds]
   nus-spec: 'win' (default: 0.25..12 step 0.25 + the z2 nodes 8.33..33.33) | 'wide' (0.25..30 step 0.25) | comma list
   -> cache/bj_eeq_D{Delta}_J{J}_d{depth}.pkl   (reference: same file name with D20000; reused when present)
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bj_ee as BJ

C = os.path.join(HERE, 'cache')
DREF = 2e4
Z2NODES = [8.332091574, 9.998509889, 11.664928204, 13.331346519, 16.664183148, 20.830228935, 24.996274722, 33.328366297]


def grid(spec):
    if spec == 'win':
        return np.array(sorted(set(np.round(np.arange(0.25, 12.001, 0.25), 6).tolist() + [4.999254944, 6.665673259] + Z2NODES)))
    if spec == 'wide':
        return np.round(np.arange(0.25, 30.001, 0.25), 6)
    return np.array([float(v) for v in spec.split(',')])


def fn(Dl, J, depth):
    return os.path.join(C, 'bj_eeq_D%g_J%d_d%d.pkl' % (Dl, J, depth))


def run(Dl, J, depth, Ns, nus, kinds=('EE2',), verbose=True):
    path = fn(Dl, J, depth)
    db = pickle.load(open(path, 'rb')) if os.path.exists(path) else dict(Delta=Dl, J=J, depth=depth, res={})
    todo = [(N, Kr) for N in Ns for Kr in (depth - 1, depth - 2) if Kr >= 0 and (N, Kr) not in db['res']]
    if not todo:
        return db
    t0 = time.time()
    P = BJ.pieces(depth, Dl, J, verbose=verbose)
    db['pieces_time'] = time.time() - t0
    db['flip'] = P['flip']
    for (N, Kr) in todo:
        db['res'][(N, Kr)] = BJ.images(P, N, Kr, 1000, nus=nus, kinds=kinds, verbose=verbose)
        pickle.dump(db, open(path + '.tmp', 'wb')); os.replace(path + '.tmp', path)
    if verbose:
        print(f'bj_eeq D={Dl:g} J={J} depth={depth}: {time.time() - t0:.0f}s', flush=True)
    return db


if __name__ == '__main__':
    Dl, J, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    Ns = [int(v) for v in sys.argv[4].split(',')]
    nus = grid(sys.argv[5] if len(sys.argv) > 5 else 'win')
    kinds = tuple(sys.argv[6].split(',')) if len(sys.argv) > 6 else ('EE2',)
    run(Dl, J, depth, Ns, nus, kinds)
