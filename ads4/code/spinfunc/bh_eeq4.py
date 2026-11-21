"""
bh_eeq4 (BELOW_FOLD3.md §1): J >= 4 EE images on the analytic-radial engine (read-only import of bg_eeq.images, bg_ee, bg_geom, bg_ser).
Same object as bg_eeq (EE k = 4 integrand pref[a][c] x (v'-u') x R, all four (a, c) entries of the 2 x 2 prefactor matrix), stored per
(state, reference) under a bh_ name so nothing of bg_* is overwritten.
Usage: python3 bh_eeq4.py Delta J K depth [kinds=EE4]  -> cache/bh_eeq_D{Delta}_J{J}_K{K}.pkl
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import bg_eeq as EQ


def fn(Dl, J, K):
    return os.path.join(HERE, 'cache', 'bh_eeq_D%g_J%d_K%d.pkl' % (Dl, J, K))


if __name__ == '__main__':
    Dl, J, K, depth = float(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
    kinds = tuple(sys.argv[5].split(',')) if len(sys.argv) > 5 else ('EE4',)
    t0 = time.time()
    rec = EQ.images(Dl, J, K, depth, kinds=kinds)
    pickle.dump(rec, open(fn(Dl, J, K) + '.tmp', 'wb')); os.replace(fn(Dl, J, K) + '.tmp', fn(Dl, J, K))
    print(f'bh_eeq4 D={Dl} J={J} K={K} depth={depth}: done {time.time() - t0:.0f}s', flush=True)
