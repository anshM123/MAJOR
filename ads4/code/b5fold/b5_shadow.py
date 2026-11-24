"""
b5_shadow (BATCH5_FOLD.md §1): the even <TTO> shadow matrices b(1000 + J, J) of the even own rows (J = 18 ... 32), which no K2 row
provides (the even chain's rows for these J are 336 / 400 / 500 / 800 / 1040 / 1600; J = 26, 30 have none).  Exactly the function
ngfill/nf_rows.py (= neargap2/ng_rows.py, w1_rows) calls when it builds an even K2 row: spinfunc/l5_shadow_matrix.shadow_matrix(Delta,
J, 0.5, ntrip=6, seed=7), cached as (b, logK, resid) under key (Delta, J) in cache/b5_shadow_k2.pkl (hc_nodes' even shadow format).
Check: the same call reproduces an existing project shadow bit for bit (e.g. (1040, 32) in spinfunc/cache/w1_shadow_k2_a.pkl), and the
EDGE_R2 agent's (1026, 26), (1030, 30) (edger2/cache/k2/w1_shadow_k2_e2.pkl, same function) are compared.
Usage: python3 b5_shadow.py check | J,J,...   -> cache/b5_shadow_k2.pkl
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import b5_guard  # noqa: F401
os.environ['L1_DEG'] = '6'
import warnings
warnings.filterwarnings('ignore')
import numpy as np
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF)
cwd = os.getcwd(); os.chdir(SF)
import l5_shadow_matrix as SM
os.chdir(cwd)
SHF = os.path.join(HERE, 'cache', 'b5_shadow_k2.pkl')


def comp(a, b):
    return max(float(np.max(np.abs(np.asarray(a[0]) - np.asarray(b[0])))), abs(complex(a[1]) - complex(b[1])))


if __name__ == '__main__':
    if sys.argv[1] == 'check':
        ref = pickle.load(open(os.path.join(SF, 'cache', 'w1_shadow_k2_a.pkl'), 'rb'))
        for key in [(1040.0, 32)]:
            t = time.time(); r = SM.shadow_matrix(key[0], key[1], 0.5, ntrip=6, seed=7)
            print(f"reproduce {key}: max |d(b, logK)| = {comp((r['b'], r['logK']), ref[key]):.3e} [{time.time() - t:.1f}s]", flush=True)
        sys.exit(0)
    db = pickle.load(open(SHF, 'rb')) if os.path.exists(SHF) else {}
    e2 = pickle.load(open(os.path.join(HERE, '..', 'edger2', 'cache', 'k2', 'w1_shadow_k2_e2.pkl'), 'rb'))
    for J in [int(v) for v in sys.argv[1].split(',')]:
        key = (1000.0 + J, J)
        if key not in db:
            t = time.time(); r = SM.shadow_matrix(key[0], J, 0.5, ntrip=6, seed=7)
            db[key] = (r['b'], r['logK'], r['resid'])
            pickle.dump(db, open(SHF + '.tmp', 'wb')); os.replace(SHF + '.tmp', SHF)
            print(f"shadow {key}: resid {np.max(np.abs(r['resid'])):.2e} [{time.time() - t:.1f}s]", flush=True)
        if key in e2:
            print(f"   vs edger2 {key}: max |d(b, logK)| = {comp(db[key], e2[key]):.3e}", flush=True)
