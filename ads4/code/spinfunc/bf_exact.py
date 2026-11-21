"""
bf_exact (BELOW_FOLD.md): exact scalar comparison records (n4_exact.record, read-only import of the CHMRS/Mellin skeleton
exactskel.skelexact.ExactState) for the states used by the below-fold validation; own cache cache/bf_exact.pkl.
Usage: python3 bf_exact.py D:J [D:J ...]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

PK = os.path.join(HERE, 'cache', 'bf_exact.pkl')


def load():
    db = {}
    for f in ('n4_exact.pkl', 'z3_exact.pkl', 'bf_exact.pkl'):
        fn = os.path.join(HERE, 'cache', f)
        if os.path.exists(fn):
            db.update(pickle.load(open(fn, 'rb')))
    return db


def get(D, J):
    db = load()
    key = (round(float(D), 9), int(J))
    if key in db:
        return db[key]
    import n4_exact as N4
    t = time.time()
    r = N4.record(float(D), int(J))
    own = pickle.load(open(PK, 'rb')) if os.path.exists(PK) else {}
    own[key] = r
    pickle.dump(own, open(PK + '.tmp', 'wb')); os.replace(PK + '.tmp', PK)
    print(f'bf_exact: D={D} J={J}: {time.time() - t:.0f}s', flush=True)
    return r


if __name__ == '__main__':
    for a in sys.argv[1:]:
        D, J = a.split(':')
        get(float(D), int(J))
