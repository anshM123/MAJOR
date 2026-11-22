"""
rs_engine (HRS_CHECK.md): the analytic-radial engine's EE k = 4 image (bg_eeq.images, read-only import, unchanged) for states that
have no bg_eeq cache (e.g. the Delta = 100 sensitivity state), written to the own cache cache/rs_eng_D<Delta>_J<J>_K<K>.pkl in
bg_eeq's {st, ref} format; and gamma_K: bg_qan.gamma_K without bg_qan's T cache side effect (T from bf_qan.T, cached in rs_T_*).
Usage: python3 rs_engine.py Delta J [K depth Kref]
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'resultq'))
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')


def fn(Dl, J, K):
    return os.path.join(C, 'rs_eng_D%g_J%d_K%d.pkl' % (Dl, J, K))


def get(Dl, J, K=4, depth=5, Kref=2):
    import bg_eeq as EQ
    if os.path.exists(EQ.fn(Dl, J, K)):
        return pickle.load(open(EQ.fn(Dl, J, K), 'rb'))
    p = fn(Dl, J, K)
    if os.path.exists(p):
        return pickle.load(open(p, 'rb'))
    st = EQ.images(Dl, J, K, depth)
    rf = EQ.images(EQ.DREF, J, Kref, Kref + (depth - K))
    db = dict(st=st, ref=rf)
    pickle.dump(db, open(p + '.tmp', 'wb')); os.replace(p + '.tmp', p)
    return db


def tval(nus, D, J, R):
    import bf_qan as BQ
    p = os.path.join(C, 'rs_T_D%g_J%d.pkl' % (D, J))
    if os.path.exists(p):
        d = pickle.load(open(p, 'rb'))
        if len(d) == len(nus):
            return d
    v = BQ.T(nus, BQ.mass(D, J)[1], R)
    pickle.dump(v, open(p, 'wb'))
    return v


def tcache_ro(nus, D, J, R):
    """bg_qan's T cache if it exists (read-only), else ours."""
    p = os.path.join(C, 'bg_qan_T_D%g_J%d.pkl' % (D, J))
    if os.path.exists(p):
        d = pickle.load(open(p, 'rb'))
        if len(d) == len(nus):
            return d
    return tval(nus, D, J, R)


def gamma_K(db, h='++', kind='EE4', Kt=None):
    """= bg_qan.gamma_K (same formulas), T read from bg_qan's cache when present, else from rs_T_*."""
    import bf_qan as BQ
    import bg_eeq as EQ
    import bg_qan as BG
    st, rf = db['st'], db['ref']
    D, J = st['Delta'], st['J']
    n = 1024
    nus = st['nus'][:n]
    m, _ = BQ.mass(D, J)
    R = BQ.Rj(nus)
    Tst = tcache_ro(nus, D, J, R); Tref = tcache_ro(nus, EQ.DREF, J, R)
    NH = sum(R.values())
    HYB = BG.GAM[h]*sum(R[j]*(nus - 1j*j)**2 for j in R)
    Ist, _ = EQ.value(st, kind, 0, 0, Kt)
    Iref, _ = EQ.value(rf, kind, 0, 0)
    Ist, Iref = Ist[:n], Iref[:n]
    tot = Tref*Ist/Iref
    cc = Tst[0]/tot[0]
    Gex = m*m*(cc*tot - Tst)
    GK = Gex - HYB
    GK = GK - (GK[0]/NH[0])*NH
    return dict(nus=nus, GK=GK, Gex=Gex, m=m, NH=NH, Tst=Tst, Tref=Tref, ratio=Ist/Iref)


if __name__ == '__main__':
    Dl, J = float(sys.argv[1]), int(sys.argv[2])
    K, depth, Kref = (int(v) for v in sys.argv[3:6]) if len(sys.argv) > 5 else (4, 5, 2)
    t = time.time()
    get(Dl, J, K, depth, Kref)
    print(f'rs_engine D={Dl} J={J} K={K}: {time.time() - t:.0f}s', flush=True)
