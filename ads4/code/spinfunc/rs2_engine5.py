"""
rs2_engine5 (HRS_CLOSE.md §4): one order more of the analytic-radial engine (K = 5, depth 6) for the EE k = 4 image of a state, to
measure the a^5 term directly (BELOW_FOLD2 used the size of the a^4 term as the truncation error).  The reference (Delta_ref = 2e4)
is read from the existing production cache bg_eeq_D{Delta}_J{J}_K4.pkl (read-only; its own a^k terms are ~ (1/2e4)^k).
Writes only cache/rs2_eeq_D{Delta}_J{J}_K5.pkl.
Usage: python3 rs2_engine5.py Delta J [full]   ('full': v'-u' and R from a deeper geometry -> rs2_eeq_D*_K5full.pkl)
"""
import os, sys, pickle, time
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_eeq as EQ

C = os.path.join(HERE, 'cache')


def fn(Dl, J, full=False):
    return os.path.join(C, 'rs2_eeq_D%g_J%d_K5%s.pkl' % (Dl, J, 'full' if full else ''))


def patch_full(Kdeep=9):
    """in-process patch of bg_ee.ee_pieces: v'-u' and Lemma-B R taken from a deeper geometry (bg_geom.scalar_pieces(Kdeep)), so that
    their relative a^4 corrections (P = a^4 q_x q_y, a^4 C^2, ...) are kept; at depth 5-6 the production pieces carry v'-u' and R only
    to relative order a^2-a^3 (the leading cancellation in v'-u' costs two orders)."""
    import bg_ee as BE, bg_geom as GG
    if getattr(BE, '_rs2_patched', False):
        return
    orig = BE.ee_pieces
    _, _, _, nm = GG.scalar_pieces(Kdeep)
    vmu = nm['vp'] - nm['up']
    R = (vmu*vmu - 1.0)*(nm['up']*nm['vp']).recip("u'v'").scale(0.25)
    def ee_pieces(K, Delta, J, kind='EE', verbose=False):
        P = orig(K, Delta, J, kind, verbose)
        P['vmu'], P['R'] = vmu, R
        return P
    BE.ee_pieces = ee_pieces
    BE._rs2_patched = True


def get(Dl, J, full=False):
    path = fn(Dl, J, full)
    if os.path.exists(path):
        return pickle.load(open(path, 'rb'))
    t0 = time.time()
    if full:
        patch_full()
    st = EQ.images(Dl, J, 5, 6, kinds=('EE4',))
    ref = pickle.load(open(EQ.fn(Dl, J, 4), 'rb'))['ref']
    db = dict(st=st, ref=ref)
    pickle.dump(db, open(path + '.tmp', 'wb')); os.replace(path + '.tmp', path)
    print(f'rs2_engine5 D={Dl} J={J} full={full}: {time.time() - t0:.0f}s', flush=True)
    return db


if __name__ == '__main__':
    get(float(sys.argv[1]), int(sys.argv[2]), len(sys.argv) > 3 and sys.argv[3] == 'full')
