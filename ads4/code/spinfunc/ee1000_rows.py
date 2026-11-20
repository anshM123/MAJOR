"""
ee1000_rows (EE_DG1000_LARGEJ.md): K2 (L1_DEG = 6) crystal rows for the Delta_g = 1000 large-J near-gap completion; a thin wrapper of
w1_rows (EE_NEARGAP.md) with its own cache names, so that nothing of w1/w2 is overwritten:
  even J: kinds EE (w0 integrand), C3H (o_H q), refb (EE + q^2 at Delta_ref, x_ref = 0.3 ... 0.8, for the EE-own helicity basis)
          -> cache/ee1000_pts_<kind>_J<J>_k2.pkl                 (l9_points engine, state + reference Delta_ref = 2e4 at the same nu)
  odd  J: kinds EE, C3H (p3_odd_action engine, 1x1 parity-odd +-, n0 = +1) -> cache/ee1000_p3_<kind>_J<J>.pkl
  shadows: looked up in every existing K2 shadow cache; new ones -> cache/ee1000_shadow_k2_<Q>.pkl / cache/ee1000_p2shadow_<Q>.pkl
  (Q = EE1000_Q, one file per queue: no two processes write the same cache file).
Resumable at point granularity (each point is pickled as soon as it is computed; existing points in w1/u5/q2/l9/ee1000 caches are
skipped).
Usage: EE1000_Q=A OMP_NUM_THREADS=1 python3 ee1000_rows.py even|odd EE,C3H J Delta x1,x2,...|X35|X39|X41|LO|HI|LOHI
       EE1000_Q=A python3 ee1000_rows.py refb J1,J2,...
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import w1_rows as W1

C = os.path.join(HERE, 'cache')
QN = os.environ.get('EE1000_Q', 'A')
DREF = 20000.0
XS19 = W1.XS19
MID = [0.0625, 0.0875, 0.1125, 0.1375, 0.225, 0.275, 0.325, 0.375, 0.425, 0.475, 0.525, 0.575, 0.625, 0.675, 0.725, 0.775]
HI = [0.825, 0.85, 0.875, 0.9]
LO = [0.025, 0.0375]
XSETS = {'X35': sorted(XS19 + MID), 'X39': sorted(XS19 + MID + HI), 'X41': sorted(LO + XS19 + MID + HI), 'MIDHI': sorted(MID + HI),
         'HI': HI, 'LO': LO, 'LOHI': LO + HI}
mphys = W1.mphys


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(C, 'ee1000_pts_%s_J%d_k2.pkl' % (kind, J))
    sfn = os.path.join(C, 'ee1000_shadow_k2_%s.pkl' % QN)

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(C, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(C, 'u5_shadow_k2.pkl'), os.path.join(C, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(C, 'w1_shadow_k2*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_shadow_k2*.pkl')))):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P.shadow = shadow
    return P


def setup_odd():
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    P._fn = lambda kind, J: os.path.join(C, 'ee1000_p3_%s_J%d.pkl' % (kind, J))
    sfn = os.path.join(C, 'ee1000_p2shadow_%s.pkl' % QN)
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(C, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(C, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(C, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(C, 'ee1000_p2shadow*.pkl')))):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        r = orig_all(float(Delta), int(J))
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        return db[key]
    P2.shadow_odd = shadow_odd
    return P


def existing_even(kind, J):
    R = W1.existing_even(kind, J)
    for fn in sorted(glob.glob(os.path.join(C, 'ee1000_pts_%s_J%d_k2*.pkl' % (kind, J)))):
        R.update(pickle.load(open(fn, 'rb')))
    return R


def existing_odd(kind, J):
    R = {}
    for p in ('p3_pts_%s_J%d.pkl', 'q2_p3_%s_J%d_*.pkl', 'w1_p3_%s_J%d*.pkl', 'ee1000_p3_%s_J%d.pkl'):
        for fn in sorted(glob.glob(os.path.join(C, p % (kind, J)))):
            R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if k[2] == 1.0})
    return R


have = W1.have

if __name__ == '__main__':
    par = sys.argv[1]
    if par == 'refb':
        P = setup_even()
        for J in [int(a) for a in sys.argv[2].split(',')]:
            t0 = time.time()
            mr = mphys(DREF, J)
            for kind in ('C2', 'EE'):
                R = existing_even(kind, J)
                nus = [x*mr for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8) if not have(R, DREF, x*mr)]
                if nus:
                    P.get(kind, J, DREF, nus, verbose=False)
            print(f"done refb J={J}: {time.time() - t0:.0f}s", flush=True)
        sys.exit(0)
    kinds = sys.argv[2].split(',')
    Js = [int(a) for a in sys.argv[3].split(',')]
    Ds = [float(a) for a in sys.argv[4].split(',')]
    xarg = sys.argv[5]
    xs = XSETS[xarg] if xarg in XSETS else [float(a) for a in xarg.split(',')]
    P = setup_even() if par == 'even' else setup_odd()
    for J in Js:
        for D in Ds:
            t0 = time.time(); n = 0
            m = mphys(D, J)
            for x in xs:
                for kind in kinds:
                    for DD in (D, DREF):
                        R = existing_even(kind, J) if par == 'even' else existing_odd(kind, J)
                        if not have(R, DD, x*m):
                            t1 = time.time()
                            P.get(kind, J, DD, [x*m], verbose=False); n += 1
                            print(f"  {par} {kind} J={J} Delta={DD:.0f} x={x}: {time.time() - t1:.0f}s", flush=True)
            print(f"done {par} {','.join(kinds)} J={J} Delta={D:.4g} m={m:.1f} J^2/m={J*J/m:.2f}: {n} new points in "
                  f"{time.time() - t0:.0f}s", flush=True)
