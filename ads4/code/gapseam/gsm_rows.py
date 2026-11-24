"""
gsm_rows (GAP_SEAM.md): own-mass K2 crystal rows (L1_DEG = 6) of a near-gap state at Delta_g = 2000 (Delta = 2000 + J), on the dense
x grid 0.025 ... 0.9, with the heavy reference at a PARAMETERISED Delta_ref (env GSM_DREF, default 5e4; GAP_NEARGAP N1: the 2e4
reference's leading partial-wave zero x0 = 2.405 m_ref/((2J+1) m_row) falls inside the 2000 rows for J >= 16).
Copy of gapnear/gn_rows.py (odd J = gapscale/gs_rows.py, even J = neargap2/ng_rows.py's even branch) with:
  * DREF a parameter; cache paths moved to gapseam/cache/rows (file names *_gsm.pkl);
  * existing state records of gapnear/cache/rows are reused (read-only) when present (J = 16, 21 re-runs: only the new reference points
    are computed); have() also looks there;
  * 'refb' mode (even J): six EE points at the reference itself, x_ref = 0.3 ... 0.8 (l8_analysis.basis needs >= 4 points with
    x_ref >= 0.15 to build the reference helicity basis; w1_rows.py's 'refb' at the new Delta_ref);
  * 'oo' mode (even J): K2 points of the 3x3 p3_odd_action action at even J (hc_oo_rows.py: kind EE = w0 E (x) E integrand, C3H) -> OO
    entry [-1, -1]; cache oo_p3_<kind>_J<J>_gsm.pkl (not seen by the EE-chain loaders).
Engines (spinfunc/w1_rows, p3_odd_action, p2_shadow_odd, l9_points, l5_shadow_matrix, w0_ee_common) are imported read-only and run
unchanged.  Shadows are looked up read-only in spinfunc/cache, neargap2/cache, gapscale/cache/rows, gapnear/cache/rows; new ones go to
gapseam/cache/rows.  gsm_guard refuses any write outside gapseam/cache.
Usage: OMP_NUM_THREADS=1 [GSM_DREF=5e4] [GSM_NOREF=1: state points only] python3 gsm_rows.py odd|even|oo KINDS J Delta x1,x2,...|dense
       OMP_NUM_THREADS=1 [GSM_DREF=5e4] python3 gsm_rows.py refb EE J
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gsm_guard  # noqa: F401
SF = os.path.join(HERE, '..', 'spinfunc')
NGC = os.path.join(HERE, '..', 'neargap2', 'cache')
GSR = os.path.join(HERE, '..', 'gapscale', 'cache', 'rows')
GNR = os.path.join(HERE, '..', 'gapnear', 'cache', 'rows')
ROWS = os.path.join(HERE, 'cache', 'rows')
os.makedirs(ROWS, exist_ok=True)
sys.path.insert(0, SF)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
os.chdir(SF)
import w1_rows as W1

SFC = os.path.join(SF, 'cache')
DREF = float(os.environ.get('GSM_DREF', '50000'))
DENSE = [round(0.025*i, 4) for i in range(1, 37)]          # x = 0.025 ... 0.9


def _odd_shadow_setup(P2):
    sfn = os.path.join(ROWS, 'gsm_p2shadow.pkl')
    orig_all = P2.shadow_all

    def shadow_odd(Delta, J):
        key = (float(Delta), int(J), P2.T.DEG)
        for fn in ([os.path.join(SFC, 'p2_shadow.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(SFC, 'w1_p2shadow*.pkl'))) + sorted(glob.glob(os.path.join(SFC, 'ee1000_p2shadow*.pkl'))) +
                   sorted(glob.glob(os.path.join(SFC, 'hc_p2shadow*.pkl'))) +
                   [os.path.join(NGC, 'ng_p2shadow.pkl'), os.path.join(GSR, 'gs_p2shadow.pkl'), os.path.join(GNR, 'gn_p2shadow.pkl'), sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        t0 = time.time()
        r = orig_all(float(Delta), int(J))
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        print(f"  new odd shadow ({Delta:g}, {J}): {time.time() - t0:.0f}s, resid {r['resid']}", flush=True)
        return db[key]
    P2.shadow_odd = shadow_odd


def setup_odd(oo=False):
    import p2_shadow_odd as P2
    import p3_odd_action as P
    import w0_ee_common as W0
    P.build = W0.make_build_odd_EE(P)
    if oo:
        P._fn = lambda kind, J: os.path.join(ROWS, 'oo_p3_%s_J%d_gsm.pkl' % (kind, J))
        P._old = lambda kind, J: []
    else:
        P._fn = lambda kind, J: os.path.join(ROWS, 'w1_p3_%s_J%d_gsm.pkl' % (kind, J))
        P._old = lambda kind, J: sorted(glob.glob(os.path.join(GNR, 'w1_p3_%s_J%d_gn.pkl' % (kind, J))))
    P2._DB = os.path.join(ROWS, 'gsm_p2shadow_db.pkl')
    _odd_shadow_setup(P2)
    return P


def setup_even():
    import l9_points as P
    import l5_shadow_matrix as SM
    import l1_tps as T, l5_spin as S, l1_crystal as Cr, l6_spin_action as A6
    import w0_ee_common as W0
    bEE = W0.make_build_even_EE(T, S, Cr, A6)
    orig = P.build_spin
    P.build_spin = lambda zv, nu, Delta, J, mons, sol, solt, kind: (bEE(zv, nu, Delta, J, mons, sol, solt) if kind == 'EE' else
                                                                  orig(zv, nu, Delta, J, mons, sol, solt, kind))
    P._fn = lambda kind, J: os.path.join(ROWS, 'w1_pts_%s_J%d_k2_gsm.pkl' % (kind, J))
    P._old = lambda kind, J: sorted(glob.glob(os.path.join(GNR, 'w1_pts_%s_J%d_k2_gn.pkl' % (kind, J))))
    sfn = os.path.join(ROWS, 'gsm_shadow_k2.pkl')

    def shadow(Delta, J):
        key = (float(Delta), int(J))
        for fn in ([os.path.join(SFC, 'l9_shadow_k2.pkl')] + sorted(glob.glob(os.path.join(SFC, 'q2_shadow_k2*.pkl'))) +
                   [os.path.join(SFC, 'u5_shadow_k2.pkl'), os.path.join(SFC, 'r1_shadow_k2.pkl')] +
                   sorted(glob.glob(os.path.join(SFC, 'w1_shadow_k2*.pkl'))) + [os.path.join(NGC, 'ng_shadow_k2.pkl'),
                                                                               os.path.join(GNR, 'gn_shadow_k2.pkl'), sfn]):
            if os.path.exists(fn):
                db = pickle.load(open(fn, 'rb'))
                if key in db:
                    return db[key]
        t0 = time.time()
        r = SM.shadow_matrix(Delta, J, 0.5, ntrip=6, seed=7)
        db = pickle.load(open(sfn, 'rb')) if os.path.exists(sfn) else {}
        db[key] = (r['b'], r['logK'], r['resid'])
        pickle.dump(db, open(sfn + '.tmp', 'wb')); os.replace(sfn + '.tmp', sfn)
        print(f"  new even shadow ({Delta:g}, {J}): {time.time() - t0:.0f}s, resid {r['resid']}", flush=True)
        return db[key]
    P.shadow = shadow
    return P


def existing(P, kind, J):
    R = {}
    for fn in P._old(kind, J) + [P._fn(kind, J)]:
        if os.path.exists(fn):
            R.update({(k[0], k[1]): v for k, v in pickle.load(open(fn, 'rb')).items() if len(k) == 2 or k[2] == 1.0})
    return R


if __name__ == '__main__':
    par = sys.argv[1]
    if par == 'refb':
        P = setup_even(); J = int(sys.argv[3]); mr = W1.mphys(DREF, J); t0 = time.time(); n = 0
        for kind in sys.argv[2].split(','):
            for x in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
                if not W1.have(existing(P, kind, J), DREF, x*mr):
                    t1 = time.time(); P.get(kind, J, DREF, [x*mr], verbose=False); n += 1
                    print(f"  refb {kind} J={J} D={DREF:g} x_ref={x} nu={x*mr:.1f}: {time.time() - t1:.1f}s", flush=True)
        print(f"done refb J={J} Dref={DREF:g}: {n} new points in {time.time() - t0:.0f}s", flush=True)
        sys.exit(0)
    kinds, J, D = sys.argv[2].split(','), int(sys.argv[3]), float(sys.argv[4])
    xs = DENSE if sys.argv[5] == 'dense' else [float(v) for v in sys.argv[5].split(',')]
    if par == 'oo':
        assert J % 2 == 0 and J >= 4
        P = setup_odd(oo=True)
    else:
        assert (par == 'odd') == (J % 2 == 1)
        P = setup_even() if par == 'even' else setup_odd()
    m = W1.mphys(D, J)
    t0 = time.time(); n = 0
    for kind in kinds:
        for x in xs:
            for DD in ((D,) if os.environ.get('GSM_NOREF') == '1' else (D, DREF)):
                if not W1.have(existing(P, kind, J), DD, x*m):
                    t1 = time.time()
                    P.get(kind, J, DD, [x*m], verbose=False); n += 1
                    print(f"  {kind} J={J} D={DD:g} x={x:.4f} nu={x*m:.3f}: {time.time() - t1:.1f}s", flush=True)
    print(f"done {par} {','.join(kinds)} J={J} Delta={D:.6g} Dref={DREF:g} m={m:.3f} J^2/m={J*J/m:.3f}: {n} new points in "
          f"{time.time() - t0:.0f}s", flush=True)
