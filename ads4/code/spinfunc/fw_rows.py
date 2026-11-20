"""
fw_rows (FOLD_WINDOW_LARGEJ.md): K2 (L1_DEG = 6) crystal records with the fw_frame precision fixes (spin-J Fourier floor; scalar Q_J by
recurrence), for (a) helicity bases at large J (state rows and the heavy reference 'refb' points x_ref = 0.3 ... 0.8), (b) the J = 48 rows
at (1048, 48) recomputed.  Same engine as ee1000_rows (l9_points + w0 EE integrand), same shadow lookup (existing caches; new ones ->
cache/fw_shadowk2.pkl).  Cache: cache/fw_k2_<kind>_J<J>.pkl keyed (Delta, round(nu, 9)).  Resumable.
Usage: OMP_NUM_THREADS=1 python3 fw_rows.py row|state KINDS J Delta xspec (xspec: X35 | B (basis set 0.15..0.8 step 0.05) | x1,x2,..;
                                                                          'state': no reference points at the same nu)
       OMP_NUM_THREADS=1 python3 fw_rows.py refb KINDS J1,J2,...           (reference Delta = 2e4 at x_ref = 0.3, ..., 0.8)
"""
import os, sys, glob, pickle, time
os.environ['L1_DEG'] = '6'
os.environ['EE1000_Q'] = 'fw'          # new shadows -> cache/ee1000_shadow_k2_fw.pkl (never another agent's file)
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

C = os.path.join(HERE, 'cache')
DREF = 20000.0
XB = [round(0.15 + 0.05*i, 4) for i in range(14)]


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def fn(kind, J):
    return os.path.join(C, 'fw_k2_%s_J%d.pkl' % (kind, J))


def load(kind, J):
    f = fn(kind, J)
    return pickle.load(open(f, 'rb')) if os.path.exists(f) else {}


def setup():
    import ee1000_rows as ER
    P = ER.setup_even()                          # EE build + ee1000 shadow lookup
    import fw_frame as FW
    import l1_crystal as Cr
    FW.patch_qrec_tps(Cr, HERE)
    FW.install_tps(P)
    FW.patch_conserved_clean()
    P._fn = fn
    sfn = os.path.join(C, 'fw_shadowk2.pkl')
    sh0 = P.shadow

    def shadow(Delta, J):
        if os.path.exists(sfn):
            db = pickle.load(open(sfn, 'rb'))
            if (float(Delta), int(J)) in db:
                return db[(float(Delta), int(J))]
        return sh0(Delta, J)
    P.shadow = shadow
    return P


def compute(P, kind, J, D, xs):
    m = mphys(D, J)
    for x in xs:
        nu = x*m
        if (float(D), round(nu, 9)) in load(kind, J):
            continue
        t = time.time()
        P.get(kind, J, D, [nu], verbose=False)
        print(f"  fw_rows {kind} J={J} Delta={D:g} x={x}: {time.time() - t:.0f}s", flush=True)


if __name__ == '__main__':
    P = setup()
    mode, kinds = sys.argv[1], sys.argv[2].split(',')
    if mode == 'refb':
        for J in [int(a) for a in sys.argv[3].split(',')]:
            for kind in kinds:
                compute(P, kind, J, DREF, [0.3, 0.4, 0.5, 0.6, 0.7, 0.8])
        sys.exit(0)
    J, D, xa = int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
    import ee1000_rows as ER
    xs = ER.XSETS[xa] if xa in ER.XSETS else (XB if xa == 'B' else [float(v) for v in xa.split(',')])
    m = mphys(D, J)
    for x in xs:                                  # state and reference at the same nu, node by node (resumable)
        for kind in kinds:
            compute(P, kind, J, D, [x])
            if mode == 'row':
                compute(P, kind, J, DREF, [x*m/mphys(DREF, J)])
    print(f"done fw_rows {kinds} J={J} Delta={D:g}", flush=True)
