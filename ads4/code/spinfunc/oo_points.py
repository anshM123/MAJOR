"""
oo_points (HYPOTHESES_CLOSE2.md, H_OO at J >= 12): DIRECT small-nu crystal points of the FULL 3x3 [even(2), odd(1)] structure action at
EVEN J -- the integrand of the OO block (X^O_same = entry [-1, -1]) -- state and heavy reference (Delta_ref = 2e4, same nu) on COMMON
RANDOM NUMBERS.  Nothing new in the method: hc_points' n3 path (z2 n3 hybrid contour, nu_f = max(nu, 15) above 13, nu_f = 13 below;
proposal adapted on the state and re-used for the reference, same Sobol seeds; Q_J overflow rescaling; Legendre recurrence for the
scalar Q_J; Fourier floor fw_frame.patch_S), with hc_points.setup(odd=True) FORCED at even J, i.e. the graviton evaluators are
n3_run_odd.OddEval on p8_odd_pw's point-wise copy of p3_odd_action.build:
  kind 'EE'  : w0_ee_common.make_build_odd_EE (y5_n4.patch) -- C^'_EE integrand, 3x3 pref (2 parity-even + 1 parity-odd structure);
  kind 'C3H' : p3's own o_H q integrand (C^_3), 3x3.
The scalar (kind 'S') is the unchanged even-J scalar on the same samples (so the scalar control variate of z3 applies).
Caches (ONE FILE PER STATE; the reference points of a state live only in that state's file):
  cache/oo_pts_<kind>_J<J>_D<Dstate>.pkl  keyed (Delta, round(nu, 9)), z2/n1_run format (E 3x3, err, reps, logN, Es, errs, repsS, logNs);
  raw: cache/oo_raw_J<J>_D<Dstate>.pkl.   Resumable (finished nodes skipped).
Usage: OMP_NUM_THREADS=1 [OO_M=10 OO_R=8] python3 oo_points.py J Dstate x1,x2,...   (x = nu/m(Dstate, J))  |  ... J Dstate nu=n1,n2,...
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ['HC_M'] = os.environ.get('OO_M', '10')
os.environ['HC_R'] = os.environ.get('OO_R', '8')
os.environ.setdefault('Z2_NUSPLIT', '5')          # n3 also at nu = 8.3, 10 (no thimbles)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import hc_points as HP
Z2 = HP.Z2

C = os.environ.get('OO_CACHE', os.path.join(HERE, 'cache'))
STATE = [None]


def fn(kind, J):
    return os.path.join(C, 'oo_pts_%s_J%d_D%g.pkl' % (kind, J, STATE[0]))


def rawfn(J):
    return os.path.join(C, 'oo_raw_J%d_D%g.pkl' % (J, STATE[0]))


def _save(path, key, val):
    if os.path.basename(path).startswith('z2') and path.endswith('_raw.pkl'):
        path = rawfn(key[2])
    return HP._save0(path, key, val)


Z2.fn = fn
Z2._save = _save


if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2])
    assert J % 2 == 0 and J >= 4
    STATE[0] = D
    m = HP.mphys(D, J)
    a = sys.argv[3]
    nus = [float(v) for v in a[3:].split(',')] if a.startswith('nu=') else [float(v)*m for v in a.split(',')]
    HP.setup(odd=True)                 # FORCED at even J: OddEval (p8 point-wise p3 build, 3x3) for the graviton kinds
    t0 = time.time()
    print(f"# oo_points [n3, 3x3 p8 integrand] J={J} Delta={D:g} m={m:.3f}: nu = {[round(v, 4) for v in nus]} "
          f"(x = {[round(v/m, 5) for v in nus]}); M={Z2.MQ} R={Z2.RQ} kinds {Z2.KINDS}", flush=True)
    Z2.run(J, nus, [D, 20000.0], Z2.MQ, Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
