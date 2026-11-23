"""
ext/point.py -- one direct crystal point (or one small-nu node), for running on an external machine.
Identical physics to ads4/code/smallnu/sn_points.py and ads4/code/highspin/hs_points.py: it calls
ads4/code/spinfunc/hc_points.py unchanged.  Each job writes into its own cache directory (no shared files).
Usage: python3 point.py METH J DELTA NU SEED OUTDIR     (METH = n1 | n3; NU explicit, not x)
"""
import os, sys, time
sys.dont_write_bytecode = True
meth, J, D, nu, seed, out = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'code', 'spinfunc')
os.makedirs(out, exist_ok=True)
os.environ['HC_CACHE'] = os.path.abspath(out)
os.environ['HC_METHOD'] = meth
os.environ['Z2_SEED'] = str(seed)
defaults = ((('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')) if meth == 'n1'
            else (('HC_M', '10'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')))
for k, v in defaults:
    os.environ.setdefault(k, v)
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
sys.path.insert(0, SF)
os.chdir(SF)
import warnings
warnings.filterwarnings('ignore')
import hc_points as HP
HP.STATE[0] = D
m = HP.mphys(D, J)
HP.setup(odd=(J % 2 == 1))
t0 = time.time()
print(f"# point [{meth}] J={J} Delta={D:g} m={m:.3f} nu={nu} x={nu/m:.6f} seed={seed}", flush=True)
if meth == 'n1':
    HP.run_n1(J, [nu], D, HP.Z2.MQ, HP.Z2.RQ)
else:
    HP.Z2.run(J, [nu], [D, 20000.0], HP.Z2.MQ, HP.Z2.RQ)
print(f"# done [{time.time() - t0:.0f}s]", flush=True)
