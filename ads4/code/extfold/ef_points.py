"""
ef_points (EXT_FOLD.md §4.1): recompute one ext job point LOCALLY (spinfunc/hc_points.py, unchanged, read-only; exactly ext/point.py's
settings) into extfold/cache/repro/<tag>/, to check that the external machine reproduces this machine's numbers.
Usage: python3 ef_points.py METH J DELTA NU SEED TAG
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_guard  # noqa: F401  (first)
meth, J, D, nu, seed, tag = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), int(sys.argv[5]), sys.argv[6]
SF = os.path.join(HERE, '..', 'spinfunc')
out = os.path.join(HERE, 'cache', 'repro', tag)
os.makedirs(out, exist_ok=True)
os.environ['HC_CACHE'] = os.path.abspath(out)
os.environ['HC_METHOD'] = meth
os.environ['Z2_SEED'] = str(seed)
defaults = ((('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')) if meth == 'n1'
            else (('HC_M', '10'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')))
for k, v in defaults:
    os.environ.setdefault(k, v)
sys.path.insert(0, SF)
os.chdir(SF)
import warnings
warnings.filterwarnings('ignore')
import hc_points as HP
HP.STATE[0] = D
m = HP.mphys(D, J)
HP.setup(odd=(J % 2 == 1))
t0 = time.time()
print(f"# ef_points [{meth}] J={J} Delta={D:g} m={m:.3f} nu={nu} x={nu/m:.6f} seed={seed} -> {out}", flush=True)
if meth == 'n1':
    HP.run_n1(J, [nu], D, HP.Z2.MQ, HP.Z2.RQ)
else:
    HP.Z2.run(J, [nu], [D, 20000.0], HP.Z2.MQ, HP.Z2.RQ)
print(f"# done [{time.time() - t0:.0f}s]", flush=True)
