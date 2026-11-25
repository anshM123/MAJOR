"""
b4_points (BATCH4_FOLD.md §1.2; a copy of b5fold/b5_points.py): LOCAL re-computation of batch-4 direct points -- ext/point.py's code path (spinfunc/hc_points.py
unchanged, imported read-only; n1, HC_ADAPT=1, HC_M=HC_MREF=11, L1_DEG=6, Z2_NUSPLIT=5, 8 replicas, fixed Sobol seeds 2000 + r), with
the cache moved to b4fold/cache/loc/ and only the kinds in B4_KINDS (default EE).  Purpose: the even-J structure basis.  hc_points'
even records are matrices over the <TTO> conserved-structure basis sol(Delta) x solt(3 - Delta) of l5_spin.conserved, an SVD basis of
a 2-dim null space, which is fixed only up to a 2x2 transformation (rotation in the degenerate subspace + unit-max column scaling); the
external machine's LAPACK may pick another one (EDGE_R2 found this on batch 3).  The n1 integral itself is deterministic, so one local
point at the same (Delta, nu) gives the transformation exactly (b4_basis.py).
Usage: python3 b4_points.py J Delta nu1,nu2,...   (explicit nu, as in ext/jobs_batch4.txt)
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b4_guard  # noqa: F401  (first)
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
OUT = os.path.join(HERE, 'cache', 'loc')
os.makedirs(OUT, exist_ok=True)

if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2]); nus = [float(v) for v in sys.argv[3].split(',')]
    os.environ['HC_CACHE'] = OUT
    os.environ['HC_METHOD'] = 'n1'
    os.environ['Z2_SEED'] = '4242'
    os.environ['Z2_KINDS'] = os.environ.get('B4_KINDS', 'EE')
    for k, v in (('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')):
        os.environ.setdefault(k, v)
    sys.path.insert(0, SF)
    os.chdir(SF)
    import warnings
    warnings.filterwarnings('ignore')
    import hc_points as HP
    HP.STATE[0] = D
    HP.setup(odd=(J % 2 == 1))
    t0 = time.time()
    for nu in nus:
        HP.run_n1(J, [nu], D, HP.Z2.MQ, HP.Z2.RQ)
    print(f"# done J={J} Delta={D:g} nus={nus} [{time.time() - t0:.0f}s]", flush=True)
