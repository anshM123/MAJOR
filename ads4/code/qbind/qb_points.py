"""
qb_points (QBIND_COROLLARY.md §1; copy of neargap2/ng_points.py, cache -> qbind/cache): DIRECT crystal points (hc_points.py, unchanged, imported read-only; n1 adaptive re-centred contour, RQMC,
state and reference Delta_ref = 2e4 on common random numbers; odd J -> parity-odd 1x1 integrand, even J -> w0 2x2 EE integrand + C3H),
written to qbind/cache/hc_pts_<kind>_J<J>_D<D>.pkl through hc_points' own HC_CACHE switch (no other change).
Settings of RESULT_D_CLEANUP §2 (edge nodes): HC_METHOD=n1 HC_ADAPT=1 HC_M=11 HC_MREF=11.
Usage: python3 ng_points.py J Delta x1,x2,...
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.join(HERE, '..', 'spinfunc')
sys.path.insert(0, SF); sys.path.insert(0, HERE)
import qb_guard  # noqa: F401
os.chdir(SF)
os.environ['HC_CACHE'] = os.path.join(HERE, 'cache')
for k, v in (('HC_METHOD', 'n1'), ('HC_ADAPT', '1'), ('HC_M', '11'), ('HC_MREF', '11'), ('L1_DEG', '6'), ('Z2_NUSPLIT', '5')):
    os.environ.setdefault(k, v)
import warnings
warnings.filterwarnings('ignore')
import hc_points as HP

if __name__ == '__main__':
    J = int(sys.argv[1]); D = float(sys.argv[2])
    HP.STATE[0] = D
    m = HP.mphys(D, J)
    nus = [float(v)*m for v in sys.argv[3].split(',')]
    HP.setup(odd=(J % 2 == 1))
    t0 = time.time()
    print(f"# qb_points [n1] J={J} Delta={D:g} m={m:.3f} x={sys.argv[3]} -> {HP.fn('EE', J)}", flush=True)
    HP.run_n1(J, nus, D, HP.Z2.MQ, HP.Z2.RQ)
    print(f"# done [{time.time() - t0:.0f}s]", flush=True)
