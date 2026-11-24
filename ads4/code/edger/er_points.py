"""
edger/er_points.py (EDGE_R.md §2, §6): DIRECT crystal points with spinfunc/hc_points.py UNCHANGED (imported read-only, run in-process)
under the er_guard write guard, cache moved to edger/cache (hc_points' own HC_CACHE switch).
Method (env ER_METHOD, default n3): 'n3' = the n3 hybrid wall-free contour (Phi-thimble flow beyond rho1 = 3.5, graviton-adapted proposal,
state and Delta_ref = 2e4 reference on common random numbers; z2_points.n3_pair) -- the method that survives nu/Delta > 0.94 (EDGE_R §2);
'n1' = the adaptive re-centred real contour (HIGHSPIN's), which fails there.
x-list as hs_points: comma list of x = nu/m(Delta, J) and/or 'gk15:a:b[:xcap]' (15 GK nodes), joined by '+'.
Env: ER_M (log2 points per replica, default 10), ER_R (replicas, 8), ER_SEED (n3 base seed; independent seeds must differ by >= 8).
Usage: python3 er_points.py J Delta x-list
"""
import er_guard  # noqa: F401  (first)
import os, sys, runpy
HERE = os.path.dirname(os.path.abspath(__file__))
SF = os.path.abspath(os.path.join(HERE, '..', 'spinfunc'))
sys.path.insert(0, SF); sys.path.insert(0, os.path.join(HERE, '..', 'highspin'))
sys.dont_write_bytecode = True

if __name__ == '__main__':
    J, D, spec = int(sys.argv[1]), float(sys.argv[2]), sys.argv[3]
    import hs_points as HSP                     # parse_x only (read-only; its main is not run)
    import numpy as np
    al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
    nus = [x*m for x in HSP.parse_x(spec)]
    meth = os.environ.get('ER_METHOD', 'n3')
    os.environ['HC_CACHE'] = os.path.join(HERE, 'cache')
    os.environ['HC_METHOD'] = meth
    os.environ['HC_M'] = os.environ.get('ER_M', '10'); os.environ['HC_MREF'] = os.environ.get('ER_MREF', os.environ['HC_M'])
    os.environ['HC_R'] = os.environ.get('ER_R', '8')
    os.environ.setdefault('Z2_SEED', os.environ.get('ER_SEED', '4242'))
    os.environ.setdefault('Z2_NUSPLIT', '5'); os.environ.setdefault('HC_ADAPT', '1')
    os.chdir(SF)
    sys.argv = [os.path.join(SF, 'hc_points.py'), str(J), repr(D), 'nu=' + ','.join(repr(float(v)) for v in nus)]
    runpy.run_path(sys.argv[0], run_name='__main__')
