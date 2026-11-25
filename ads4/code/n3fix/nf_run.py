"""n3fix/nf_run.py: run ext/point_patched.py (the exact batch-6 code path) locally under the n3fix write guard, output in
n3fix/cache/<sub>.  Usage: python3 nf_run.py SUB METH J DELTA x1,x2,.. [SEED]   (x = nu/m(DELTA, J));  env HC_M etc. as point.py
NF_ORIG=1: the ORIGINAL ext/point.py (unpatched) instead."""
import nf_guard  # noqa: F401
import os, sys, runpy
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path = [p for p in sys.path if os.path.abspath(p or '.') != HERE]
EXT = os.path.abspath(os.path.join(HERE, '..', '..', 'ext'))
sub, meth, J, D, xs = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), sys.argv[5]
seed = sys.argv[6] if len(sys.argv) > 6 else '4242'
al = D - J - 2.0; m = np.sqrt(al*(al + 2*J + 1))
nus = [float(v) for v in xs[3:].split(',')] if xs.startswith('nu=') else [float(x)*m for x in xs.split(',')]
out = os.path.join(HERE, 'cache', sub)
if os.environ.get('NF_ORIG') == '1':
    for nu in nus:
        sys.argv = [os.path.join(EXT, 'point.py'), meth, str(J), repr(D), repr(float(nu)), seed, out]
        runpy.run_path(sys.argv[0], run_name='__main__')
else:
    sys.argv = [os.path.join(EXT, 'point_patched.py'), meth, str(J), repr(D), ','.join(repr(float(v)) for v in nus), seed, out]
    runpy.run_path(sys.argv[0], run_name='__main__')
