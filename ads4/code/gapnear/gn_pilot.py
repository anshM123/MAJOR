"""
gn_pilot (GAP_NEARGAP.md §4): one batch-4 job run locally, exactly as ext/run_ext.py runs it (ext/point.py, unchanged, executed
in-process), with gn_guard installed first and the output directory under gapnear/cache/pilot/.  Validates the n1 settings at
Delta = 2000 + J (cost per node, replica errors, failure at small nu / near the light-cone edge).
Usage: python3 gn_pilot.py METH J DELTA NU SEED TAG
"""
import os, sys, runpy
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gn_guard  # noqa: F401
meth, J, D, nu, seed, tag = sys.argv[1:7]
out = os.path.join(HERE, 'cache', 'pilot', tag, f'{meth}_J{J}_D{float(D):g}_nu{float(nu):.6f}_s{seed}')
sys.argv = [os.path.join(HERE, '..', '..', 'ext', 'point.py'), meth, J, D, nu, seed, out]
runpy.run_path(sys.argv[0], run_name='__main__')
