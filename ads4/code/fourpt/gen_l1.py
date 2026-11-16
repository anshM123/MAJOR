"""Generate (and cache) the numerical code for the 25 word pairs of <JJJJ> (l = 1)."""
import time, sys
from spinblock2 import phi_funcs, _phi_exprs
import sympy as sp
t = time.time()
ex = _phi_exprs(1)
print("symbolic word pairs built: %.0f s" % (time.time() - t), flush=True)
for k, v in ex.items():
    print(k, len(sp.Add.make_args(v)), flush=True)
phi_funcs(1)
print("code generated: %.0f s" % (time.time() - t), flush=True)
