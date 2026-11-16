"""Generate (and cache) the numerical code for the 25 word pairs of <JJJJ> (l = 1): coefficients of the z-monomials."""
import time
import sympy as sp
from spinblock2 import phi_funcs, _phi_exprs, _monomials
t = time.time()
ex = _phi_exprs(1)
print("symbolic word pairs: %.0f s" % (time.time() - t), {k: len(sp.Add.make_args(v)) for k, v in ex.items()}, flush=True)
monos, polys = _monomials(1)
print("monomials: %d  (%.0f s)" % (len(monos), time.time() - t), flush=True)
phi_funcs(1)
print("code generated: %.0f s" % (time.time() - t), flush=True)
