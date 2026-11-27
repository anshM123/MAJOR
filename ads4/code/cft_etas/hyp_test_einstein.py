"""Step 1: Einstein TT shear mode on the hyperbolic black hole is a massless scalar (no alpha_0), and eta/s = 1/4pi."""
import time, sympy as sp
from hyp_engine import *
t0 = time.time()
f = background_f()
g = hyp_metric(f)
gi, Gam, Rl, Ric, Rs, W = se.curvature(g)
print('curvature done %.1fs' % (time.time() - t0))
Rss = eps_series(Rs + 6)   # Lagrangian R + 6
print('background R+6 =', sp.simplify(Rss[0]), ' (Einstein: R = -12 => R+6 = -6)')
L2 = reduced_lagrangian(Rss[2])
print('reduced Einstein L2 =', sp.factor(L2))
al = reduce_quadratic(L2)
print('alphas:', {k: sp.factor(v) for k, v in al.items()})
eta0 = flux_at_horizon(al, 1/f)
print('eta0 (membrane flux) =', sp.factor(eta0), ' ; per horizon area r_h^2:', sp.simplify(eta0/rh**2))
I, J = se.invariants(g, gi, W)
Is, Js = eps_series(I), eps_series(J)
print('background I, J:', sp.factor(Is[0]), sp.factor(Js[0]), ' first order:', sp.simplify(Is[1]), sp.simplify(Js[1]))
print('done %.1fs' % (time.time() - t0))
