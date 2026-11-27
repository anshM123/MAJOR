"""Wald entropy ratios s1/s0 for theta*I and theta*J on the hyperbolic black hole (symbolic r_h), background only."""
import sympy as sp
import hyp_engine as he
import static_engine as se_mod
from static_engine import r
rh = he.rh; u = he.u; w = 1 - u**2
f = r**2 - 1 - (rh**3 - rh)/r
gB = sp.diag(-f, 1/f, r**2/w**2, r**2/w)
gi, Gam, Rl, Ric, Rs, W = se_mod.curvature(gB)
dd = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: Rl[k] + dd*Xt.get(k, 0) for k in Rl}
Wp, Rsp = se_mod.weyl_from(Rp, gB, gi)
Ip, Jp = se_mod.invariants(gB, gi, Wp)
dR = sp.simplify(sp.diff(Rsp, dd).subs(dd, 0))
x = sp.Symbol('x')
for name, inv in (('I', Ip), ('J', Jp)):
    ratio = sp.simplify(sp.diff(inv, dd).subs(dd, 0)/dR)
    val = sp.factor(sp.simplify(ratio.subs(r, rh)))
    print('[theta %s] s1/s0 = theta(r_h) *' % name, val, '  = in x:', sp.factor(sp.simplify(val.subs(rh, 1/sp.sqrt(1 - x)))))
print('background I, J at horizon in x:', sp.factor(sp.simplify(se_mod.invariants(gB, gi, W)[0].subs(r, rh).subs(rh, 1/sp.sqrt(1-x)))),
      sp.factor(sp.simplify(se_mod.invariants(gB, gi, W)[1].subs(r, rh).subs(rh, 1/sp.sqrt(1-x)))))
