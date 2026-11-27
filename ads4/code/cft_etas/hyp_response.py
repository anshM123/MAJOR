"""Hyperbolic responses from cache_hyp_quad.pkl (hyp_fast.py): reduced Lagrangians after the H^2 integral, membrane eta,
Wald entropy, and delta(x) = e1 - sigma1 and V(x) = delta + (1-x) delta' for theta(r) I and theta(r) J.
Planar check: x -> 1 must give dV_E = -6 theta_I'(r_h) r_h  and  -108 theta_J + 18 theta_J' r_h."""
import pickle, sympy as sp
import hyp_engine  # sets sys.path
import static_engine as se_mod
from static_engine import reduce_quadratic, P
from hyp_engine import r, rh, u, integrate_u, flux_at_horizon
data = pickle.load(open('cache_hyp_quad.pkl', 'rb'))
loc = {'r': r, 'r_h': rh, 'u': u}
mu = rh**3 - rh
f = r**2 - 1 - mu/r

def reduced(name, weight=1):
    L2 = 0
    for (a, b), s in data[name]['quad'].items():
        q = sp.sympify(s, locals=loc)
        L2 += integrate_u(sp.together(q*r**2))*weight*P[a]*P[b]
    return L2

print('background R, I, J:', [sp.factor(sp.sympify(data[k]['c0'], locals=loc)) for k in ('R', 'I', 'J')])
print('first-order R, I, J zero:', [data[k]['c1'] == '0' for k in ('R', 'I', 'J')])
alE = reduce_quadratic(reduced('R'))
print('Einstein alphas:', {k: sp.factor(v) for k, v in alE.items()})
eta0 = flux_at_horizon(alE, 1/f)
print('eta0 =', sp.factor(eta0))
th = sp.Function('theta')(r)
Tj = sp.symbols('t0:5')
thpoly = sum(Tj[k]*(r - rh)**k/sp.factorial(k) for k in range(5))
e1 = {}
for name in ('I', 'J'):
    al = reduce_quadratic(sp.expand(reduced(name)*th))
    al = {k: sp.simplify(v.subs(th, thpoly).doit()) for k, v in al.items()}
    if 0 in al:
        print('[%s] mass term alpha_0 =' % name, sp.factor(al[0]))
    e1[name] = sp.simplify(flux_at_horizon(al, 1/f)/eta0)
    print('[theta %s] e1 = ' % name, sp.collect(sp.expand(e1[name]), Tj))
# Wald on the background
t, Y = sp.symbols('t Y')
w = 1 - u**2
gB = sp.diag(-f, 1/f, r**2/w**2, r**2/w)
gi, Gam, Rl, Ric, Rs, W = se_mod.curvature(gB)
dd = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: Rl[k] + dd*Xt.get(k, 0) for k in Rl}
Wp, Rsp = se_mod.weyl_from(Rp, gB, gi)
Ip, Jp = se_mod.invariants(gB, gi, Wp)
dR = sp.diff(Rsp, dd).subs(dd, 0)
sig = {}
for name, inv in (('I', Ip), ('J', Jp)):
    ratio = sp.factor(sp.simplify((sp.diff(inv, dd).subs(dd, 0)/dR).subs(r, rh)))
    sig[name] = ratio*Tj[0]
    print('[theta %s] sigma1 = s1/s0 = t0 *' % name, ratio)
x = sp.Symbol('x')
for name in ('I', 'J'):
    dl = sp.expand(sp.simplify(e1[name] - sig[name]))
    print('[theta %s] delta = e1 - sigma1 =' % name, sp.collect(dl, Tj))
    print('    in terms of x (rh^2 = 1/(1-x)) and theta jets t1*rh etc.:',
          sp.collect(sp.simplify(sp.expand(dl.subs({Tj[1]: Tj[1]/rh, Tj[2]: Tj[2]/rh**2}).subs(rh, 1/sp.sqrt(1 - x)))), Tj))
pickle.dump({'e1': {k: str(v) for k, v in e1.items()}, 'sig': {k: str(v) for k, v in sig.items()}}, open('cache_hyp_resp.pkl', 'wb'))
