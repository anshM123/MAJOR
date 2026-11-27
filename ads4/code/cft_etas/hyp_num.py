"""Hyperbolic shear response at a FIXED rational horizon radius r_h (argv[1]), sympy with cancellation.
Outputs e1 (= eta1/eta0) and sigma1 (= s1/s0) for theta(r) I and theta(r) J in terms of horizon jets t0..t2 of theta.
x = 1 - 1/r_h^2.  Conventions: hyp_engine.py.  Appends a line to out_hyp_num.txt and caches cache_hyp_num_<rh>.pkl."""
import sys, time, itertools, pickle
import sympy as sp
import hyp_engine as he
import static_engine as se_mod
from static_engine import reduce_quadratic, jets, P, phi, r, eps
RH = sp.Rational(sys.argv[1])
T0 = time.time()
u = he.u; w = 1 - u**2
mu = RH**3 - RH
f = sp.factor(r**2 - 1 - mu/r)
g = he.hyp_metric(f).subs(he.rh, RH)
gi, Gam, Rl, Ric, Rs, W = se_mod.curvature(g)
print('curvature %.0fs' % (time.time() - T0)); sys.stdout.flush()
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def ser(e):
    e = jets(e)
    return [sp.cancel(sp.diff(e, eps, k).subs(eps, 0)/sp.factorial(k)) for k in range(3)]
Mk = [[None]*6 for _ in range(6)]
for i, (a, b) in enumerate(PAIRS):
    for j, (c, d) in enumerate(PAIRS):
        Mk[i][j] = ser(W[a, b, c, d]*gi[c, c]*gi[d, d])
print('bivector series %.0fs' % (time.time() - T0)); sys.stdout.flush()
off = [(i, j) for i in range(6) for j in range(6) if i != j and Mk[i][j][0] != 0]
print('background off-diagonal entries:', off)
assert not off
d0 = [Mk[i][i][0] for i in range(6)]
M11d = [sp.cancel(sum(Mk[i][k][1]*Mk[k][i][1] for k in range(6))) for i in range(6)]
I2 = sp.cancel(4*(2*sum(d0[i]*Mk[i][i][2] for i in range(6)) + sum(M11d)))
J2 = sp.cancel(24*(sum(d0[i]**2*Mk[i][i][2] for i in range(6)) + sum(d0[i]*M11d[i] for i in range(6))))
I0 = sp.factor(4*sum(x**2 for x in d0)); J0 = sp.factor(8*sum(x**3 for x in d0))
I1 = sp.cancel(8*sum(d0[i]*Mk[i][i][1] for i in range(6)))
print('I0 =', I0, ' J0 =', J0, ' I1 =', I1); sys.stdout.flush()
R2 = ser(Rs)[2]
def reduced(dens):
    e = sp.expand(dens*r**2)
    out = 0
    for a in range(5):
        for b in range(a, 5):
            if a == b:
                c = e.coeff(P[a], 2)
                for o in range(5):
                    if o != a: c = c.subs(P[o], 0)
            else:
                c = e.coeff(P[a], 1).coeff(P[b], 1)
            if c != 0:
                out += he.integrate_u(c)*P[a]*P[b]
    return out
alE = reduce_quadratic(reduced(R2))
flux = lambda al: sp.limit(sp.simplify(-sum((-1)**(k - 1)*sp.diff(2*a*sp.diff(-1/f, r, k - 1), r, k - 1) for k, a in al.items() if k > 0)), r, RH)
eta0 = flux(alE)
print('Einstein alphas', alE, ' eta0 =', eta0); sys.stdout.flush()
th = sp.Function('theta')(r)
Tj = sp.symbols('t0:5')
thpoly = sum(Tj[k]*(r - RH)**k/sp.factorial(k) for k in range(5))
res = {'rh': str(RH), 'I0': str(I0), 'J0': str(J0)}
for name, d2 in (('I', I2), ('J', J2)):
    al = reduce_quadratic(sp.expand(reduced(d2)*th))
    al = {k: sp.simplify(v.subs(th, thpoly).doit()) for k, v in al.items()}
    if 0 in al:
        print('[%s] alpha_0 =' % name, sp.factor(al[0]))
    e1 = sp.expand(sp.simplify(flux(al)/eta0))
    res['e1' + name] = str(e1)
    print('[theta %s] e1 =' % name, e1); sys.stdout.flush()
# Wald on the background
gB = g.subs(eps, 0)
giB, GamB, RlB, RicB, RsB, WB = se_mod.curvature(gB)
dd = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
Rp = {k: RlB[k] + dd*Xt.get(k, 0) for k in RlB}
Wp, Rsp = se_mod.weyl_from(Rp, gB, giB)
Ip, Jp = se_mod.invariants(gB, giB, Wp)
dR = sp.diff(Rsp, dd).subs(dd, 0)
for name, inv in (('I', Ip), ('J', Jp)):
    ratio = sp.simplify((sp.diff(inv, dd).subs(dd, 0)/dR).subs(r, RH))
    res['s' + name] = str(ratio)
    print('[theta %s] sigma1/t0 =' % name, ratio)
pickle.dump(res, open('cache_hyp_num_%s.pkl' % str(RH).replace('/', '_'), 'wb'))
print('done %.0fs' % (time.time() - T0))
