"""
Hyperbolic shear response for theta(r) I and theta(r) J (fast engine: eps-truncated diagonal curvature + bivector traces).
See hyp_engine.py for conventions.  Output: e1 = eta1/eta0, sigma1 = s1/s0 as functions of the horizon jets of theta,
for the hyperbolic black hole with horizon r_h; x = 1 - 1/r_h^2.  Planar check: r_h -> infinity must reproduce
e1 = 4 t0 - 6 t1, sigma1 = 4 t0 (theta I) and e1 = -102 t0 + 18 t1, sigma1 = 6 t0 (theta J) with t1 in units r_h = 1.
"""
import itertools, time, sys
import sympy as sp
from hyp_engine import r, eps, P, phi, rh, w, u, reduced_lagrangian, flux_at_horizon, background_f
from static_engine import reduce_quadratic, jets

T0 = time.time()
t, Y = sp.symbols('t Y')
Xc = [t, r, u, Y]
ORD = 2


def trunc(e):
    e = sp.expand(e)
    return sum(e.coeff(eps, k)*eps**k for k in range(ORD + 1))


A = phi*w
E = 1 + eps*A + eps**2*A**2/2
Ei = 1 - eps*A + eps**2*A**2/2
f = background_f()
g = [-f, 1/f, r**2*E/w**2, r**2*Ei/w]
gi = [-1/f, f, trunc(w**2/r**2*Ei), trunc(w*E/r**2)]


def d(e, c):
    return sp.diff(e, Xc[c])

# Christoffels (diagonal metric): Gam^a_bc = gi_a/2 (d_c g_ab + d_b g_ac - d_a g_bc)
Gam = {}
for a, b, c in itertools.product(range(4), repeat=3):
    v = 0
    if a == b: v += d(g[a], c)
    if a == c: v += d(g[a], b)
    if b == c: v -= d(g[b], a)
    Gam[a, b, c] = trunc(gi[a]*v/2) if v != 0 else 0
print('Gamma done %.1fs' % (time.time() - T0)); sys.stdout.flush()
# Riemann R^a_bcd then lower
Rl = {}
for a, b, c, dd in itertools.product(range(4), repeat=4):
    if c >= dd:
        continue
    v = d(Gam[a, b, dd], c) - d(Gam[a, b, c], dd) + sum(Gam[a, c, h]*Gam[h, b, dd] - Gam[a, dd, h]*Gam[h, b, c] for h in range(4))
    v = trunc(g[a]*v)
    Rl[a, b, c, dd] = v
    Rl[a, b, dd, c] = -v
for a, b, c in itertools.product(range(4), repeat=3):
    Rl[a, b, c, c] = 0
print('Riemann done %.1fs' % (time.time() - T0)); sys.stdout.flush()
Ric = [[trunc(sum(gi[a]*Rl[a, b, a, dd] for a in range(4))) for dd in range(4)] for b in range(4)]
Rs = trunc(sum(gi[b]*Ric[b][b] for b in range(4)))
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def Wl(a, b, c, dd):
    gm = lambda i, j: g[i] if i == j else 0
    return (Rl[a, b, c, dd] - sp.Rational(1, 2)*(gm(a, c)*Ric[b][dd] - gm(a, dd)*Ric[b][c] - gm(b, c)*Ric[a][dd] + gm(b, dd)*Ric[a][c])
            + Rs/6*(gm(a, c)*gm(b, dd) - gm(a, dd)*gm(b, c)))

# bivector matrix M_AB = W_{ab}^{cd}, A=(a<b), B=(c<d) ; I = 4 tr M^2, J = 8 tr M^3
M = sp.zeros(6, 6)
for i, (a, b) in enumerate(PAIRS):
    for j, (c, dd) in enumerate(PAIRS):
        M[i, j] = trunc(Wl(a, b, c, dd)*gi[c]*gi[dd])
print('Weyl bivector done %.1fs' % (time.time() - T0)); sys.stdout.flush()
M2 = (M*M).applyfunc(trunc)
I = trunc(4*M2.trace())
J = trunc(8*(M2*M).trace())
print('I, J done %.1fs' % (time.time() - T0)); sys.stdout.flush()


def coeffs(e):
    e = jets(sp.expand(e))
    return [sp.simplify(e.coeff(eps, k)) for k in range(ORD + 1)]

Ic, Jc, Rc = coeffs(I), coeffs(J), coeffs(Rs)
print('background I =', sp.factor(Ic[0]), ' J =', sp.factor(Jc[0]), ' R =', sp.factor(Rc[0]))
print('first order I, J (must vanish):', sp.simplify(Ic[1]), sp.simplify(Jc[1]))
sys.stdout.flush()

# Einstein normalisation
alE = reduce_quadratic(reduced_lagrangian(Rc[2]))
eta0 = flux_at_horizon(alE, 1/f)
print('Einstein alphas', {k: sp.factor(v) for k, v in alE.items()}, ' eta0 =', sp.factor(eta0))

# theta-weighted operators with local horizon jets
th = sp.Function('theta')(r)
Tj = sp.symbols('t0:5')
thpoly = sum(Tj[k]*(r - rh)**k/sp.factorial(k) for k in range(5))
res = {}
for name, ser in (('I', Ic), ('J', Jc)):
    L2 = reduced_lagrangian(sp.expand(th*ser[2]))
    al = reduce_quadratic(L2)
    al = {k: sp.simplify(v.subs(th, thpoly).doit()) for k, v in al.items()}
    if 0 in al:
        print('[%s] WARNING mass term alpha_0 =' % name, sp.factor(al[0]))
    e1 = sp.simplify(flux_at_horizon(al, 1/f)/eta0)
    res[name] = e1
    print('[theta %s] e1 = eta1/eta0 =' % name, sp.collect(sp.expand(e1), Tj), '  (%.1fs)' % (time.time() - T0))
    sys.stdout.flush()

# Wald entropy ratio: perturb Riemann along the t-r binormal, Riemann independent; background (eps = 0) only
ddv = sp.Symbol('dd')
Xt = {(0, 1, 0, 1): 1, (1, 0, 1, 0): 1, (0, 1, 1, 0): -1, (1, 0, 0, 1): -1}
gb = [x.subs(eps, 0) for x in g]
gib = [x.subs(eps, 0) for x in gi]
Rlb = {k: (v.subs(eps, 0) + ddv*Xt.get(k, 0)) for k, v in Rl.items()}
Ricb = [[sum(gib[a]*Rlb[a, b, a, dd] for a in range(4)) for dd in range(4)] for b in range(4)]
Rsb = sum(gib[b]*Ricb[b][b] for b in range(4))


def Wlb(a, b, c, dd):
    gm = lambda i, j: gb[i] if i == j else 0
    return (Rlb[a, b, c, dd] - sp.Rational(1, 2)*(gm(a, c)*Ricb[b][dd] - gm(a, dd)*Ricb[b][c] - gm(b, c)*Ricb[a][dd] + gm(b, dd)*Ricb[a][c])
            + Rsb/6*(gm(a, c)*gm(b, dd) - gm(a, dd)*gm(b, c)))
Mb = sp.Matrix(6, 6, lambda i, j: Wlb(*PAIRS[i], *PAIRS[j])*gib[PAIRS[j][0]]*gib[PAIRS[j][1]])
Ib = 4*(Mb*Mb).trace(); Jb = 8*(Mb*Mb*Mb).trace()
dR = sp.diff(Rsb, ddv).subs(ddv, 0)
for name, inv in (('I', Ib), ('J', Jb)):
    ratio = sp.simplify((sp.diff(inv, ddv).subs(ddv, 0)/dR).subs(r, rh))
    print('[theta %s] sigma1 = s1/s0 = t0 *' % name, sp.factor(ratio))
    res['s' + name] = ratio*Tj[0]
print('total %.1fs' % (time.time() - T0))
xs = sp.Symbol('x')
for name in ('I', 'J'):
    dl = sp.simplify(res[name] - res['s' + name])
    print('[theta %s] delta = e1 - sigma1 =' % name, sp.collect(sp.expand(dl), Tj))
