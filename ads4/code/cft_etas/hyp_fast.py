"""
Fast hyperbolic shear engine (symengine): O(eps^2) Lagrangian densities of R, I = C^2, J = C^3 for the TT shear mode on the
AdS4 hyperbolic black hole (conventions of hyp_engine.py).  phi(r) is handled through jets p0..p4 with the total derivative
Dr = d_r + p1 d_p0 + p2 d_p1 + ... ; the O(eps^2) parts are returned as quadratic forms in the jets and cached (pickle of
sympy srepr strings) for hyp_response.py.
"""
import itertools, time, pickle, sys, math
import symengine as se
T0 = time.time()
r, u, rh, eps = se.symbols('r u r_h epsilon')
p = se.symbols('p0:6')
w = 1 - u**2
mu = rh**3 - rh
f = r**2 - 1 - mu/r


def trunc(e, n=2):
    e = se.expand(e)
    out = 0
    for k in range(n + 1):
        out += e.diff(eps, k).subs({eps: 0})/math.factorial(k)*eps**k
    return se.expand(out)


def Dr(e):
    return se.diff(e, r) + sum(p[k + 1]*se.diff(e, p[k]) for k in range(5))


def D(e, c):
    if c == 1: return Dr(e)
    if c == 2: return se.diff(e, u)
    return 0

A = p[0]*w
E = 1 + eps*A + eps**2*A**2/2
Ei = 1 - eps*A + eps**2*A**2/2
g = [-f, 1/f, r**2*E/w**2, r**2*Ei/w]
gi = [-1/f, f, trunc(w**2/r**2*Ei), trunc(w*E/r**2)]
Gam = {}
for a, b, c in itertools.product(range(4), repeat=3):
    v = 0
    if a == b: v += D(g[a], c)
    if a == c: v += D(g[a], b)
    if b == c: v -= D(g[b], a)
    Gam[a, b, c] = trunc(gi[a]*v/2) if v != 0 else 0
print('Gamma %.1fs' % (time.time() - T0)); sys.stdout.flush()
Rl = {}
for a, b, c, d in itertools.product(range(4), repeat=4):
    if c >= d: continue
    v = D(Gam[a, b, d], c) - D(Gam[a, b, c], d) + sum(Gam[a, c, h]*Gam[h, b, d] - Gam[a, d, h]*Gam[h, b, c] for h in range(4))
    v = trunc(g[a]*v)
    Rl[a, b, c, d] = v; Rl[a, b, d, c] = -v
for a, b, c in itertools.product(range(4), repeat=3):
    Rl[a, b, c, c] = 0
print('Riemann %.1fs' % (time.time() - T0)); sys.stdout.flush()
Ric = [[trunc(sum(gi[a]*Rl[a, b, a, d] for a in range(4))) for d in range(4)] for b in range(4)]
Rs = trunc(sum(gi[b]*Ric[b][b] for b in range(4)))
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
gm = lambda i, j: g[i] if i == j else 0


def Wl(a, b, c, d):
    return (Rl[a, b, c, d] - (gm(a, c)*Ric[b][d] - gm(a, d)*Ric[b][c] - gm(b, c)*Ric[a][d] + gm(b, d)*Ric[a][c])/2
            + Rs/6*(gm(a, c)*gm(b, d) - gm(a, d)*gm(b, c)))
M = [[trunc(Wl(*PAIRS[i], *PAIRS[j])*gi[PAIRS[j][0]]*gi[PAIRS[j][1]]) for j in range(6)] for i in range(6)]
print('Weyl %.1fs' % (time.time() - T0)); sys.stdout.flush()
M2 = [[trunc(sum(M[i][k]*M[k][j] for k in range(6))) for j in range(6)] for i in range(6)]
I = trunc(4*sum(M2[i][i] for i in range(6)))
J = trunc(8*sum(M2[i][k]*M[k][i] for i in range(6) for k in range(6)))
print('I, J %.1fs' % (time.time() - T0)); sys.stdout.flush()
out = {}
for name, e in (('R', Rs), ('I', I), ('J', J)):
    c0 = se.expand(e.subs({eps: 0}))
    c1 = se.expand(e.diff(eps).subs({eps: 0}))
    c2 = se.expand(e.diff(eps, 2).subs({eps: 0})/2)
    quad = {}
    for a in range(5):
        for b in range(a, 5):
            q = se.diff(c2, p[a], p[b])
            if a == b: q = q/2
            q = se.expand(q.subs({pk: 0 for pk in p}))
            if q != 0:
                quad[(a, b)] = str(q)
    out[name] = {'c0': str(c0), 'c1': str(c1), 'quad': quad}
    print(name, 'c0 terms:', len(str(c0)), ' c1 =', 'zero' if c1 == 0 else 'NONZERO(len %d)' % len(str(c1)), ' quad keys:', list(quad))
pickle.dump(out, open('cache_hyp_quad.pkl', 'wb'))
print('saved cache_hyp_quad.pkl  total %.1fs' % (time.time() - T0))
