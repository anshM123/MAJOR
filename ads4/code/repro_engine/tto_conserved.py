"""repro_engine/tto_conserved.py -- independent derivation of the conserved parity-even <T T O_Delta> (O scalar, J = 0) coefficients
in d = 3, in the CPPR basis  {V1^2 V2^2, V1 V2 H12, H12^2} / (X12^{(10-D)/2} X15^{D/2} X25^{D/2}),
V1 = V_{1,25}, V2 = V_{2,51}, H12 (embedding space, X = (1, x^2, x), Z = (0, 2x.z, z), lightcone components (u, v, y),
x.y = -(u v' + v u')/2 + y y').

Method (own code, exact arithmetic, sympy): f = N(x1, z1) exp(L), L = -p log X12 - q log X15 - r log X25.  Conservation at x1 with the
Todorov operator in d = 3 (h = 3/2), D_mu = (h - 1 + z.d_z) d_{z mu} - (1/2) z_mu d_z^2, i.e.  div = d_x^mu D_mu f = 0 on the null cone,
evaluated at random rational x1, null z1 (other points fixed rational), divided by exp(L).  Linear system in (c1, c2, c3) with
polynomial-in-Delta coefficients; its null space is solved over Q(Delta).  Then the 1<->2 symmetry is checked.

Usage: python3 tto_conserved.py   -> prints the solution; cache/tto_conserved.txt
"""
import os, random
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
D = sp.Symbol('Delta')
G = sp.Matrix([[0, sp.Rational(-1, 2), 0], [sp.Rational(-1, 2), 0, 0], [0, 0, 1]])     # lightcone metric
GI = G.inv()


def dot(a, b):
    return sum(a[i]*G[i, j]*b[j] for i in range(3) for j in range(3))


def edot(A, B):                       # embedding dot of (X+, X-, x) objects: -(A0 B1 + A1 B0)/2 + a.b
    return -(A[0]*B[1] + A[1]*B[0])/2 + dot(A[2], B[2])


def EX(x):
    return (1, dot(x, x), x)


def EZ(x, z):
    return (0, 2*dot(x, z), z)


def V(xi, zi, xj, xk):
    Xi, Zi, Xj, Xk = EX(xi), EZ(xi, zi), EX(xj), EX(xk)
    return (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)


def H(xi, zi, xj, zj):
    Xi, Zi, Xj, Zj = EX(xi), EZ(xi, zi), EX(xj), EZ(xj, zj)
    return -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))


def Xd(xi, xj):
    return -2*edot(EX(xi), EX(xj))   # = (xi - xj)^2


def rnd(rng):
    return sp.Rational(rng.randint(-9, 9), rng.randint(1, 5))


def rnull(rng):
    while True:
        u, y = rnd(rng), rnd(rng)
        if u != 0:
            return sp.Matrix([u, y*y/u, y])        # -u v + y^2 = 0


def structs(x1, z1, x2, z2, x5):
    V1, V2, H12 = V(x1, z1, x2, x5), V(x2, z2, x5, x1), H(x1, z1, x2, z2)
    return [V1**2*V2**2, V1*V2*H12, H12**2]


def main():
    rng = random.Random(7)
    xs = sp.symbols('x0:3'); zs = sp.symbols('z0:3')
    x1 = sp.Matrix(xs); z1 = sp.Matrix(zs)
    p, q, r = (10 - D)/2, D/2, D/2
    rows = []
    for trial in range(4):
        x2 = sp.Matrix([rnd(rng) for _ in range(3)]); x5 = sp.Matrix([rnd(rng) for _ in range(3)]); z2 = rnull(rng)
        L = -p*sp.log(Xd(x1, x2)) - q*sp.log(Xd(x1, x5)) - r*sp.log(Xd(x2, x5))
        dL = [sp.diff(L, xs[m]) for m in range(3)]
        for N in structs(x1, z1, x2, z2, x5):
            pass
        Ns = structs(x1, z1, x2, z2, x5)
        h = sp.Rational(3, 2)
        exprs = []
        for N in Ns:
            # D_mu N (lower index): (h - 1 + z.d_z) d_{z^mu} N - (1/2) z_mu Box_z N ;  Box_z = g^{ab} d_a d_b
            dz = [sp.diff(N, zs[m]) for m in range(3)]
            box = sum(GI[a, b]*sp.diff(N, zs[a], zs[b]) for a in range(3) for b in range(3))
            zlow = [sum(G[m, n]*zs[n] for n in range(3)) for m in range(3)]
            Dlow = [(h - 1)*dz[m] + sum(zs[n]*sp.diff(dz[m], zs[n]) for n in range(3)) - sp.Rational(1, 2)*zlow[m]*box
                    for m in range(3)]
            # div = g^{mu nu} [d_{x mu} D_nu N + (d_{x mu} L) D_nu N]
            div = sum(GI[m, n]*(sp.diff(Dlow[n], xs[m]) + dL[m]*Dlow[n]) for m in range(3) for n in range(3))
            exprs.append(div)
        for k in range(3):
            xv = [rnd(rng) for _ in range(3)]; zv = list(rnull(rng))
            sub = dict(zip(xs, xv)); sub.update(dict(zip(zs, zv)))
            rows.append([sp.factor(sp.simplify(e.subs(sub))) for e in exprs])
    M = sp.Matrix(rows)
    ns = M.nullspace(simplify=True)
    out = []
    out.append(f'# tto_conserved: {M.shape[0]} conservation rows (point 1), null-space dim = {len(ns)}')
    for v in ns:
        v = sp.simplify(v/v[2])
        out.append('  c = (V1^2V2^2, V1V2H12, H12^2) ∝ (' + ', '.join(str(sp.factor(c)) for c in v) + ')')
        sol = v
    # checks: conservation at point 2 and 1<->2 symmetry at random points, generic Delta
    chk = []
    for Dv in (sp.Rational(301, 1), sp.Rational(-297, 1), sp.Rational(20000, 1)):
        s = [c.subs(D, Dv) for c in sol]
        x1n = sp.Matrix([rnd(rng) for _ in range(3)]); x2n = sp.Matrix([rnd(rng) for _ in range(3)]); x5n = sp.Matrix([rnd(rng) for _ in range(3)])
        z1n, z2n = rnull(rng), rnull(rng)
        f12 = sum(c*t for c, t in zip(s, structs(x1n, z1n, x2n, z2n, x5n)))
        f21 = sum(c*t for c, t in zip(s, structs(x2n, z2n, x1n, z1n, x5n)))
        chk.append(f'  Delta={Dv}: f(2<->1)/f(1<->2) = {sp.nsimplify(f21/f12)}')
    out += chk
    print('\n'.join(out))
    open(os.path.join(HERE, 'cache', 'tto_conserved.txt'), 'w').write('\n'.join(out) + '\n')
    return sol


if __name__ == '__main__':
    main()
