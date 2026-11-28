"""Independent exact check of SM Prop. S18.2 (Eq. T4 / Letter Eq. KKric).

Claim: for g_D = g_ab dx^a dx^b + G_ij (dy^i + A^i_a dx^a)(dy^j + A^j_b dx^b) on a 3-dim Lorentzian base,
with u = (1/2) log det G, ghat = e^{2u} g, P_a = G^{-1} d_a G, F^i = dA^i, and k null for g,
K = k^a (d_a - A^i_a d_i):

  Ric_ghat(k,k) = Ric_D(K,K) + (k u)^2 + 1/4 Tr(P_k^2) + 1/2 G_ij (i_k F^i)_a (i_k F^j)^a.

We verify the stronger tensor statement: the symmetric base tensor
  Delta_ab = Ric_ghat_ab - [Ric_D(E_a,E_b) + du_a du_b + 1/4 Tr(P_a P_b) + 1/2 G_ij F^i_ac F^j_bd g^cd]
is pure trace (proportional to g_ab) at a point, so that it vanishes on every null k.
Ricci is computed from exact 2-jets at the origin (random rational polynomial data).
"""
import random
import sympy as sp

random.seed(12345)


def rnd():
    return sp.Rational(random.randint(-9, 9), random.randint(5, 13))


def ricci_at0(gfun, X):
    n = len(X)
    zero = {x: 0 for x in X}
    g0 = sp.Matrix(n, n, lambda i, j: gfun[i, j].subs(zero))
    gi = g0.inv()
    dg = [[[sp.diff(gfun[i, j], X[k]).subs(zero) for k in range(n)] for j in range(n)] for i in range(n)]
    ddg = [[[[sp.diff(gfun[i, j], X[k], X[l]).subs(zero) for l in range(n)] for k in range(n)]
             for j in range(n)] for i in range(n)]
    # Christoffel Gamma^a_bc at 0 and its derivative d_d Gamma^a_bc at 0
    Gam = [[[sum(gi[a, e] * (dg[e][b][c] + dg[e][c][b] - dg[b][c][e]) for e in range(n)) / 2
             for c in range(n)] for b in range(n)] for a in range(n)]
    dgi = [[[-sum(gi[a, p] * dg[p][q][d] * gi[q, b] for p in range(n) for q in range(n))
             for d in range(n)] for b in range(n)] for a in range(n)]
    dGam = [[[[sum(dgi[a][e][d] * (dg[e][b][c] + dg[e][c][b] - dg[b][c][e]) / 2
                   + gi[a, e] * (ddg[e][b][c][d] + ddg[e][c][b][d] - ddg[b][c][e][d]) / 2 for e in range(n))
               for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    # R_bd = d_a Gam^a_bd - d_d Gam^a_ba + Gam^a_ae Gam^e_bd - Gam^a_de Gam^e_ba
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            s = 0
            for a in range(n):
                s += dGam[a][b][d][a] - dGam[a][b][a][d]
                for e in range(n):
                    s += Gam[a][a][e] * Gam[e][b][d] - Gam[a][d][e] * Gam[e][b][a]
            Ric[b, d] = sp.nsimplify(s)
    return Ric, g0


def run(m, cP=sp.Rational(1,4), cF=sp.Rational(1,2)):
    t, x, z = sp.symbols('t x z')
    base = [t, x, z]
    ys = sp.symbols('y0:%d' % m)
    def poly2():
        mon = [1, t, x, z, t * t, x * x, z * z, t * x, t * z, x * z]
        return sum(rnd() * mm for mm in mon[1:])
    eta = sp.diag(-1, 1, 1)
    g = sp.Matrix(3, 3, lambda i, j: 0)
    for i in range(3):
        for j in range(i, 3):
            v = eta[i, j] + poly2() / 3
            g[i, j] = v
            g[j, i] = v
    G = sp.Matrix(m, m, lambda i, j: 0)
    for i in range(m):
        for j in range(i, m):
            v = (1 if i == j else 0) + poly2() / 4
            G[i, j] = v
            G[j, i] = v
    A = sp.Matrix(m, 3, lambda i, a: poly2())
    n = 3 + m
    X = base + list(ys)
    gD = sp.zeros(n, n)
    for a in range(3):
        for b in range(3):
            gD[a, b] = g[a, b] + sum(G[i, j] * A[i, a] * A[j, b] for i in range(m) for j in range(m))
    for a in range(3):
        for i in range(m):
            v = sum(G[i, j] * A[j, a] for j in range(m))
            gD[a, 3 + i] = v
            gD[3 + i, a] = v
    for i in range(m):
        for j in range(m):
            gD[3 + i, 3 + j] = G[i, j]
    RicD, gD0 = ricci_at0(gD, X)
    zero = {s: 0 for s in X}
    A0 = A.subs(zero)
    # horizontal frame E_a = d_a - A^i_a d_i
    E = sp.zeros(n, 3)
    for a in range(3):
        E[a, a] = 1
        for i in range(m):
            E[3 + i, a] = -A0[i, a]
    RicDh = E.T * RicD * E
    detG = G.det()
    ghat = sp.Matrix(3, 3, lambda i, j: sp.expand(detG * g[i, j]))
    Rich, g0 = ricci_at0(ghat, base)
    detG0 = detG.subs(zero)
    du = [sp.diff(detG, s).subs(zero) / (2 * detG0) for s in base]
    G0 = G.subs(zero)
    G0i = G0.inv()
    P = [G0i * G.diff(s).subs(zero) for s in base]
    F = [[[sp.diff(A[i, b], base[a]).subs(zero) - sp.diff(A[i, a], base[b]).subs(zero) for b in range(3)]
          for a in range(3)] for i in range(m)]
    gi = g0.inv()
    extra = sp.zeros(3, 3)
    for a in range(3):
        for b in range(3):
            fterm = sum(G0[i, j] * F[i][a][c] * gi[c, d] * F[j][b][d]
                        for i in range(m) for j in range(m) for c in range(3) for d in range(3))
            extra[a, b] = du[a] * du[b] + cP*(P[a] * P[b]).trace() + cF*fterm
    Delta = (Rich - RicDh - extra).applyfunc(sp.nsimplify)
    c = sp.nsimplify(Delta[0, 0] / g0[0, 0])
    resid = (Delta - c * g0).applyfunc(sp.simplify)
    return resid, c, Delta


if __name__ == '__main__':
    for m in (1, 2):
        resid, c, Delta = run(m)
        print('m =', m, ' max |Delta - c g| =', max(abs(v) for v in resid), '  trace coefficient c =', c)
        r2,_,_ = run(m, cF=1)
        r3,_,_ = run(m, cP=sp.Rational(1,2))
        print('   controls (F coeff 1; P coeff 1/2): max resid =', max(abs(v) for v in r2), max(abs(v) for v in r3))
    print('done')
