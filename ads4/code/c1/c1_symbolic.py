r"""
c1_symbolic (C1_COMPUTE.md §3.1, §3.3, §3.4): exact symbolic checks used in the c1 derivation.

(A) d = 3 Hofman-Maldacena one-point matrix in the helicity basis: the t4 structure |eps_ij n_i n_j|^2/(eps*eps) - 2/(d^2-1) is purely
    helicity-FLIP (diagonal identically zero) and the t2 structure vanishes identically in d = 3.  Hence H = 1 + (t4/4)(e^{-4i phi} s+ + h.c.).
(B) AdS4 shock geometry ds^2 = [-du dv + dx^2 + dz^2 + g(u) F(x,z) du^2]/z^2 (g arbitrary, so delta(u) is covered):
    (B1) the only non-zero Weyl components W_{u i u j} (i, j in {x, z}) are g(u) x (traceless H2-Hessian of f), f = F/z, up to a factor;
    (B2) R_uu  = g(u) x (Box_H2 - 2) f up to a factor (so the profile propagator is (2 - Box)^{-1}, harmonic transform 1/(nu^2 + 9/4)).
(C) On H2 (unit radius): for the traceless Hessian D_ij f = nabla_i nabla_j f - (1/2) g_ij Box f:
    (C1) nabla^i nabla^j D_ij f = (1/2) Box (Box - 2) f   (identity used for the preserving vertex);
    (C2) D_ij f = 0 for f = 1 and for f = X^A (embedding coordinates, Box X^A = 2 X^A);
    (C3) the zonal harmonic P_{-1/2+i nu}(cosh L) has Box-eigenvalue -(nu^2 + 1/4).
(D) Algebra of the preserving O(t4) leading-Regge amplitude: P(Box) with P(0) = P(2) = 0 and deg P <= 2 gives
    alpha_pres(nu) = P(Box)/(2 - Box) |_{Box = -(nu^2+1/4)} = a (nu^2 + 1/4): no pole, and constant/nu^2 ratio locked to 1/4.

Usage: python3 c1_symbolic.py  -> cache/c1_symbolic.txt   (seconds)
"""
import os, sys, io
os.environ['OMP_NUM_THREADS'] = '1'
sys.dont_write_bytecode = True
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'cache', 'c1_symbolic.txt')
buf = io.StringIO()


def pr(s=""):
    print(s, flush=True); buf.write(str(s) + "\n")


# ---------------------------------------------------------------- (A) Hofman-Maldacena in d = 3
def part_A():
    pr("(A) d = 3 HM one-point structures in the helicity basis")
    phi = sp.symbols('phi', real=True)
    n = sp.Matrix([sp.cos(phi), sp.sin(phi)])
    ep = sp.Matrix([1, sp.I])/sp.sqrt(2); em = sp.Matrix([1, -sp.I])/sp.sqrt(2)
    E = {+1: ep*ep.T, -1: em*em.T}                      # spin-2 polarisation tensors, traceless symmetric, eps*eps = 1
    d = 3
    for h in (+1, -1):
        e = E[h]
        assert sp.simplify(e.trace()) == 0
        norm = sp.simplify(sum(sp.conjugate(e[i, j])*e[i, j] for i in range(2) for j in range(2)))
        pr(f"   helicity {h:+d}: Tr = 0, eps*eps = {norm}")
    # matrix elements <h'| structure |h> = eps_{h'}^* ... eps_h / norm
    def t2(hp, h):
        a, b = E[hp], E[h]
        return sp.simplify(sum(sp.conjugate(a[i, j])*b[i, l]*n[j]*n[l] for i in range(2) for j in range(2) for l in range(2))
                           - (1 if hp == h else 0)*sp.Rational(1, d - 1))
    def t4(hp, h):
        a, b = E[hp], E[h]
        an = sum(sp.conjugate(a[i, j])*n[i]*n[j] for i in range(2) for j in range(2))
        bn = sum(b[i, j]*n[i]*n[j] for i in range(2) for j in range(2))
        return sp.simplify(sp.expand(an*bn - (1 if hp == h else 0)*sp.Rational(2, d**2 - 1)))
    for hp in (+1, -1):
        for h in (+1, -1):
            pr(f"   <{hp:+d}| t2-structure |{h:+d}> = {t2(hp, h)};   <{hp:+d}| t4-structure |{h:+d}> = {sp.simplify(sp.expand_complex(t4(hp, h)).rewrite(sp.exp))}")
    pr("   => t2 absent in d = 3; t4 purely off-diagonal (flip), modulus 1/4: H = 1 + (t4/4)(e^{-+4 i phi} flips) [agrees with Q8 / LIGHT_INPUTS t2]")


# ---------------------------------------------------------------- (B) AdS4 shock Weyl tensor
def christoffel(g, X):
    n = len(X); gi = g.inv()
    return [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d])) for d in range(n))/2)
              for c in range(n)] for b in range(n)] for a in range(n)]


def riemann(g, X, G):
    n = len(X)
    R = [[[[0]*n for _ in range(n)] for _ in range(n)] for _ in range(n)]   # R^a_{bcd}
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    e = sp.diff(G[a][b][d], X[c]) - sp.diff(G[a][b][c], X[d])
                    e += sum(G[a][c][k]*G[k][b][d] - G[a][d][k]*G[k][b][c] for k in range(n))
                    R[a][b][c][d] = e
    Rl = [[[[sp.expand(sum(g[a, k]*R[k][b][c][d] for k in range(n))) for d in range(n)] for c in range(n)] for b in range(n)] for a in range(n)]
    return Rl


def part_B():
    pr("\n(B) AdS4 shock: ds^2 = [-du dv + dx^2 + dz^2 + g(u) F(x,z) du^2]/z^2")
    u, v, x, z = sp.symbols('u v x z', real=True)
    gf = sp.Function('g')(u); F = sp.Function('F')(x, z)
    X = [u, v, x, z]
    g = sp.zeros(4)
    g[0, 1] = g[1, 0] = -sp.Rational(1, 2)/z**2
    g[2, 2] = g[3, 3] = 1/z**2
    g[0, 0] = gf*F/z**2
    G = christoffel(g, X)
    Rl = riemann(g, X, G)
    gi = g.inv()
    Ric = sp.Matrix(4, 4, lambda b, d: sp.simplify(sum(gi[a, c]*Rl[a][b][c][d] for a in range(4) for c in range(4))))
    Rs = sp.simplify(sum(gi[b, d]*Ric[b, d] for b in range(4) for d in range(4)))
    pr(f"   Ricci scalar = {Rs}   (AdS4 unit radius: -12)")
    # Weyl in 4D
    def W(a, b, c, d):
        val = Rl[a][b][c][d]
        val -= sp.Rational(1, 2)*(g[a, c]*Ric[b, d] - g[a, d]*Ric[b, c] - g[b, c]*Ric[a, d] + g[b, d]*Ric[a, c])
        val += Rs/6*(g[a, c]*g[b, d] - g[a, d]*g[b, c])
        return sp.simplify(val)
    # H2 with metric (dx^2 + dz^2)/z^2 and f = F/z
    f = F/z
    h2 = sp.diag(1/z**2, 1/z**2); Y = [x, z]; h2i = h2.inv()
    Gh = [[[sp.simplify(sum(h2i[a, d]*(sp.diff(h2[d, b], Y[c]) + sp.diff(h2[d, c], Y[b]) - sp.diff(h2[b, c], Y[d])) for d in range(2))/2)
            for c in range(2)] for b in range(2)] for a in range(2)]
    Hess = sp.Matrix(2, 2, lambda i, j: sp.diff(f, Y[i], Y[j]) - sum(Gh[k][i][j]*sp.diff(f, Y[k]) for k in range(2)))
    Box = sp.simplify(sum(h2i[i, j]*Hess[i, j] for i in range(2) for j in range(2)))
    D = sp.simplify(Hess - h2*Box/2)
    idx = {'x': 2, 'z': 3}
    ratios = []
    for (i, j) in (('x', 'x'), ('x', 'z'), ('z', 'z')):
        w = W(0, idx[i], 0, idx[j])
        r = sp.simplify(w/D[Y.index(sp.Symbol(i, real=True)), Y.index(sp.Symbol(j, real=True))]) if D[Y.index(sp.Symbol(i, real=True)), Y.index(sp.Symbol(j, real=True))] != 0 else None
        ratios.append(r)
        pr(f"   W_u{i}u{j} / D_{i}{j}f = {r}")
    # other Weyl components (sample) vanish
    others = []
    for (a, b, c, d) in [(0, 2, 2, 3), (2, 3, 2, 3), (0, 1, 0, 1), (0, 1, 2, 3), (1, 2, 1, 2), (0, 2, 1, 2)]:
        others.append(W(a, b, c, d))
    pr(f"   sample of other Weyl components (u v x z = 0 1 2 3): {others}")
    Ruu = sp.simplify(Ric[0, 0] - (-3)*g[0, 0])          # Einstein eq. with Lambda: R_mn = -3 g_mn in vacuum
    pr(f"   R_uu + 3 g_uu = {sp.factor(Ruu)}")
    pr(f"   (Box_H2 - 2) f * g(u)/z = {sp.factor(sp.simplify((Box - 2*f)*gf/z))}")
    pr(f"   ratio = {sp.simplify(Ruu/((Box - 2*f)*gf/z))}")


# ---------------------------------------------------------------- (C) H2 identities
def part_C():
    pr("\n(C) H2 identities (unit radius, coordinates (x, z), metric (dx^2+dz^2)/z^2)")
    x, z = sp.symbols('x z', positive=True)
    f = sp.Function('f')(x, z)
    Y = [x, z]
    h = sp.diag(1/z**2, 1/z**2); hi = h.inv()
    Gm = [[[sp.simplify(sum(hi[a, d]*(sp.diff(h[d, b], Y[c]) + sp.diff(h[d, c], Y[b]) - sp.diff(h[b, c], Y[d])) for d in range(2))/2)
            for c in range(2)] for b in range(2)] for a in range(2)]
    def hess(s):
        return sp.Matrix(2, 2, lambda i, j: sp.diff(s, Y[i], Y[j]) - sum(Gm[k][i][j]*sp.diff(s, Y[k]) for k in range(2)))
    def box(s):
        H = hess(s); return sp.simplify(sum(hi[i, j]*H[i, j] for i in range(2) for j in range(2)))
    def divdiv(T):    # nabla^i nabla^j T_ij for a symmetric covariant 2-tensor
        # nabla_k T_ij
        def cov(i, j, k):
            return sp.diff(T[i, j], Y[k]) - sum(Gm[l][k][i]*T[l, j] + Gm[l][k][j]*T[i, l] for l in range(2))
        V = [sp.simplify(sum(hi[k, i]*cov(i, j, k) for i in range(2) for k in range(2))) for j in range(2)]   # V_j = nabla^i T_ij
        return sp.simplify(sum(hi[k, j]*(sp.diff(V[j], Y[k]) - sum(Gm[l][k][j]*V[l] for l in range(2))) for j in range(2) for k in range(2)))
    D = hess(f) - h*box(f)/2
    lhs = divdiv(D)
    rhs = sp.Rational(1, 2)*box(box(f) - 2*f)
    pr(f"   (C1) nabla^i nabla^j D_ij f - (1/2) Box(Box - 2) f = {sp.simplify(lhs - rhs)}")
    # embedding coordinates of H2 in Poincare coords: X^0 = (1 + x^2 + z^2)/(2z), X^1 = x/z, X^2 = (1 - x^2 - z^2)/(2z)
    Xs = [(1 + x**2 + z**2)/(2*z), x/z, (1 - x**2 - z**2)/(2*z)]
    for k, Xa in enumerate(Xs):
        Dx = sp.simplify(hess(Xa) - h*box(Xa)/2)
        pr(f"   (C2) X^{k}: Box X - 2X = {sp.simplify(box(Xa) - 2*Xa)},  D_ij X = {list(Dx)}")
    pr(f"   (C2) f = 1: D_ij 1 = {list(sp.simplify(hess(sp.Integer(1)) - h*box(sp.Integer(1))/2))}")
    # (C3) zonal: in geodesic polar coordinates Box = d^2/dL^2 + coth L d/dL on functions of L
    L, nu = sp.symbols('L nu', positive=True)
    import mpmath as mp
    for nv in (0.3, 2.0):
        for Lv in (0.4, 1.3):
            P = lambda LL: mp.legenp(-0.5 + 1j*nv, 0, mp.cosh(LL), type=3)
            bx = mp.diff(P, Lv, 2) + mp.coth(Lv)*mp.diff(P, Lv, 1)
            pr(f"   (C3) nu={nv}, L={Lv}: Box P / P = {mp.nstr(bx/P(Lv), 12)}   vs -(nu^2+1/4) = {-(nv**2 + 0.25)}")


def part_D():
    pr("\n(D) preserving O(t4) leading-Regge amplitude: alpha_pres(nu) = P(Box)/(2 - Box), Box -> -(nu^2+1/4)")
    B, a, b, c, nu = sp.symbols('B a b c nu')
    P = a*B**2 + b*B + c
    sol = sp.solve([P.subs(B, 0), P.subs(B, 2)], [b, c], dict=True)[0]
    Pf = sp.factor(P.subs(sol))
    alpha = sp.simplify((Pf/(2 - B)).subs(B, -(nu**2 + sp.Rational(1, 4))))
    pr(f"   P with P(0) = P(2) = 0: {Pf};  alpha_pres = {sp.expand(alpha)}")
    Pg = sp.expand((B - 2)*(a*B + b))
    alpha_g = sp.expand(sp.simplify((Pg/(2 - B)).subs(B, -(nu**2 + sp.Rational(1, 4)))))
    pr(f"   for comparison, P(2) = 0 only: P = (B-2)(aB+b) -> alpha = {alpha_g}  (a constant -b survives if P(0) != 0)")


if __name__ == '__main__':
    part_A(); part_B(); part_C(); part_D()
    with open(OUT, 'w') as fo:
        fo.write(buf.getvalue())
