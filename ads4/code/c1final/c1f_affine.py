#!/usr/bin/env python3
"""C1_FINAL step A1 (exact algebra, sympy).

(a) The k = 2 CLS-D subtraction factor f = v' - u' (J_e = 2), with
        u' = x12^2 x34^2/(x13^2 x24^2),  v' = x14^2 x23^2/(x13^2 x24^2),
    evaluated with the two light-ray points x1, x3 on future null infinity (light rays based at x0 = infinity),
        x1 = (u1 + r1, r1 n1),  x3 = (u3 + r3, r3 n3),  r1, r3 -> infinity,
    is AFFINE in the retarded times (u1, u3) and has NO u1*u3 term:
        f = 2/((1 - n1.n3) x24^2) * [ u1 nb3.x42 + u3 nb1.x24 + (nb1.x4)(nb3.x2) - (nb1.x2)(nb3.x4) ],
    with nb.x := x^0 - n.xvec  (so that x_{1j}^2 -> -2 r1 (u1 - nb1.x_j)).
(b) The (2,4) part of the EE dual structure, G = H24 - 2 V_{2,04} V_{4,02} (CLS-D building blocks), tends to
    x24^2 (z2.z4) (times a z-independent constant) as x0 -> infinity: the (2,4) polarisations are contracted by the identity.
Metric (-,+,+). Output: cache/c1f_affine.txt
"""
import sympy as sp, os, sys

out = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

eta = sp.diag(-1, 1, 1)
def dot(a, b): return (sp.Matrix(a).T * eta * sp.Matrix(b))[0]

# ---------------- (a) f at x0 = infinity ----------------
u1, u3, e1, e3 = sp.symbols('u1 u3 epsilon1 epsilon3', real=True)
c1, s1, c3, s3 = sp.symbols('c1 s1 c3 s3', real=True)          # n_i = (c_i, s_i), c^2+s^2 = 1 imposed later
x2 = sp.symbols('t2 a2 b2', real=True); x4 = sp.symbols('t4 a4 b4', real=True)
r1, r3 = 1/e1, 1/e3
X1 = [u1 + r1, r1*c1, r1*s1]; X3 = [u3 + r3, r3*c3, r3*s3]
def sq(a, b): v = [a[i] - b[i] for i in range(3)]; return dot(v, v)
sub = {s1**2: 1 - c1**2, s3**2: 1 - c3**2}
def red(e): return sp.expand(sp.expand(e).subs({s1**2: 1 - c1**2}).subs({s3**2: 1 - c3**2}))
x12, x14 = red(sq(X1, x2)), red(sq(X1, x4)); x32, x34 = red(sq(X3, x2)), red(sq(X3, x4))
x13, x24 = red(sq(X1, X3)), red(sq(x2, x4))
# leading behaviour: x1j^2 = (1/e1)*A1j + O(1);  x13^2 = (1/(e1 e3))*A13 + ...
A12 = sp.limit(x12*e1, e1, 0); A14 = sp.limit(x14*e1, e1, 0)
A32 = sp.limit(x32*e3, e3, 0); A34 = sp.limit(x34*e3, e3, 0)
A13 = sp.simplify(sp.limit(sp.limit(x13*e1*e3, e1, 0), e3, 0))
P("x12^2 * eps1 ->", sp.factor(A12)); P("x13^2 * eps1*eps3 ->", sp.factor(A13))
# exact limit of the cross ratios (the eps powers cancel identically)
up = sp.simplify(A12*A34/(A13*x24)); vp = sp.simplify(A14*A32/(A13*x24))
f = sp.simplify(vp - up)
P("f = v'-u' at x0=inf:", f)
P("d^2 f/du1du3 =", sp.simplify(sp.diff(f, u1, u3)), "; d^2f/du1^2 =", sp.simplify(sp.diff(f, u1, 2)),
  "; d^2f/du3^2 =", sp.simplify(sp.diff(f, u3, 2)))
nb = lambda c, s, x: x[0] - c*x[1] - s*x[2]
x42 = [x4[i] - x2[i] for i in range(3)]; x24v = [x2[i] - x4[i] for i in range(3)]
claim = 2/((1 - c1*c3 - s1*s3)*x24) * (u1*nb(c3, s3, x42) + u3*nb(c1, s1, x24v)
                                      + nb(c1, s1, x4)*nb(c3, s3, x2) - nb(c1, s1, x2)*nb(c3, s3, x4))
diffa = sp.simplify(red(sp.together(f - claim).as_numer_denom()[0]))
P("f - closed form (numerator, after n_i^2 = 1):", diffa)
# subleading check: the finite-r correction vanishes as r -> infinity (numerical spot check, exact rationals)
import random
random.seed(3)
def spot(R):
    vals = {u1: sp.Rational(3, 7), u3: sp.Rational(-2, 5), c1: sp.Rational(3, 5), s1: sp.Rational(4, 5),
            c3: sp.Rational(-5, 13), s3: sp.Rational(12, 13), x2[0]: sp.Rational(1, 3), x2[1]: sp.Rational(1, 7),
            x2[2]: sp.Rational(-1, 4), x4[0]: sp.Rational(-2, 3), x4[1]: sp.Rational(1, 5), x4[2]: sp.Rational(2, 9),
            e1: sp.Rational(1, R), e3: sp.Rational(1, 2*R)}
    ff = (x14*x32 - x12*x34)/(x13*x24)
    return sp.N(ff.subs(vals), 30), sp.N(claim.subs(vals), 30)
for R in (10**3, 10**6, 10**9):
    a, b = spot(R); P(f"  finite r = {R}: exact f = {a},  limit = {b},  diff = {sp.N(a-b, 5)}")
okA = (sp.simplify(sp.diff(f, u1, u3)) == 0 and sp.simplify(sp.diff(f, u1, 2)) == 0
       and sp.simplify(sp.diff(f, u3, 2)) == 0 and diffa == 0)
P("CHECK (a) f affine in (u1,u3), no u1*u3 term, closed form exact:", okA)

# ---------------- (b) EE (2,4) structure at x0 -> infinity ----------------
L = sp.symbols('L', positive=True)
z2 = sp.symbols('p0 p1 p2'); z4 = sp.symbols('q0 q1 q2'); z0 = sp.symbols('w0 w1 w2')
m = [sp.Rational(1, 3), sp.Rational(2, 3), sp.Rational(-2, 3)]   # direction of x0 = L*m (spacelike: -1/9+4/9+4/9 > 0)
def emb(x):   # X = (X+, X-, x),  X.Y = -(X+Y- + X-Y+)/2 + x.y
    return (1, dot(x, x), list(x))
def embZ(x, z): return (0, 2*dot(x, z), list(z))
def edot(A, B): return -(A[0]*B[1] + A[1]*B[0])/2 + dot(A[2], B[2])
X2, X4 = emb(list(x2)), emb(list(x4)); Z2, Z4 = embZ(list(x2), list(z2)), embZ(list(x4), list(z4))
x0 = [L*mi for mi in m]; X0 = emb(x0); Z0 = embZ(x0, list(z0))
def V(Zi, Xi, Xj, Xk): return (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)
def H(Zi, Xi, Zj, Xj): return -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))
G = H(Z2, X2, Z4, X4) - 2*V(Z2, X2, X0, X4)*V(Z4, X4, X0, X2)
Ginf = sp.limit(sp.simplify(G), L, sp.oo)
target = dot([x2[i]-x4[i] for i in range(3)], [x2[i]-x4[i] for i in range(3)])*dot(list(z2), list(z4))
ratio = sp.simplify(Ginf/target)
P("G = H24 - 2 V_{2,04} V_{4,02} at x0 -> inf:", sp.factor(Ginf))
P("G_inf / (x24^2 z2.z4) =", ratio)
okB = ratio.free_symbols == set()
P("CHECK (b) EE (2,4) contraction -> const * x24^2 (z2.z4) [identity contraction]:", okB)
P("ALL CHECKS PASS:", okA and okB)
os.makedirs(os.path.join(os.path.dirname(__file__), 'cache'), exist_ok=True)
open(os.path.join(os.path.dirname(__file__), 'cache', 'c1f_affine.txt'), 'w').write("\n".join(out) + "\n")
