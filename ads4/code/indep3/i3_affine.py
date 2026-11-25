"""INDEP_CHECK3, claim 1 step (i) and (iii)-geometry: independent re-derivation.

(a) f = v' - u' with x1, x3 on future null infinity (light rays based at x0 = spatial infinity):
    x_{1,3} = (u + r, r n), r -> oo.  Derive the limit symbolically and check it is affine in (u1, u3).
(b) EE (2,4) building block G = H24 - 2 V_{2,04} V_{4,02} at X0 -> I (point at infinity) equals x24^2 (z2.z4),
    with my own embedding conventions X = (X+, X-, x) = (1, x^2, x), X.Y = -(X+Y- + X-Y+)/2 + x.y,
    Z = (0, 2 x.z, z).  Checked symbolically and with random numbers at large finite x0.
Mostly-plus signature in d = 3.
"""
import sympy as sp
import random

out = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s); out.append(s)

# ---------- (a) ----------
r, u1, u3 = sp.symbols('r u1 u3', real=True)
c1, s1, c3, s3 = sp.symbols('c1 s1 c3 s3', real=True)    # n1 = (c1, s1), n3 = (c3, s3), unit
x2 = sp.symbols('t2 a2 b2', real=True); x4 = sp.symbols('t4 a4 b4', real=True)
def sq(a, b):  # (a-b)^2 mostly plus
    return -(a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2
X1 = (u1 + r, r*c1, r*s1); X3 = (u3 + r, r*c3, r*s3)
x12, x34, x13, x24 = sq(X1, x2), sq(X3, x4), sq(X1, X3), sq(x2, x4)
x14, x23 = sq(X1, x4), sq(x2, X3)
unit = {s1**2: 1 - c1**2, s3**2: 1 - c3**2}
up = x12*x34/(x13*x24); vp = x14*x23/(x13*x24)
f = sp.simplify(vp - up)
# leading large-r limit
sub = {s1: sp.sqrt(1 - c1**2), s3: sp.sqrt(1 - c3**2)}
numer = sp.Poly(sp.expand((x14*x23 - x12*x34).subs(sub)), r); den = sp.Poly(sp.expand(x13.subs(sub)), r)
log("(a) degree in r: numerator", numer.degree(), "x13", den.degree())
flim = sp.simplify(numer.LC()/(den.LC()*x24))
log("(a) f(oo) =", sp.factor(flim))
def nbar_dot(c, s, x):  # nbar.x := x^0 - n.xvec  (= -(1,n).x in mostly plus)
    return x[0] - c*x[1] - s*x[2]
x42 = tuple(x4[i] - x2[i] for i in range(3)); x24v = tuple(-v for v in x42)
closed = 2/((1 - (c1*c3 + s1*s3))*x24) * (u1*nbar_dot(c3, s3, x42) + u3*nbar_dot(c1, s1, x24v)
                                          + nbar_dot(c1, s1, x4)*nbar_dot(c3, s3, x2) - nbar_dot(c1, s1, x2)*nbar_dot(c3, s3, x4))
diff = sp.simplify(sp.expand(flim - closed.subs(sub)))
log("(a) limit r->oo of v'-u' minus C1_FINAL closed form (symbolic):", diff)
d11 = sp.simplify(sp.diff(flim, u1, 2)); d33 = sp.simplify(sp.diff(flim, u3, 2)); d13 = sp.simplify(sp.diff(flim, u1, u3))
log("(a) d2f/du1^2 =", d11, "; d2f/du3^2 =", d33, "; d2f/du1du3 =", d13)
# finite-r convergence, numeric
random.seed(3)
num = {u1: sp.Rational(37,100), u3: sp.Rational(-121,100), c1: sp.Rational(3,5), s1: sp.Rational(4,5), c3: sp.Rational(-7,25), s3: sp.Rational(24,25),
       x2[0]: sp.Rational(3,10), x2[1]: sp.Rational(-1,5), x2[2]: sp.Rational(1,2), x4[0]: sp.Rational(-11,10), x4[1]: sp.Rational(2,5), x4[2]: sp.Rational(1,10)}
lim_val = float(flim.subs(num))
num_f = {k: v for k, v in num.items() if k not in (s1, s3)}
for R in [1e3, 1e6, 1e9]:
    val = sp.N((vp - up).subs(num).subs(r, sp.Integer(int(R))), 30)
    log(f"(a) r = {R:.0e}: |f(r) - f(oo)| = {abs(float(val) - lim_val):.3e}")

# ---------- (b) ----------
def edot(X, Y):
    return -(X[0]*Y[1] + X[1]*Y[0])/2 + sum(X[2+i]*Y[2+i]*(-1 if i == 0 else 1) for i in range(3))
def Xpt(x):
    return (1, -x[0]**2 + x[1]**2 + x[2]**2, x[0], x[1], x[2])
def Zpol(x, z):
    return (0, 2*(-x[0]*z[0] + x[1]*z[1] + x[2]*z[2]), z[0], z[1], z[2])
def H(Zi, Xi, Zj, Xj):
    return -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))
def V(Zi, Xi, Xj, Xk):  # V_{i,jk}
    return (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)
zs2 = sp.symbols('z20 z21 z22'); zs4 = sp.symbols('z40 z41 z42')
X2s, X4s = Xpt(x2), Xpt(x4); Z2s, Z4s = Zpol(x2, zs2), Zpol(x4, zs4)
I = (0, 1, 0, 0, 0)
G = H(Z2s, X2s, Z4s, X4s) - 2*V(Z2s, X2s, I, X4s)*V(Z4s, X4s, I, X2s)
z2dz4 = -zs2[0]*zs4[0] + zs2[1]*zs4[1] + zs2[2]*zs4[2]
log("(b) G(X0 = I) - x24^2 (z2.z4) =", sp.simplify(sp.expand(G - sq(x2, x4)*z2dz4)))
# numeric: finite x0 far along a spatial direction (spacelike to x2, x4) -> ratio -> 1
for L in [1e2, 1e4, 1e6]:
    rnd = lambda: [random.uniform(-1, 1) for _ in range(3)]
    a2, a4 = rnd(), rnd()
    def null_pol():
        v = [random.uniform(-1, 1), random.uniform(-1, 1)]
        nv = (v[0]**2 + v[1]**2)**0.5
        return [1.0, v[0]/nv, v[1]/nv]   # null
    p2, p4 = null_pol(), null_pol()
    x0 = [0.0, L*0.6, L*0.8]
    X0 = Xpt(x0); X2n, X4n = Xpt(a2), Xpt(a4); Z2n, Z4n = Zpol(a2, p2), Zpol(a4, p4)
    Gn = H(Z2n, X2n, Z4n, X4n) - 2*V(Z2n, X2n, X0, X4n)*V(Z4n, X4n, X0, X2n)
    ref = (-(a2[0]-a4[0])**2 + (a2[1]-a4[1])**2 + (a2[2]-a4[2])**2)*(-p2[0]*p4[0] + p2[1]*p4[1] + p2[2]*p4[2])
    log(f"(b) x0 at spatial distance {L:.0e}: G/(x24^2 z2.z4) - 1 = {Gn/ref - 1:.3e}")
# invariance of H, V under (X4, Z4) -> -(X4, Z4) (the x4^+ reading)
X4m = tuple(-v for v in X4s); Z4m = tuple(-v for v in Z4s)
Gm = H(Z2s, X2s, Z4m, X4m) - 2*V(Z2s, X2s, I, X4m)*V(Z4m, X4m, I, X2s)
log("(b) G(-X4,-Z4) - G(X4,Z4) =", sp.simplify(sp.expand(Gm - G)))
open('/home/user/MAJOR/ads4/code/indep3/cache/i3_affine.txt', 'w').write("\n".join(out) + "\n")
