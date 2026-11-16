"""
Parity-odd structure for l = 1 (checks used in SPINNING_MFT.md sect. 4.4):
 (i)   D~ S3 = -J (Delta-1) eps(Z1,P1,Z5,P2,P5) S3 / v5,   D~ = eps(Z1, P1, d/dP1, P2, d/dP2)   (numerical point, exact algebra)
 (ii)  infinity frame (Euclidean):  eps(Z1,P1,Z5,P2,P5) -> -L^2 eps3(z1, zeta, x12),  v5 -> -(1/2) L^2 zeta.x12
       for x5 = L e, z5 = I(e) zeta, L -> infinity
 (iii) heavy frame of SPINNING_E2E (Lorentzian; x1 = -y, x2 = x/x^2, x5 = e, z5 = n = (1,1,0); x = X/m, y = Y/m):
       eps/v5 -> -2 z1.e_perp   (and the same at the right vertex x3 = y, x4 = -x/x^2)
eps(A,...,E) := det of the upper components (P+, P-, x^0, x^1, x^2); P.Q = -(P+Q- + P-Q+)/2 + x.y.
"""
import itertools, sympy as sp
from sympy import LeviCivita
R = sp.Rational

def emb(sig):
    eta = sp.diag(sp.Matrix([[0, -R(1, 2)], [-R(1, 2), 0]]), sp.diag(*sig))
    d = lambda A, B: (A.T*eta*B)[0]
    P = lambda x: sp.Matrix([1, sum(sig[i]*x[i]**2 for i in range(3))] + list(x))
    Zv = lambda x, z: sp.Matrix([0, 2*sum(sig[i]*x[i]*z[i] for i in range(3))] + list(z))
    return eta, d, P, Zv

def eps5(*vs): return sp.Matrix.hstack(*vs).det()

def check_i():
    eta, d, P, Zv = emb((-1, 1, 1))
    P1 = sp.Matrix(sp.symbols('p0:5')); P2 = sp.Matrix(sp.symbols('q0:5'))
    x1 = [R(1, 3), R(-2, 5), R(1, 7)]; x2 = [R(-3, 2), R(1, 4), R(2, 9)]; x5 = [R(5, 3), R(1, 2), R(-1, 3)]
    z1 = [R(2, 3), R(5, 7), R(-1, 4)]; n = [1, R(3, 5), R(4, 5)]
    P5 = P(x5); Z5 = Zv(x5, n); Z1 = Zv(x1, z1)
    for (J, Dl, D2) in [(3, R(73, 10), R(7, 4)), (5, R(121, 10), R(5, 2))]:
        D1 = 2
        al = (D1 + D2 - Dl + J)/2; be = (Dl + J + D1 - D2)/2; ga = (Dl + J - D1 + D2)/2
        s = lambda A, B: 2*d(A, B)
        v5 = d(Z5, P1)*d(P2, P5) - d(Z5, P2)*d(P1, P5)
        S3 = v5**J/((-s(P1, P2))**al*(-s(P1, P5))**be*(-s(P2, P5))**ga)
        ei = eta.inv()
        g1 = ei*sp.Matrix([sp.diff(S3, P1[i]) for i in range(5)])
        sub = {P1[i]: P(x1)[i] for i in range(5)}; sub.update({P2[i]: P(x2)[i] for i in range(5)})
        G2 = [[(ei*sp.Matrix([sp.diff(g1[c], P2[i]) for i in range(5)]))[e].subs(sub) for e in range(5)] for c in range(5)]
        P1v, P2v = P(x1), P(x2)
        tot = sum(LeviCivita(*p)*Z1[p[0]]*P1v[p[1]]*G2[p[2]][p[4]]*P2v[p[3]] for p in itertools.permutations(range(5)))
        rhs = (-J*(Dl - 1)*eps5(Z1, P1v, Z5, P2v, P5)*S3/v5).subs(sub)
        print("(i)  J=%d Delta=%s Dphi=%s:  D~S3 / [-J(Delta-1) eps S3/v5] = %s" % (J, Dl, D2, sp.N(tot/rhs, 20)))

def check_ii():
    eta, d, P, Zv = emb((1, 1, 1))
    L = sp.Symbol('L')
    x1 = sp.Matrix([R(1, 3), R(-2, 5), R(3, 7)]); x2 = sp.Matrix([R(-1, 2), R(1, 4), R(2, 9)]); z1 = sp.Matrix([R(2, 3), R(5, 7), R(-1, 4)])
    zeta = sp.Matrix([1, sp.I, 0]); e = sp.Matrix([1, 0, 0]); z5 = (sp.eye(3) - 2*e*e.T)*zeta
    E = sp.Poly(sp.expand(eps5(Zv(x1, z1), P(x1), Zv(L*e, z5), P(x2), P(L*e))), L)
    v5 = sp.Poly(sp.expand(d(Zv(L*e, z5), P(x1))*d(P(x2), P(L*e)) - d(Zv(L*e, z5), P(x2))*d(P(x1), P(L*e))), L)
    x = x1 - x2
    print("(ii) eps: degree %d in L, leading/eps3(z1,zeta,x) = %s;  v5: degree %d, leading/(zeta.x) = %s" % (
        E.degree(), sp.nsimplify(sp.expand(E.coeffs()[0]/sp.Matrix.hstack(z1, zeta, x).det())), v5.degree(),
        sp.nsimplify(sp.expand(v5.coeffs()[0]/(zeta.T*x)[0]))))

def check_iii():
    eta, d, P, Zv = emb((-1, 1, 1))
    mdot = lambda a, b: -a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    X = sp.Matrix([R(13, 10), R(4, 10), R(-5, 10)]); Y = sp.Matrix([R(9, 10), R(-3, 10), R(2, 10)])
    zs = sp.symbols('z0:3'); z = sp.Matrix(zs)
    e = sp.Matrix([1, 0, 0]); n = sp.Matrix([1, 1, 0])
    for side in ('left', 'right'):
        for m in (10, 100, 1000, 10000):
            x = X/m; y = Y/m
            xa, xb = (-y, x/mdot(x, x)) if side == 'left' else (y, -x/mdot(x, x))
            E = eps5(Zv(xa, z), P(xa), Zv(e, n), P(xb), P(e))
            v5 = d(Zv(e, n), P(xa))*d(P(xb), P(e)) - d(Zv(e, n), P(xb))*d(P(xa), P(e))
            r = sp.expand(E/v5)
            print("(iii) %s vertex m=%5d: eps/v5 = %s z^t + %s z^1 + %s z^2   (-> -2 z.e_perp)" % (
                side, m, sp.N(r.coeff(zs[0]), 6), sp.N(r.coeff(zs[1]), 6), sp.N(r.coeff(zs[2]), 8)))

if __name__ == "__main__":
    check_i(); check_ii(); check_iii()
