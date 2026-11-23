"""
p2spin/p2_core.py (P2_SPINNING.md): exact building blocks for <T T O_{Delta,J}> in d = 3.

Conventions are those of repro_engine/tto_conserved.py and repro_engine/integrand.py (read, not imported here, to keep this module
dependency-free; tto_conserved.py is imported read-only by s1_structures.py for the J = 0 validation):
  * 3-vectors in light-cone components (u, v, y), x.y = -(u v' + v u')/2 + y y'  (signature (-,+,+));
  * embedding P = (1, x^2, x), Z = (0, 2 x.z, z), P.Q = -(P+ Q- + P- Q+)/2 + p.q,  X_ij = -2 P_i.P_j = (x_i - x_j)^2;
  * V_{i,jk} = [(Z_i.P_j)(P_i.P_k) - (Z_i.P_k)(P_i.P_j)]/(P_j.P_k),   H_ij = -2[(Z_i.Z_j)(P_i.P_j) - (Z_i.P_j)(Z_j.P_i)];
  * eps(A,B,C,D,E) = det of the five embedding vectors in components (P+, P-, u, v, y) (a fixed constant times the invariant eps).
CPPR structures (1107.3554 / 1109.6321), T at 1, 2 (Delta = 3, l = 2), O at 5 (Delta, J):
  even:  V1^{2-a-b} V2^{2-a-c} V5^{J-b-c} H12^a H15^b H25^c / (X12^{(10-Delta-J)/2} X15^{(Delta+J)/2} X25^{(Delta+J)/2}),
         V1 = V_{1,25}, V2 = V_{2,51}, V5 = V_{5,12};  14 monomials (a,b,c), a+b <= 2, a+c <= 2;
  odd:   e12 x deg(1,1,J), e15 x deg(1,2,J-1), e25 x deg(2,1,J-1) even monomials, with
         e12 = eps(Z1,Z2,P1,P2,P5) X12, e15 = eps(Z1,Z5,P1,P2,P5) X15, e25 = eps(Z2,Z5,P1,P2,P5) X25, and the common factor
         1/sqrt(X12 X15 X25) moved into the X-prefactor (so every odd monomial has the homogeneity of an even one); 21 monomials.
'N' below is always the monomial divided by V5^J (rational in the points); the factor V5^J X^{-alpha} is carried as exp(L).

Exact arithmetic: python Fractions; first derivatives in x1 by forward-mode dual numbers (Dual), exact.
"""
from fractions import Fraction as Fr
import random

# ---------------------------------------------------------------- dual numbers (value + gradient w.r.t. x1, 3 components)


class Dual:
    __slots__ = ('v', 'g')

    def __init__(self, v, g=None):
        self.v = v
        self.g = g if g is not None else (0, 0, 0)

    @staticmethod
    def lift(a):
        return a if isinstance(a, Dual) else Dual(a)

    def __add__(self, o):
        o = Dual.lift(o); return Dual(self.v + o.v, tuple(a + b for a, b in zip(self.g, o.g)))
    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v, tuple(-a for a in self.g))

    def __sub__(self, o):
        return self + (-Dual.lift(o))

    def __rsub__(self, o):
        return Dual.lift(o) - self

    def __mul__(self, o):
        o = Dual.lift(o)
        return Dual(self.v*o.v, tuple(self.v*b + o.v*a for a, b in zip(self.g, o.g)))
    __rmul__ = __mul__

    def inv(self):
        iv = 1/Fr(self.v) if not isinstance(self.v, float) else 1.0/self.v
        return Dual(iv, tuple(-a*iv*iv for a in self.g))

    def __truediv__(self, o):
        return self*Dual.lift(o).inv()

    def __rtruediv__(self, o):
        return Dual.lift(o)*self.inv()

    def __pow__(self, n):
        assert isinstance(n, int)
        if n < 0:
            return self.inv()**(-n)
        r = Dual(1)
        for _ in range(n):
            r = r*self
        return r

    def dlog(self):
        return tuple(a/self.v for a in self.g)


def V0(v):
    return v.v if isinstance(v, Dual) else v


# ---------------------------------------------------------------- geometry

def dot(a, b):
    return -(a[0]*b[1] + a[1]*b[0])/2 + a[2]*b[2]


def EX(x):
    return (1, dot(x, x), x[0], x[1], x[2])


def EZ(x, z):
    return (0, 2*dot(x, z), z[0], z[1], z[2])


def edot(A, B):
    return -(A[0]*B[1] + A[1]*B[0])/2 + dot(A[2:], B[2:])


def Vf(xi, zi, xj, xk):
    Xi, Zi, Xj, Xk = EX(xi), EZ(xi, zi), EX(xj), EX(xk)
    return (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)


def Hf(xi, zi, xj, zj):
    Xi, Zi, Xj, Zj = EX(xi), EZ(xi, zi), EX(xj), EZ(xj, zj)
    return -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))


def Xd(xi, xj):
    return -2*edot(EX(xi), EX(xj))


def det(M):
    """exact determinant (Laplace for n <= 5 is fine; Dual-compatible)."""
    n = len(M)
    if n == 1:
        return M[0][0]
    s = 0
    for j in range(n):
        if (isinstance(M[0][j], (int, Fr)) and M[0][j] == 0):
            continue
        minor = [row[:j] + row[j + 1:] for row in M[1:]]
        term = M[0][j]*det(minor)
        s = s + term if j % 2 == 0 else s - term
    return s


def eps5(A, B, C, D, E):
    return det([list(A), list(B), list(C), list(D), list(E)])


# ---------------------------------------------------------------- monomial lists

EVEN = [(a, b, c) for a in range(3) for b in range(3) for c in range(3) if a + b <= 2 and a + c <= 2]          # 14
ODD = ([('12', a, b, c) for a in range(2) for b in range(2) for c in range(2) if a + b <= 1 and a + c <= 1] +   # 5
       [('15', a, b, c) for a in range(2) for b in range(2) for c in range(3) if a + b <= 1 and a + c <= 2] +   # 8
       [('25', a, b, c) for a in range(2) for b in range(3) for c in range(2) if a + b <= 2 and a + c <= 1])    # 8


def blocks(x1, z1, x2, z2, x5, z5):
    """all building blocks at one configuration (works with Fractions, floats, Duals, sympy)."""
    P1, P2, P5 = EX(x1), EX(x2), EX(x5)
    Z1, Z2, Z5 = EZ(x1, z1), EZ(x2, z2), EZ(x5, z5)
    X12, X15, X25 = Xd(x1, x2), Xd(x1, x5), Xd(x2, x5)
    B = dict(V1=Vf(x1, z1, x2, x5), V2=Vf(x2, z2, x5, x1), V5=Vf(x5, z5, x1, x2),
             H12=Hf(x1, z1, x2, z2), H15=Hf(x1, z1, x5, z5), H25=Hf(x2, z2, x5, z5),
             X12=X12, X15=X15, X25=X25)
    B['e12'] = eps5(Z1, Z2, P1, P2, P5)*X12
    B['e15'] = eps5(Z1, Z5, P1, P2, P5)*X15
    B['e25'] = eps5(Z2, Z5, P1, P2, P5)*X25
    return B


def mono_even(B, abc):
    a, b, c = abc
    return B['V1']**(2 - a - b)*B['V2']**(2 - a - c)*B['V5']**(-b - c)*B['H12']**a*B['H15']**b*B['H25']**c


def mono_odd(B, t):
    typ, a, b, c = t
    if typ == '12':
        return B['e12']*B['V1']**(1 - a - b)*B['V2']**(1 - a - c)*B['V5']**(-b - c)*B['H12']**a*B['H15']**b*B['H25']**c
    if typ == '15':
        return B['e15']*B['V1']**(1 - a - b)*B['V2']**(2 - a - c)*B['V5']**(-1 - b - c)*B['H12']**a*B['H15']**b*B['H25']**c
    return B['e25']*B['V1']**(2 - a - b)*B['V2']**(1 - a - c)*B['V5']**(-1 - b - c)*B['H12']**a*B['H15']**b*B['H25']**c


# z1-factorisation: each monomial = (z1-independent scalar) x (linear form) x (linear form)
def mono_even_split(abc):
    a, b, c = abc
    lin = ['V1']*(2 - a - b) + ['H12']*a + ['H15']*b
    sc = lambda B: B['V2']**(2 - a - c)*B['V5']**(-b - c)*B['H25']**c
    return lin, sc


def mono_odd_split(t):
    typ, a, b, c = t
    if typ == '12':
        lin = ['e12'] + ['V1']*(1 - a - b) + ['H12']*a + ['H15']*b
        sc = lambda B: B['V2']**(1 - a - c)*B['V5']**(-b - c)*B['H25']**c
    elif typ == '15':
        lin = ['e15'] + ['V1']*(1 - a - b) + ['H12']*a + ['H15']*b
        sc = lambda B: B['V2']**(2 - a - c)*B['V5']**(-1 - b - c)*B['H25']**c
    else:
        lin = ['V1']*(2 - a - b) + ['H12']*a + ['H15']*b
        sc = lambda B: B['e25']*B['V2']**(1 - a - c)*B['V5']**(-1 - b - c)*B['H25']**c
    assert len(lin) == 2
    return lin, sc


# ---------------------------------------------------------------- random rational kinematics

def rnd(rng, lo=-9, hi=9, den=5):
    return Fr(rng.randint(lo, hi), rng.randint(1, den))


def rnull(rng):
    while True:
        u, y = rnd(rng), rnd(rng)
        if u != 0:
            return (u, y*y/u, y)


def rpoint(rng):
    return tuple(rnd(rng) for _ in range(3))


GL = ((0, Fr(-1, 2), 0), (Fr(-1, 2), 0, 0), (0, 0, 1))
GI = ((0, -2, 0), (-2, 0, 0), (0, 0, 1))
