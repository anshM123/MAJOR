"""
d = 3 conformal three-point structures <O(x4,z4) T(x1,z1) T(x3,z3)> (O: dimension Delta, spin J) in the
embedding formalism, and the solution of conservation + Bose symmetry.

Embedding: X = (X^+, X^-, x^mu) = (1, x^2, x^mu), Z = (0, 2 x.z, z^mu), X.Y = -(X^+Y^- + X^-Y^+)/2 + x.y,
metric (-,+,+).  Building blocks (Costa-Penedones-Poland-Rychkov):
  V1 = V_{1,34}, V3 = V_{3,41}, V4 = V_{4,13},  V_{i,jk} = (Z_i.X_j X_i.X_k - Z_i.X_k X_i.X_j)/(X_j.X_k)
  H_ij = -2[(Z_i.Z_j)(X_i.X_j) - (Z_i.X_j)(Z_j.X_i)],   X_ij = -2 X_i.X_j = x_ij^2
  eps13 = det[Z1,Z3,X1,X3,X4], eps14 = det[Z1,Z4,X1,X3,X4], eps34 = det[Z3,Z4,X1,X3,X4]   (5x5, components (+,-,0,1,2))
A monomial is  mono = (m1, m3, m4, n13, n14, n34, eps)  with eps in ('', '13', '14', '34'):
  parity even:  S = V1^m1 V3^m3 V4^m4 H13^n13 H14^n14 H34^n34 / (X13^{e13} X14^{e14} X34^{e34}),
                e13 = 5 - s, e14 = e34 = s,  s := (Delta + J)/2
  parity odd :  S = eps_ij X_ij * (same, with the Z-degrees at i and j lowered by one)
                    / (X13^{e13+1/2} X14^{e14+1/2} X34^{e34+1/2})
                (V_i and H_ij have X-degree one at i (and j), eps_ij has X-degree one at all three points, so the
                 homogeneity at the third point is restored by (X_ij/(X_ik X_jk))^{1/2}.)
Two-point normalisation: <O(x,z)O(0,z')> = (z.I(x).z')^J/(x^2)^Delta  (unit), <TT> = C_T (z.I.z')^2/x^6.
"""
import itertools, random
from fractions import Fraction
import sympy as sp

ETA = (-1, 1, 1)
def dot3(u, v):
    return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

def emb_X(x):
    return [1, dot3(x, x), x[0], x[1], x[2]]
def emb_Z(x, z):
    return [0, 2*dot3(x, z), z[0], z[1], z[2]]
def edot(A, B):
    return -(A[0]*B[1] + A[1]*B[0])/sp.Integer(2) + dot3(A[2:], B[2:])

def det5(rows):
    return sp.Matrix(rows).det(method='berkowitz')

def blocks_from_embedding(X1, Z1, X3, Z3, X4, Z4):
    V = lambda Zi, Xi, Xj, Xk: (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)
    H = lambda Zi, Xi, Zj, Xj: -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))
    return dict(
        V1=V(Z1, X1, X3, X4), V3=V(Z3, X3, X4, X1), V4=V(Z4, X4, X1, X3),
        H13=H(Z1, X1, Z3, X3), H14=H(Z1, X1, Z4, X4), H34=H(Z3, X3, Z4, X4),
        X13=-2*edot(X1, X3), X14=-2*edot(X1, X4), X34=-2*edot(X3, X4),
        eps13=det5([Z1, Z3, X1, X3, X4]), eps14=det5([Z1, Z4, X1, X3, X4]), eps34=det5([Z3, Z4, X1, X3, X4]))

def blocks(x1, z1, x3, z3, x4, z4):
    return blocks_from_embedding(emb_X(x1), emb_Z(x1, z1), emb_X(x3), emb_Z(x3, z3), emb_X(x4), emb_Z(x4, z4))

def mono_value(mono, B):
    m1, m3, m4, n13, n14, n34, eps = mono
    r = B['V1']**m1*B['V3']**m3*B['V4']**m4*B['H13']**n13*B['H14']**n14*B['H34']**n34
    if eps:
        r = r*B['eps'+eps]*B['X'+eps]
    return r

def monomials(J, parity):
    out = []
    if parity == 'even':
        specs = [('', 2, 2, J)]
    else:
        specs = [('13', 1, 1, J), ('14', 1, 2, J-1), ('34', 2, 1, J-1)]
    for eps, j1, j3, j4 in specs:
        if j4 < 0:
            continue
        for n13 in range(0, 3):
            for n14 in range(0, 3):
                for n34 in range(0, 3):
                    m1, m3, m4 = j1-n13-n14, j3-n13-n34, j4-n14-n34
                    if min(m1, m3, m4) >= 0:
                        out.append((m1, m3, m4, n13, n14, n34, eps))
    return out

# ------------------------------------------------------------------ sample points
NULLS = [(5, 3, 4), (5, 4, 3), (13, 5, 12), (13, 12, 5), (17, 8, 15), (17, 15, 8), (5, -3, 4), (5, 4, -3),
         (13, -5, 12), (25, 7, 24), (25, 24, -7), (29, 20, 21), (1, 1, 0), (1, 0, -1), (37, 12, 35), (41, 9, -40)]
def rnd_point(rng):
    return [sp.Rational(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(3)]
def rnd_null(rng):
    v = rng.choice(NULLS)
    return [sp.Integer(c) for c in v]

def independent_subset(monos, seed=1, npts=None):
    """Select a maximal linearly independent subset of the monomials (as functions, null polarisations)."""
    rng = random.Random(seed)
    npts = npts or (len(monos) + 6)
    rows = []
    for _ in range(npts):
        B = blocks(rnd_point(rng), rnd_null(rng), rnd_point(rng), rnd_null(rng), rnd_point(rng), rnd_null(rng))
        rows.append([mono_value(m, B) for m in monos])
    M = sp.Matrix(rows)
    _, piv = M.rref()
    return [monos[i] for i in piv], M.nullspace()

# ------------------------------------------------------------------ conservation + Bose symmetry
_x1 = sp.symbols('x1_0:3'); _z1 = sp.symbols('z1_0:3')
def divergence_rows(monos, parity, s, rng, n_conf=4):
    """Rows of the linear system  sum_a c_a Div_1[S_a] = 0, where Div = d/dx1^mu D_{z1}^mu (Todorov operator, d=3).
    The non-integer power prefactor is treated analytically.  s may be a sympy symbol or a number."""
    half = sp.Rational(1, 2) if parity == 'odd' else 0
    e13, e14 = 5 - s + half, s + half
    rows = []
    for _ in range(n_conf):
        x3, z3, x4, z4 = rnd_point(rng), rnd_null(rng), rnd_point(rng), rnd_null(rng)
        x1v = rnd_point(rng)
        B = blocks(list(_x1), list(_z1), x3, z3, x4, z4)
        sub = dict(zip(_x1, x1v))
        dlog = [-(e13*sp.diff(B['X13'], _x1[mu])/B['X13'] + e14*sp.diff(B['X14'], _x1[mu])/B['X14']) for mu in range(3)]
        cols = []
        for m in monos:
            R = mono_value(m, B)
            dR = [sp.diff(R, _z1[mu]) for mu in range(3)]
            box = sum(ETA[mu]*sp.diff(dR[mu], _z1[mu]) for mu in range(3))
            # D^mu R = (3/2) d^mu R - (1/2) z^mu box R   (R homogeneous of degree 2 in z1); upper index
            DR = [sp.Rational(3, 2)*ETA[mu]*dR[mu] - sp.Rational(1, 2)*_z1[mu]*box for mu in range(3)]
            div = sum(sp.diff(DR[mu], _x1[mu]) + dlog[mu]*DR[mu] for mu in range(3))
            div = div.subs(sub)
            comps = [sp.expand(sp.diff(div, _z1[nu])) for nu in range(3)]   # linear in z1
            cols.append(comps)
        for nu in range(3):
            rows.append([cols[a][nu] for a in range(len(monos))])
    return rows

def symmetry_rows(monos, rng, n_pts=None):
    """Rows of  sum_a c_a [S_a(1,3) - S_a(3,1)] = 0  (the power prefactor is symmetric because e14 = e34)."""
    rows = []
    for _ in range(n_pts or len(monos) + 4):
        x1, z1, x3, z3, x4, z4 = (rnd_point(rng), rnd_null(rng), rnd_point(rng), rnd_null(rng),
                                  rnd_point(rng), rnd_null(rng))
        B = blocks(x1, z1, x3, z3, x4, z4)
        Bs = blocks(x3, z3, x1, z1, x4, z4)
        rows.append([mono_value(m, B) - mono_value(m, Bs) for m in monos])
    return rows

def conserved_structures(J, parity, s=None, seed=7, verbose=False, symmetric=True):
    """Return (monos, basis): basis = list of coefficient vectors (sympy, rational functions of s) spanning the
    conserved, Bose-symmetric structures.  s: sympy symbol (generic Delta) or a number."""
    if s is None:
        s = sp.Symbol('s')
    rng = random.Random(seed)
    monos, rel = independent_subset(monomials(J, parity), seed=seed)
    rows = divergence_rows(monos, parity, s, rng)
    if symmetric:
        rows += symmetry_rows(monos, rng)
    M = sp.Matrix(rows)
    if M.free_symbols:
        from sympy.polys.matrices import DomainMatrix
        dM = DomainMatrix.from_Matrix(M).to_field()
        ns = dM.nullspace().to_Matrix()
        basis = [sp.simplify(ns.row(i).T) for i in range(ns.rows)]
        basis = [b for b in basis if any(e != 0 for e in b)]
    else:
        basis = M.nullspace()
    if verbose:
        print(f"J={J} {parity}: {len(monomials(J, parity))} monomials, {len(monos)} independent, "
              f"{len(basis)} conserved+symmetric")
    return monos, basis
