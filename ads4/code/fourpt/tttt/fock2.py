"""
Exact mean-field OPE data of double-twist primaries [O_1 O_2]_{N,J} of two conserved operators of spin l (Delta_O = l + 1), from the
two-particle Fock space (method of SPINNING_MFT.md sect. 2, extended to two spinning particles; one-particle Gram matrices from
ads4/code/mft/fock.py, imported unchanged).

 * one-particle reproducing kernel of a conserved spin-l operator: fock.kernel_matrix('cons', a, l) (degree a in x, polarisation mod z^2);
 * a primary X of dimension 2 Delta_O + N, spin J has wavefunction q_X(x1, z1; x2, z2) = <0|X(inf, zeta) O(x1,z1) O(x2,z2)|0> = T(x1 - x2),
   a polynomial three-point structure conserved at x1 and x2; conversely every such polynomial is a primary;
 * two-particle inner product: tensor product of the one-particle ones:  <q, q'> = sum_{a1 + a2 = N} Tr(Q_P^T  M1_PP^{-1}  Q'_P  M2_PP^{-1})
   (P: pivot rows; range checked);
 * for degenerate primaries with structures S_a:  G_ab = <S_a(inf), S_b(inf)>, and the OPE-coefficient matrix for unit-normalised
   primaries is  P = (zeta.zetabar)^J G^{-1}  (contribution  sum_ab P_ab <..>_a <..>_b).
"""
import os, sys, functools
import sympy as sp
from flint import fmpq, fmpq_mat
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'mft'))
import fock                                    # noqa: E402
from fock import gram_inverse, pol_basis, reduce_null

U = sp.symbols('u0:3'); V = sp.symbols('v0:3'); A = sp.symbols('a0:3'); B = sp.symbols('b0:3')

def expand_wavefunction(T, X, l):
    """T: polynomial in X (= x1 - x2), A (pol of O_1), B (pol of O_2), complex rational coefficients.
    returns dict (a1, a2) -> {((m1, p1), (m2, p2)): coefficient}"""
    q = sp.expand(T.subs({X[i]: U[i] - V[i] for i in range(3)}, simultaneous=True))
    q = reduce_null(reduce_null(q, A, l), B, l)
    P = sp.Poly(q, *(U + V + A + B))
    data = {}
    for mon, c in P.terms():
        a1, a2 = sum(mon[0:3]), sum(mon[3:6])
        data.setdefault((a1, a2), {})[((mon[0:3], mon[6:9]), (mon[3:6], mon[9:12]))] = c
    return data

def _mat(entries, rows1, rows2, part):
    i1 = {r: i for i, r in enumerate(rows1)}; i2 = {r: i for i, r in enumerate(rows2)}
    Q = fmpq_mat(len(rows1), len(rows2))
    nz = False
    for (r1, r2), c in entries.items():
        v = sp.Rational(part(c))
        if v != 0:
            Q[i1[r1], i2[r2]] = fmpq(int(v.p), int(v.q)); nz = True
    return Q, nz

def _restrict(Q, piv1, piv2):
    R = fmpq_mat(len(piv1), len(piv2))
    for i, p in enumerate(piv1):
        for j, q in enumerate(piv2):
            R[i, j] = Q[p, q]
    return R

def _check_range(Q, M1, piv1, Minv1, M2, piv2, Minv2):
    """Q must lie in range(M1) (x) range(M2): Q = M1[:,P1] Y M2[P2,:] with Y = Minv1 Q_P Minv2"""
    QP = _restrict(Q, piv1, piv2)
    Y = Minv1*QP*Minv2
    M1P = fmpq_mat(M1.nrows(), len(piv1)); M2P = fmpq_mat(len(piv2), M2.ncols())
    for i in range(M1.nrows()):
        for j, p in enumerate(piv1): M1P[i, j] = M1[i, p]
    for i, p in enumerate(piv2):
        for j in range(M2.ncols()): M2P[i, j] = M2[p, j]
    assert M1P*Y*M2P == Q, "wavefunction not in the two-particle range (conservation violated?)"
    return QP

def bilinear(Ta, Tb, X, l):
    """<T_a, T_b> (sesquilinear: linear in T_a, antilinear in T_b)"""
    da = expand_wavefunction(Ta, X, l); db = expand_wavefunction(Tb, X, l)
    tot = sp.Integer(0)
    for key in set(da) | set(db):
        a1, a2 = key
        rows1, piv1, Minv1, M1 = gram_inverse('cons', a1, l)
        rows2, piv2, Minv2, M2 = gram_inverse('cons', a2, l)
        parts = {}
        for name, d in (('a', da), ('b', db)):
            for part, pn in ((sp.re, 'r'), (sp.im, 'i')):
                Q, nz = _mat(d.get(key, {}), rows1, rows2, part)
                parts[name + pn] = _check_range(Q, M1, piv1, Minv1, M2, piv2, Minv2) if nz else None
        def real_bil(Pa, Pb):
            if Pa is None or Pb is None: return sp.Integer(0)
            T = Pa.transpose()*Minv1*Pb*Minv2
            tr = sum((T[i, i] for i in range(T.nrows())), fmpq(0))
            return sp.Rational(int(tr.p), int(tr.q))
        re_ = real_bil(parts['ar'], parts['br']) + real_bil(parts['ai'], parts['bi'])
        im_ = real_bil(parts['ai'], parts['br']) - real_bil(parts['ar'], parts['bi'])
        tot += re_ + sp.I*im_
    return sp.nsimplify(tot)
