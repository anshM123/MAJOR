"""
Conserved parity-even three-point structures <O_1 O_2 O_5>, O_1, O_2 conserved of spin l (Delta_i = l + 1), O_5 = (Delta, J),
in the differential basis (embed2.word3), with conservation imposed at x_1 AND at x_2.

Words (operators applied right to left; see embed2):
   l = 1:  AC, AE, BC, BE, H            (5 = number of parity-even structures for (1,1,J), J >= 2)
   l = 2:  14 words  {AA,AB,BB} x {CC,CE,EE},  {A,B} x {C,E} x H,  HH   (one d = 3 linear relation among them)
Reference prefactor:  M_ref = S3(D1 + l, D2 + l) (all words have shifts summing to 2l);  W_w = M_ref * Qrel_w.
"""
import sympy as sp, itertools, random
from embed2 import word3, Delta, J, D1, D2, shifts, cached

eta = sp.diag(-1, 1, 1)

def words(l):
    if l == 1:
        return ['AC', 'AE', 'BC', 'BE', 'H']
    out = []
    for n in range(l + 1):                       # number of H's
        k = l - n
        P1 = [''.join(w) for w in itertools.combinations_with_replacement('AB', k)]
        P2 = [''.join(w) for w in itertools.combinations_with_replacement('CE', k)]
        for a in P1:
            for b in P2:
                out.append(a + b + 'H'*n)
    return out

def rel_structures(l, Dext=None):
    """list of (word, Qrel) with Qrel relative to M_ref = S3(D1+l, D2+l); Dext: substitute D1 = D2 = Dext (default l+1)."""
    Dext = l + 1 if Dext is None else Dext
    res = []
    for w in words(l):
        Q, S, v5, _ = word3(w)
        s1, s2 = shifts(w)
        assert s1 + s2 == 2*l
        h = sp.Rational(s1 - s2, 2)
        Q = Q.subs({D1: Dext, D2: Dext})
        res.append((w, sp.together(Q*S['s15']**(-h)*S['s25']**h)))
    return res, S, v5

def logMref(l, S, v5, Dext=None):
    Dext = l + 1 if Dext is None else Dext
    d1 = d2 = Dext + l
    al = (d1 + d2 - Delta + J)/2; be = (Delta + J + d1 - d2)/2; ga = (Delta + J - d1 + d2)/2
    return J*sp.log(v5) - al*sp.log(S['s12']) - be*sp.log(S['s15']) - ga*sp.log(S['s25'])

def _phys(S):
    """physical-space substitution with symbolic x1, x2 (for differentiation), numeric x5, n."""
    X1 = sp.Matrix(sp.symbols('x1_0:3')); X2 = sp.Matrix(sp.symbols('x2_0:3'))
    Z1 = sp.Matrix(sp.symbols('z1_0:3')); Z2 = sp.Matrix(sp.symbols('z2_0:3'))
    x5 = sp.Matrix([sp.Rational(9, 4), sp.Rational(-2, 3), sp.Rational(3, 7)]); n = sp.Matrix([5, 3, 4])
    def md(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
    sub = {S['s12']: -md(X1-X2, X1-X2), S['s15']: -md(X1-x5, X1-x5), S['s25']: -md(X2-x5, X2-x5),
           S['z12']: md(Z1, X2-X1), S['z15']: md(Z1, x5-X1), S['z21']: md(Z2, X1-X2), S['z25']: md(Z2, x5-X2),
           S['z51']: md(n, X1-x5), S['z52']: md(n, X2-x5), S['z1z2']: md(Z1, Z2), S['z1z5']: md(Z1, n), S['z2z5']: md(Z2, n)}
    return X1, X2, Z1, Z2, sub

def _div_spin1(F, logM, X, Z):
    Qmu = [sp.diff(F, Z[i]) for i in range(3)]
    return sum(eta[i, i]*(sp.diff(logM, X[i])*Qmu[i] + sp.diff(Qmu[i], X[i])) for i in range(3))

def _div_spin2(F, logM, X, Z):
    """divergence of the traceless symmetric tensor encoded by F(z) (quadratic in z; trace removed), returns vector in the free index."""
    Qmn = [[sp.diff(F, Z[i], Z[j])/2 for j in range(3)] for i in range(3)]
    tr = sum(eta[i, i]*Qmn[i][i] for i in range(3))
    T = [[Qmn[i][j] - sp.Rational(1, 3)*eta[i, j]*tr for j in range(3)] for i in range(3)]
    return [sum(eta[i, i]*(sp.diff(logM, X[i])*T[i][j] + sp.diff(T[i][j], X[i])) for i in range(3)) for j in range(3)]

def conservation_matrix(l, pts=None, Dext=None, Jval=None, Dval=None):
    """rows: conservation equations (at x1 and at x2) evaluated at rational points, columns: words."""
    rs, S, v5 = rel_structures(l, Dext)
    X1, X2, Z1, Z2, sub = _phys(S)
    lm = logMref(l, S, v5, Dext).subs(sub)
    num = {}
    if Jval is not None: num[J] = Jval
    if Dval is not None: num[Delta] = Dval
    lm = lm.subs(num)
    random.seed(7)
    pts = pts or [({X1[0]: sp.Rational(1, 3), X1[1]: sp.Rational(2, 7), X1[2]: sp.Rational(-3, 5),
                    X2[0]: sp.Rational(-7, 3), X2[1]: sp.Rational(1, 2), X2[2]: sp.Rational(-1, 5)},
                   {Z1[0]: sp.Rational(3, 2), Z1[1]: sp.Rational(1, 3), Z1[2]: sp.Rational(-2, 7),
                    Z2[0]: sp.Rational(2, 5), Z2[1]: sp.Rational(-3, 4), Z2[2]: sp.Rational(5, 3)})]
    rows = []
    for (xp, zp) in pts:
        cols = []
        for (w, Q) in rs:
            F = Q.subs(sub).subs(num)
            eqs = []
            for (X, Z, Zo) in ((X1, Z1, Z2), (X2, Z2, Z1)):
                if l == 1:
                    dv = [_div_spin1(F, lm, X, Z)]
                else:
                    dv = _div_spin2(F, lm, X, Z)
                for d in dv:
                    d = d.subs(xp)
                    # the divergence is a polynomial in the other polarisation: substitute a few numeric values
                    if l == 1:
                        for zz in ([1, 0, 0], [0, 1, 0], [0, 0, 1]):
                            eqs.append(d.subs({Zo[i]: zz[i] for i in range(3)}))
                    else:
                        for zz in ([1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], [1, 0, 1], [0, 1, 1]):
                            eqs.append(d.subs({Zo[i]: zz[i] for i in range(3)}))
            cols.append([sp.together(e) for e in eqs])
        for i in range(len(cols[0])):
            rows.append([cols[k][i] for k in range(len(rs))])
    return sp.Matrix(rows), [w for (w, _) in rs]
