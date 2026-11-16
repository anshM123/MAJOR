"""
Conserved parity-even three-point structures <O_1 phi_2 O_5> (O_1 = J: l=1, Delta_1=2;  O_1 = T: l=2, Delta_1=3)
in the differential basis, and their expansion in the standard basis.

Standard basis (embedding; P_i.P_j = s_ij/2):
   V_1 = [(Z1.P2)(P1.P5) - (Z1.P5)(P1.P2)]/(P2.P5),   V_5 = [(Z5.P1)(P2.P5) - (Z5.P2)(P1.P5)]/(P1.P2),
   H_15 = (Z1.Z5)(P1.P5) - (Z1.P5)(Z5.P1),
   S_n = V_1^{l-n} H_15^n V_5^{J-n} / (s12^{k12} s15^{k15} s25^{k25}),
   k12 = (Delta1+l+Delta2-Delta-J)/2, k15 = (Delta1+l+Delta+J-Delta2)/2, k25 = (Delta2+Delta+J-Delta1-l)/2.
Differential basis:  W_w = D-word(w) applied to  S3(Delta1+#a, Delta2+#b) = (V_5 s12/2)^J/(s12^al s15^be s25^ga).
"""
import sympy as sp, itertools
from embed import word3, Delta, J, D1, D2, cached

Dphi = sp.Symbol('Dphi')
eta = sp.diag(-1, 1, 1)
def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]

WORDS = {1: ['a', 'b'], 2: ['aa', 'ab', 'bb']}

def _coords():
    t, a, b = sp.symbols('t a b'); z = sp.symbols('zz0 zz1 zz2')
    x1 = sp.Matrix([t, a, b]); x2 = sp.Matrix([sp.Rational(-7, 3), sp.Rational(1, 2), sp.Rational(-1, 5)])
    x5 = sp.Matrix([sp.Rational(9, 4), sp.Rational(-2, 3), sp.Rational(3, 7)]); n = sp.Matrix([5, 3, 4])
    zz = sp.Matrix(z)
    sub = dict(s12=-mdot(x1-x2, x1-x2), s15=-mdot(x1-x5, x1-x5), s25=-mdot(x2-x5, x2-x5),
               z12=mdot(zz, x2-x1), z15=mdot(zz, x5-x1), z51=mdot(n, x1-x5), z52=mdot(n, x2-x5), z1z5=mdot(zz, n))
    return (t, a, b), z, sub

def rel_structure(l):
    """returns list of (word, Q_w * relative prefactor) with common prefactor M_ref = v5^J s12^{-al} (s15 s25)^{-..}"""
    res = []
    for w in WORDS[l]:
        Q, sy, V5 = word3(w)
        na, nb = w.count('a'), w.count('b')
        Q = Q.subs({D1: l + 1, D2: Dphi})
        # S3(D1+na, D2+nb) = v5^J s12^{-(al0 + l/2)} s15^{-(be0 + (na-nb)/2)} s25^{-(ga0 - (na-nb)/2)};
        # relative to M_ref := v5^J s12^{-(al0+l/2)} s15^{-(be0 + l/2)} s25^{-(ga0 + l/2)}:   s15^{nb} s25^{na}
        res.append((w, sp.together(Q*sy['s15']**nb*sy['s25']**na), sy, V5))
    return res

def logMref(l, sy, V5):
    al0 = (l + 1 + Dphi - Delta + J)/2; be0 = (Delta + J + l + 1 - Dphi)/2; ga0 = (Delta + J - l - 1 + Dphi)/2
    return J*sp.log(V5) - (al0 + sp.Rational(l, 2))*sp.log(sy['s12']) - (be0 + sp.Rational(l, 2))*sp.log(sy['s15']) - (ga0 + sp.Rational(l, 2))*sp.log(sy['s25'])

def conserved_coeffs(l):
    """coefficients c_w (rational functions of Delta, J, Dphi) such that sum_w c_w W_w is conserved at point 1."""
    (t, a, b), z, sub = _coords()
    rs = rel_structure(l)
    sy, V5 = rs[0][2], rs[0][3]
    symsub = {sy[k]: v for k, v in sub.items()}
    xs = (t, a, b)
    lm = logMref(l, sy, V5).subs(symsub)
    dlogM = [sp.diff(lm, q) for q in xs]
    pt = {t: sp.Rational(1, 3), a: sp.Rational(2, 7), b: sp.Rational(-3, 5)}
    cols = []
    for (w, Q, _, _) in rs:
        Qx = Q.subs(symsub)
        if l == 1:
            Qmu = [sp.diff(Qx, z[i]) for i in range(3)]                       # covariant components Q_mu
            div = sum(eta[i, i]*(dlogM[i]*Qmu[i] + sp.diff(Qmu[i], xs[i])) for i in range(3))
            eqs = [div]
        else:
            Qmn = [[sp.diff(Qx, z[i], z[j])/2 for j in range(3)] for i in range(3)]
            tr = sum(eta[i, i]*Qmn[i][i] for i in range(3))
            T = [[Qmn[i][j] - sp.Rational(1, 3)*eta[i, j]*tr for j in range(3)] for i in range(3)]
            eqs = [sum(eta[i, i]*(dlogM[i]*T[i][j] + sp.diff(T[i][j], xs[i])) for i in range(3)) for j in range(3)]
        cols.append([sp.factor(sp.together(e.subs(pt))) for e in eqs])
    M = sp.Matrix([[cols[w][i] for w in range(len(rs))] for i in range(len(cols[0]))])
    ns = M.nullspace()
    assert len(ns) == 1, len(ns)
    v = ns[0]
    v = sp.simplify(v/v[len(v)-1])
    # clear denominators
    den = sp.lcm([sp.fraction(sp.together(x))[1] for x in v])
    v = [sp.factor(sp.cancel(x*den)) for x in v]
    return dict(zip(WORDS[l], v))

def standard_expansion(l, coeffs=None):
    """expansion of W_w (and of the conserved combination) in the standard basis S_n."""
    rs = rel_structure(l)
    sy, V5 = rs[0][2], rs[0][3]
    v1 = (sy['z12']*sy['s15'] - sy['z15']*sy['s12'])/2
    h15 = sy['z1z5']*sy['s15']/2 - sy['z15']*sy['z51']
    V1s, V5s = 2*v1/sy['s25'], 2*V5/sy['s12']
    # W_w = M_ref * Qrel_w ;   S_n = V1^{l-n} H^n V5s^{J-n}/(s12^k12 s15^k15 s25^k25)
    # M_ref = v5^J s12^{-(al0+l/2)} s15^{-(be0+l/2)} s25^{-(ga0+l/2)},  v5^J = (V5s s12/2)^J
    # => M_ref/[V5s^J/(s12^k12 s15^k15 s25^k25)] = 2^{-J} s12^{J-al0-l/2+k12} s15^{k15-be0-l/2} s25^{k25-ga0-l/2}
    #    k12 - al0 = (l-1+... ) : compute explicitly
    al0 = (l + 1 + Dphi - Delta + J)/2; be0 = (Delta + J + l + 1 - Dphi)/2; ga0 = (Delta + J - l - 1 + Dphi)/2
    k12 = (2*l + 1 + Dphi - Delta - J)/2; k15 = (2*l + 1 + Delta + J - Dphi)/2; k25 = (Dphi + Delta + J - 2*l - 1)/2
    e12 = sp.simplify(J - al0 - sp.Rational(l, 2) + k12); e15 = sp.simplify(k15 - be0 - sp.Rational(l, 2)); e25 = sp.simplify(k25 - ga0 - sp.Rational(l, 2))
    pref = sy['s12']**e12*sy['s15']**e15*sy['s25']**e25          # times 2^{-J}
    lam = sp.symbols('lam0:%d' % (l+1))
    out = {}
    import random
    random.seed(1)
    for (w, Q, _, _) in rs:
        expr = sp.together(Q*pref - sum(lam[n]*V1s**(l-n)*h15**n*V5s**(-n) for n in range(l+1)))
        num = sp.numer(expr)
        # solve by evaluating at random rational points
        eqs = []
        for _ in range(l + 3):
            pt = {s: sp.Rational(random.randint(1, 40), random.randint(1, 9)) for s in sy.values()}
            eqs.append(num.subs(pt))
        sol = sp.solve(eqs, lam, dict=True)
        assert len(sol) == 1
        out[w] = [sp.factor(sol[0][q]) for q in lam]
        chk = sp.simplify(expr.subs(sol[0]))
        assert chk == 0, chk
    return out      # W_w = 2^{-J} * sum_n out[w][n] S_n

def get(l):
    return cached('struct_l%d' % l, lambda: (conserved_coeffs(l), standard_expansion(l)))

if __name__ == "__main__":
    for l in (1, 2):
        c, se = get(l)
        print("l =", l)
        print("  conserved combination  sum_w c_w W_w,  c_w =")
        for w, v in c.items(): print("    ", w, ":", v)
        print("  standard-basis expansion  W_w = 2^{-J} sum_n e_{w,n} S_n:")
        for w, v in se.items(): print("    ", w, ":", v)
        tot = [sp.factor(sp.simplify(sum(c[w]*se[w][n] for w in c))) for n in range(l+1)]
        print("  conserved structure = 2^{-J} sum_n L_n S_n,  L_n =", tot)
        print("  normalised to L_0 = 1:", [sp.factor(sp.simplify(x/tot[0])) for x in tot])
