"""
Three-point structures <J_1 J_2 O_{Delta,J}> (two conserved currents, Delta_1 = Delta_2 = 2), parity even.
 (a) the 5 words of the differential basis are linearly independent and transverse (gauge invariant in Z_1, Z_2);
 (b) conservation at x_1 and x_2: two-dimensional solution space, polynomial basis (v_1, v_2);
 (c) EXACT Lorentzian shadow map in the differential basis: the right structure of the shadow representation is
        R[c] = sum_w c_w(Delta,J)/b_w(Delta,J) D_w S3^S,     b_w = b^{(Delta_34 of the seed of w)}_{Delta,J};
     test whether  c_w(Delta,J) b^{(0)}/b_w = c_w(3-Delta,-1-J)  for each basis vector (shadow = (1/b^{(0)}) x Weyl reflection);
 (d) identification with the Caron-Huot--Li parity-even helicity structures T^{E,opp}, T^{E,same} (arXiv:2102.08160, eq. VVstcor),
     by evaluating the differential-basis structures in CHL's frame (O_1(0), O_2(x), O(inf)) and expanding in their monomials B_V.
Output: out_t1_structures_l1.txt;  cache/basis_l1.pkl (v_1, v_2, CHL map).
"""
import sys, random, pickle, os
import sympy as sp
from structures2 import rel_structures, conservation_matrix, words
from embed2 import Delta, J, cached, CACHE
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)

def conserved_basis(l):
    def build():
        M, ws = conservation_matrix(l)
        ns = M.nullspace(simplify=True)
        vecs = []
        for v in ns:
            den = sp.lcm([sp.fraction(sp.together(x))[1] for x in v])
            vecs.append([sp.factor(sp.cancel(x*den)) for x in v])
        return ws, vecs
    return cached('consbasis_l%d' % l, build)

def b_ratio(D34):
    """b^{(D34)}/b^{(0)} for even integer D34, as a rational function of (Delta, J)"""
    k = abs(int(D34))//2
    assert 2*k == abs(int(D34))
    u, w = (J + Delta)/2, (Delta - 1 - J)/2
    # Gamma(u+k) Gamma(u-k) / Gamma(u)^2 = prod_{i<k} (u+i) / prod_{1<=i<=k} (u-i)
    r = sp.Integer(1)
    for v in (u, w):
        for i in range(k):
            r *= (v + i)/(v - i - 1)
    return sp.factor(r)

def seed_D34(word, l):
    from embed2 import shifts
    s1, s2 = shifts(word)
    return s1 - s2

if __name__ == "__main__":
    l = 1
    pr("# <J J O> parity-even structures in the differential basis (fourpt/t1_structures_l1.py)")
    rs, S, v5 = rel_structures(l)
    ws = [w for (w, _) in rs]
    random.seed(3)
    rows = []
    for k in range(8):
        pt = {s: sp.Rational(random.randint(1, 30), random.randint(1, 7)) for s in S.values()}
        pt[Delta] = sp.Rational(97, 10); pt[J] = 4
        rows.append([Q.subs(pt) for (w, Q) in rs])
    pr("(a) words:", ws, " rank at random kinematics (Delta, J) = (9.7, 4):", sp.Matrix(rows).rank())
    lam = sp.Symbol('lam')
    g1 = {S['z12']: S['z12'] + lam*S['s12']/2, S['z15']: S['z15'] + lam*S['s15']/2, S['z1z2']: S['z1z2'] + lam*S['z21'], S['z1z5']: S['z1z5'] + lam*S['z51']}
    g2 = {S['z21']: S['z21'] + lam*S['s12']/2, S['z25']: S['z25'] + lam*S['s25']/2, S['z1z2']: S['z1z2'] + lam*S['z12'], S['z2z5']: S['z2z5'] + lam*S['z52']}
    gi = [(sp.simplify(Q.subs(g1, simultaneous=True) - Q) == 0 and sp.simplify(Q.subs(g2, simultaneous=True) - Q) == 0) for (w, Q) in rs]
    pr("    transverse under Z_1 -> Z_1 + lam P_1 and Z_2 -> Z_2 + lam P_2 (exact):", gi)
    ws2, vecs = conserved_basis(l)
    assert ws2 == ws
    pr("(b) conserved at x_1 and x_2 (exact, symbolic in Delta, J): solution space of dimension", len(vecs))
    for i, v in enumerate(vecs):
        pr("    v_%d =" % (i+1), dict(zip(ws, v)))
    Ap, Am = sp.symbols('A_p A_m')
    weyl_sub = {Delta: Ap + sp.Rational(3, 2), J: Am - sp.Rational(1, 2)}
    pr("    in A_+ = Delta - 3/2, A_- = J + 1/2:")
    for i, v in enumerate(vecs):
        pr("    v_%d =" % (i+1), [sp.factor(sp.expand(x.subs(weyl_sub))) for x in v])
    pr("(c) Lorentzian shadow map  c_w -> c_w b^{(0)}/b_w  versus Weyl reflection c_w(3-Delta,-1-J):")
    refl = {Delta: 3 - Delta, J: -1 - J}
    for w in ws:
        pr("    word %s: seed Delta_34 = %d,  b_w/b^(0) = %s" % (w, seed_D34(w, l), b_ratio(seed_D34(w, l))))
    for i, v in enumerate(vecs):
        img = [sp.factor(sp.cancel(x/b_ratio(seed_D34(w, l)))) for x, w in zip(v, ws)]
        rv = [sp.factor(x.subs(refl, simultaneous=True)) for x in v]
        diff = [sp.simplify(a - b) for a, b in zip(img, rv)]
        pr("    v_%d: image - reflected =" % (i+1), diff)
    # (d) CHL frame
    x = sp.Matrix(sp.symbols('X0:3')); e1 = sp.Matrix(sp.symbols('e1_0:3')); e2 = sp.Matrix(sp.symbols('e2_0:3')); e3 = sp.Matrix(sp.symbols('e3_0:3'))
    d = lambda a, b: (a.T*b)[0]
    fsub = {S['s12']: -d(x, x), S['s15']: -1, S['s25']: -1, S['z12']: d(e1, x), S['z15']: 0, S['z21']: -d(e2, x), S['z25']: 0,
            S['z51']: 0, S['z52']: d(e3, x), S['z1z2']: d(e1, e2), S['z1z5']: d(e1, e3), S['z2z5']: d(e2, e3)}
    BV0 = [d(e1, e2), d(e1, x)*d(e2, x)/d(x, x), d(e1, x)*d(e2, e3)/d(e3, x), d(e1, e3)*d(e2, x)/d(e3, x), d(e1, e3)*d(e2, e3)*d(x, x)/d(e3, x)**2]
    BV = [b_*d(x, x) for b_ in BV0]          # the differential-basis prefactor M_ref carries an extra 1/x^2 relative to CHL's
    Jn = sp.Symbol('Jn')
    nsym = (Delta - J - 2)/2; Jt = -1 - J
    T_opp = [2*(nsym - Jt), 2*(nsym - 1), 3*Jt - 4*nsym + 1, 3*Jt - 4*nsym + 1, (Jt**2 - (8*nsym + 1)*Jt + 8*nsym**2)/(2*nsym)]
    T_same = [2*(-nsym - J), 2*(nsym - 1), J, J, (J - 1)*J/(2*nsym)]
    random.seed(11)
    def numvec():
        return sp.Matrix([sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(3)])
    def null():
        # rational null vector (Euclidean null vectors are complex; use a complex rational one)
        a, b = sp.Rational(random.randint(1, 9), random.randint(1, 5)), sp.Rational(random.randint(1, 9), random.randint(1, 5))
        return sp.Matrix([a, b, sp.I*sp.sqrt(a*a + b*b)])
    # coefficients of each conserved structure in B_V, symbolic in Delta, J: solve with numeric kinematics
    pts = []
    for k in range(7):
        pts.append({**{x[i]: v for i, v in enumerate(numvec())}, **{e1[i]: v for i, v in enumerate(null())},
                    **{e2[i]: v for i, v in enumerate(null())}, **{e3[i]: v for i, v in enumerate(null())}})
    Qs = {w: sp.together(Q.subs(fsub)) for (w, Q) in rs}
    coeffs = []
    for v in vecs:
        F = sum(c*Qs[w] for c, w in zip(v, ws))
        # F = (common) * sum_i k_i BV_i : the common factor is a power of (e3.x) and x^2: J-independent here since M_ref is factored out
        # find the power: F is homogeneous; divide by (e3.x)^0 -- test by solving with an extra overall unknown is not needed:
        ks = sp.symbols('k0:5')
        eqs = [sp.expand((F - sum(ks[i]*BV[i] for i in range(5))).subs(p)) for p in pts]
        sol = sp.solve(eqs, ks, dict=True)
        assert len(sol) == 1, sol
        coeffs.append([sp.factor(sol[0][q]) for q in ks])
    pr("(d) conserved structures in CHL's frame (O_1(0), O_2(x), O(inf)), coefficients on the monomials B_V of CHL (VVstcor),")
    pr("    up to the common factor (e3.x)^J |x|^{..} of M_ref:")
    for i, c in enumerate(coeffs):
        pr("    v_%d -> " % (i+1), c)
    # T_opp, T_same as combinations of v_1, v_2 (up to a common factor)
    a, b, c_, d_ = sp.symbols('a b c d')
    res = {}
    for name, T in (('T_opp', T_opp), ('T_same', T_same)):
        eqs = [a*coeffs[0][i] + b*coeffs[1][i] - T[i] for i in range(5)]
        sol = sp.solve(eqs, [a, b], dict=True)
        pr("    %s = a v_1 + b v_2 (times the common factor):" % name, "no solution" if not sol else {k: sp.factor(v) for k, v in sol[0].items()})
        res[name] = sol[0] if sol else None
    pickle.dump(dict(words=ws, vecs=vecs, coeffs_BV=coeffs, chl=res), open(os.path.join(CACHE, 'basis_l1.pkl'), 'wb'))
    open('out_t1_structures_l1.txt', 'w').write("\n".join(out) + "\n")
