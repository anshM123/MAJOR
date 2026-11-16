"""
Three-point structures <J_1 J_2 O_{Delta,J}> (two conserved currents, Delta_1 = Delta_2 = 2), parity even.
 (a) the 5 words of the differential basis are linearly independent and transverse (gauge invariant in Z_1, Z_2);
 (b) conservation at x_1 and x_2: two-dimensional solution space, polynomial basis (v_1, v_2);
 (c) EXACT Lorentzian shadow map in the differential basis: the right structure of the shadow representation is
        R[c] = sum_w c_w(Delta,J)/b_w(Delta,J) D_w S3^S,     b_w = b^{(Delta_34 of the seed of w)}_{Delta,J};
     test whether  c_w(Delta,J) b^{(0)}/b_w = c_w(3-Delta,-1-J)  for each basis vector (shadow = (1/b^{(0)}) x Weyl reflection);
 (d) the identification with the Caron-Huot--Li helicity structures is in t1b_helicity_l1.py.
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
    # (d) the identification with CHL's helicity structures is done at exact rational points in chl_map.py / t1b_helicity_l1.py
    pickle.dump(dict(words=ws, vecs=vecs), open(os.path.join(CACHE, 'basis_l1.pkl'), 'wb'))
    open('out_t1_structures_l1.txt', 'w').write("\n".join(out) + "\n")
