"""Exact identification of Caron-Huot--Li's parity-even helicity structures T^{E,opp}, T^{E,same} (arXiv:2102.08160, eq. VVstcor)
inside the conserved <JJO> space spanned by (v_1, v_2) (t1_structures_l1), at exact rational (Delta, J):
     T^{E,h} = kappa(Delta,J) [ a_h v_1 + b_h v_2 ]    (kappa common: frame normalisation of M_ref).
Method: evaluate in CHL's frame (O_1(0), O_2(x), O(inf)) [embedding: P_inf = (0,1,0), Z_inf = (0,0,zeta)], expand on CHL's monomials B_V."""
import random, pickle, os, sys
import sympy as sp
from structures2 import rel_structures
from embed2 import Delta, J, CACHE

def frame_Q():
    rs, S, v5 = rel_structures(1)
    x = sp.Matrix(sp.symbols('X0:3')); e1 = sp.Matrix(sp.symbols('e1_0:3')); e2 = sp.Matrix(sp.symbols('e2_0:3')); e3 = sp.Matrix(sp.symbols('e3_0:3'))
    d = lambda a, b: (a.T*b)[0]
    fsub = {S['s12']: -d(x, x), S['s15']: -1, S['s25']: -1, S['z12']: d(e1, x), S['z15']: 0, S['z21']: -d(e2, x), S['z25']: 0,
            S['z51']: 0, S['z52']: d(e3, x), S['z1z2']: d(e1, e2), S['z1z5']: d(e1, e3), S['z2z5']: d(e2, e3)}
    BV0 = [d(e1, e2), d(e1, x)*d(e2, x)/d(x, x), d(e1, x)*d(e2, e3)/d(e3, x), d(e1, e3)*d(e2, x)/d(e3, x), d(e1, e3)*d(e2, e3)*d(x, x)/d(e3, x)**2]
    BV = [b*d(x, x) for b in BV0]
    return {w: Q.subs(fsub) for (w, Q) in rs}, BV, (x, e1, e2, e3)

TRIP = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (9, 40, 41)]
def points(syms, n=7, seed=11):
    random.seed(seed)
    x, e1, e2, e3 = syms
    pts = []
    for k in range(n):
        p = {}
        for i in range(3): p[x[i]] = sp.Rational(random.randint(-9, 9), random.randint(1, 5)) or 1
        for e in (e1, e2, e3):
            a, b, c = random.choice(TRIP); v = [sp.Integer(a), sp.Integer(random.choice([1, -1])*b), sp.I*c]; random.shuffle(v)
            for i in range(3): p[e[i]] = v[i]
        pts.append(p)
    return pts

def chl_vectors(Dv, Jv):
    n = (Dv - Jv - 2)/2; Jt = -1 - Jv
    T_opp = [2*(n - Jt), 2*(n - 1), 3*Jt - 4*n + 1, 3*Jt - 4*n + 1, (Jt**2 - (8*n + 1)*Jt + 8*n**2)/(2*n)]
    T_same = [2*(-n - Jv), 2*(n - 1), Jv, Jv, (Jv - 1)*Jv/(2*n)]
    return T_opp, T_same

def identify(Dv, Jv, Qs, BV, syms, ws, vecs):
    pts = points(syms)
    num = {Delta: Dv, J: Jv}
    Bmat = sp.Matrix([[sp.expand(b.subs(p)) for b in BV] for p in pts[:5]])
    Binv = Bmat.inv()
    cv = []
    for v in vecs:
        F = sum(c.subs(num)*Qs[w].subs(num) for c, w in zip(v, ws))
        vals = [sp.expand(F.subs(p)) for p in pts]
        k = Binv*sp.Matrix(vals[:5])
        k = k.applyfunc(lambda e: sp.nsimplify(sp.expand_complex(sp.simplify(e))))
        for p, val in zip(pts[5:], vals[5:]):
            assert sp.simplify(val - sum(k[i]*sp.expand(BV[i].subs(p)) for i in range(5))) == 0
        cv.append(k)
    res = {}
    for name, T in zip(('opp', 'same'), chl_vectors(Dv, Jv)):
        A = sp.Matrix.hstack(cv[0], cv[1]).applyfunc(lambda e: sp.nsimplify(sp.expand_complex(sp.expand(e))))
        Tm = sp.Matrix(T)
        # exact solve on two independent rows, then check all five
        for (i, j) in ((0, 1), (0, 4), (1, 4), (0, 2)):
            sub = A.extract([i, j], [0, 1])
            if sub.det() != 0: break
        sol = sub.LUsolve(Tm.extract([i, j], [0])).applyfunc(lambda e: sp.nsimplify(sp.expand_complex(sp.simplify(e))))
        chk = (A*sol - Tm).applyfunc(lambda e: sp.simplify(sp.expand_complex(e)))
        assert chk == sp.zeros(5, 1), (name, chk)
        res[name] = (sol[0], sol[1])
    return res

if __name__ == "__main__":
    ws, vecs = pickle.load(open(os.path.join(CACHE, 'consbasis_l1.pkl'), 'rb'))
    Qs, BV, syms = frame_Q()
    for (Dv, Jv) in ((sp.Rational(97, 10), 4), (sp.Rational(31, 3), 2), (sp.Rational(53, 4), 5), (sp.Rational(23, 2), 3)):
        r = identify(Dv, Jv, Qs, BV, syms, ws, vecs)
        Ap, Am = Dv - sp.Rational(3, 2), Jv + sp.Rational(1, 2); m2 = Ap**2 - Am**2
        print(Dv, Jv, {k: (a, b, sp.nsimplify(b/a)) for k, (a, b) in r.items()}, ' Am^2 m^2 =', Am**2*m2, ' Ap^2 m^2 =', Ap**2*m2, flush=True)
