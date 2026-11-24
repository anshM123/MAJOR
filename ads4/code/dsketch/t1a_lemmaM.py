"""
t1a_lemmaM (DSKETCH_UP.md §1.4): mode completeness for free self-dual spin-2 fields at fixed boundary momentum K, and the
boundary projection b_+ of the explicit ASD second-order solution.

(i) Flat R^4 = (x1, x2, x3, z), eps_{x1 x2 x3 z} = +1.  A self-dual Weyl-type tensor is Phi = sum_ab phi_ab S^a (x) S^b with S^a the
    three SD 2-forms and phi symmetric traceless (5 components).  For Phi = phi exp(i K.x + lam z) the divergence
    d^mu Phi_{mu nu rho sigma} = 0 is M(lam) phi = 0.  The solution space of the ODE system at fixed K has dimension
    deg gcd(5x5 minors of M(lam)); we find gcd = lam^2 - |K|^2 (dimension 2, lam = -|K| decaying, +|K| growing).
    By conformal invariance of the spin-2 equation (Penrose-Rindler I, 5.7: psi_AdS = z psi_flat for g_AdS = z^-2 delta) the same
    holds on AdS4: every free SD field with x-dependence e^{iK.x} that decays as z -> oo is a multiple of ONE field.
(ii) That field is W+ of the regular label-(+) graviton mode at K: in the conformal frame z^-1 W+[reg +] = e^{-|K| z} x (constant
    tensor in ker M(-|K|)); the irregular label-(-) mode gives the e^{+|K| z} solution.  Checked exactly for K = (4,0,3) and (3,1,4).
(iii) b_+ := label-(+) TT projection of the boundary value Gam_p(z=0) of t1a_solve2's ASD solution (gauge-invariant under boundary
    diffeomorphisms and Weyl rescalings).  Phi_34 = - b_+ W+[reg + mode at K] (normalised as in mode()).  Non-zero => the remnant is
    non-zero, as C1_VALUE §3.3 found from JJ21.
Writes cache/t1a_lemmaM.txt
"""
import os, sys, pickle, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from sympy import *  # noqa: F401,F403
from asd_alg import Alg, z, N, christoffel, riemann_up, ricci, inverse_conformal, mode, helicity_basis

lines = []
def log(x=''):
    print(x, flush=True); lines.append(str(x))

lam = sp.symbols('lam')
ok = True
# SD 2-forms of flat R^4 with eps_{0123} = +1 (index order x1 x2 x3 z): eigenvectors of * with eigenvalue +1
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def star_form(F):
    out = sp.zeros(4, 4)
    for a, b in itertools.product(range(4), repeat=2):
        out[a, b] = sum(sp.Rational(1, 2)*sp.LeviCivita(a, b, c, d)*F[c, d] for c in range(4) for d in range(4))
    return out
S = []
for (a, b), (c, d) in [((0, 1), (2, 3)), ((0, 2), (3, 1)), ((0, 3), (1, 2))]:
    F = sp.zeros(4, 4); F[a, b] = 1; F[b, a] = -1; F[c, d] += 1; F[d, c] -= 1
    if star_form(F) != F:
        F2 = sp.zeros(4, 4); F2[a, b] = 1; F2[b, a] = -1; F2[c, d] -= 1; F2[d, c] += 1
        F = F2
    assert star_form(F) == F
    S.append(F)
# symmetric traceless phi: basis of 5
ph = sp.symbols('f0:5')
phi = sp.Matrix([[ph[0], ph[1], ph[2]], [ph[1], ph[3], ph[4]], [ph[2], ph[4], -ph[0] - ph[3]]])


def M_of(Kv):
    dvec = [sp.I*Kv[0], sp.I*Kv[1], sp.I*Kv[2], lam]
    eqs = []
    for nu, rho, sg in itertools.product(range(4), repeat=3):
        e = 0
        for mu in range(4):
            e += dvec[mu]*sum(phi[a, b]*S[a][mu, nu]*S[b][rho, sg] for a in range(3) for b in range(3))
        e = sp.expand(e)
        if e != 0:
            eqs.append(e)
    return sp.Matrix([[sp.diff(e, f) for f in ph] for e in eqs])


def gcd_minors(M, Kn):
    """gcd of all 5x5 minors.  Every minor vanishes at lam = +-|K| iff rank M(+-|K|) < 5 (then (lam^2 - K^2) divides each
    minor, the roots being simple); so if two minors already have gcd c (lam^2 - K^2) and rank M(+-|K|) = 4, the gcd of ALL
    minors is (lam^2 - K^2) (Smith form: solution space of M(d/dz) phi = 0 has dimension deg gcd = 2)."""
    import random
    rows = list(range(M.shape[0]))
    combos = list(itertools.combinations(rows, 5))
    random.Random(1).shuffle(combos)
    g = None
    for idx in combos:
        dmin = sp.expand(M.extract(list(idx), list(range(5))).det())
        if dmin == 0:
            continue
        g = dmin if g is None else sp.gcd(g, dmin)
        if sp.degree(g, lam) <= 2:
            break
    r_p = M.subs(lam, Kn).rank(simplify=True); r_m = M.subs(lam, -Kn).rank(simplify=True)
    return sp.factor(g), r_p, r_m


for Kv in ([4, 0, 3], [3, 1, 4]):
    Kn = sp.sqrt(sum(sp.Integer(v)**2 for v in Kv))
    M = M_of(Kv)
    g, r_p, r_m = gcd_minors(M, Kn)
    g_monic = sp.factor(sp.Poly(g, lam).monic().as_expr())
    good = sp.simplify(g_monic - (lam**2 - Kn**2)) == 0 and r_p == 4 and r_m == 4
    log(f'(i) K = {Kv}: {M.shape[0]} equations; rank M(lam) generic = {M.rank()}; rank M(+|K|) = {r_p}, M(-|K|) = {r_m}; '
        f'gcd of 5x5 minors (monic) = {g_monic}; == lam^2 - |K|^2: {good}')
    ok &= good and M.rank() == 5
    # (ii) W+ of the regular (+) and irregular (-) modes at K, conformal frame, first order (order (1,1) slot, no legs)
    A = Alg([0, 0, 0], Kv)            # put the mode in the (0,1) slot: momentum K
    for lab, reg, sgn in ((+1, True, -1), (-1, False, +1)):
        gam = mode(Kv, lab, regular=reg)
        H = [[{(0, 1): gam[i, j]} for j in range(N)] for i in range(N)]
        g_ = [[A.add({(0, 0): sp.Integer(1 if i == j else 0)}, H[i][j]) for j in range(N)] for i in range(N)]
        gi = inverse_conformal(A, H)
        G = christoffel(A, g_, gi)
        Rup = riemann_up(A, G)
        Ric = ricci(A, Rup)
        # first-order Weyl, lowered with delta (W^a_bcd is conformally invariant); Schouten with flat background
        Rsc = sum(Ric[i][i].get((0, 1), 0) for i in range(N))
        Sch = [[(Ric[i][j].get((0, 1), 0) - Rsc*(1 if i == j else 0)/6)/2 for j in range(N)] for i in range(N)]
        W = {}
        for a, b, c, d in itertools.product(range(N), repeat=4):
            dl = lambda i, j: 1 if i == j else 0
            W[(a, b, c, d)] = (Rup[(a, b, c, d)].get((0, 1), 0) - dl(a, c)*Sch[b][d] + dl(a, d)*Sch[b][c]
                               - dl(b, d)*Sch[a][c] + dl(b, c)*Sch[a][d])
        # SD part (background star; first order) and its components phi_ab = (1/4) S^a.W.S^b
        Wp = {}
        for a, b, c, d in itertools.product(range(N), repeat=4):
            st = sum(sp.Rational(1, 2)*sp.LeviCivita(a, b, e, f)*W[(e, f, c, d)] for e in range(4) for f in range(4))
            Wp[(a, b, c, d)] = (W[(a, b, c, d)] + st)/2
        phis = sp.zeros(3, 3)
        for A1, B1 in itertools.product(range(3), repeat=2):
            phis[A1, B1] = sp.simplify(sum(S[A1][a, b]*Wp[(a, b, c, d)]*S[B1][c, d] for a in range(4) for b in range(4)
                                           for c in range(4) for d in range(4))/16)
        # check structure z^-1 phis = exp(sgn*|K| z) * constant
        const = sp.simplify(phis/z*sp.exp(-sgn*Kn*z))
        isconst = all(sp.diff(const[i, j], z) == 0 for i in range(3) for j in range(3))
        nonzero = any(const[i, j] != 0 for i in range(3) for j in range(3))
        vec = sp.Matrix([const[0, 0], const[0, 1], const[0, 2], const[1, 1], const[1, 2]])
        inker = sp.simplify(M.subs(lam, sgn*Kn)*vec) == sp.zeros(M.shape[0], 1)
        tr = sp.simplify(const.trace())
        # reconstruct check: Wp == sum phis_ab S^a S^b (Wp is purely SD-SD)
        recon = all(sp.simplify(Wp[(a, b, c, d)] - sum(phis[A1, B1]*S[A1][a, b]*S[B1][c, d] for A1 in range(3)
                                                       for B1 in range(3))/1) == 0 for a, b, c, d in itertools.product(range(4), repeat=4))
        lab_s = '+' if lab > 0 else '-'
        log(f'(ii) K = {Kv}: {"regular" if reg else "irregular"} label-({lab_s}) mode: z^-1 W+ = exp({sgn}|K| z) x const: {isconst}; '
            f'nonzero: {nonzero}; traceless: {tr == 0}; W+ = sum phi S S: {recon}; const in ker M({sgn}|K|): {inker}')
        ok &= isconst and nonzero and inker and tr == 0 and recon

# (iii) b_+ of the explicit ASD solution(s)
for tag in ('0_0_3_4_0_0', '1_2_2_2_-1_2'):
    fn = os.path.join(HERE, 'cache', f't1a_solve2_{tag}.pkl')
    if not os.path.exists(fn):
        log(f'(iii) {tag}: no t1a_solve2 output'); continue
    GamP = eval(pickle.load(open(fn, 'rb'))['GamP'])
    kk = [int(v) for v in tag.split('_')]
    Kv = [kk[0] + kk[3], kk[1] + kk[4], kk[2] + kk[5]]
    Kn, kh, n1, n2 = helicity_basis(Kv)
    G0 = sp.Matrix(3, 3, lambda i, j: sp.simplify(GamP[i, j].subs(z, 0)))
    out = {}
    for lab in (+1, -1):
        e = n1 + lab*sp.I*n2
        ec = n1 - lab*sp.I*n2               # conjugate: e.ec = 2
        out[lab] = sp.nsimplify(sp.simplify((ec.T*G0*ec)[0, 0]/4))   # coefficient of e e in G0 (e e . ec ec = 4)
    log(f'(iii) {tag}: K = {Kv}; TT projections of Gam_p(0): b_+ = {out[1]}, b_- = {out[-1]}  '
        f'(Phi_34 = -b_+ W+[reg + at K] != 0: {out[1] != 0})')
log(f'ALL CHECKS PASS: {ok}')
open(os.path.join(HERE, 'cache', 't1a_lemmaM.txt'), 'w').write('\n'.join(lines) + '\n')
