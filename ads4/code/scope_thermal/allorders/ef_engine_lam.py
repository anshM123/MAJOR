"""
ef_engine.py on the O(lam) J-corrected brane (N = 1 - 6 lam/r^6), everything truncated at O(lam).

Hall viscosity of the AdS4 black brane from parity-odd curvature terms, at first order in the coupling.  [computed]

Independent of any bivector shortcut: all Weyl contractions are explicit index sums over sparse component
dictionaries; every quantity is carried as an exact eps-truncated triple (c0, c1, c2) (eps = perturbation amplitude).

Setup (ingoing EF, L_AdS = 1, r_h = 1):
    ds^2 = -f dv^2 + 2 dv dr + r^2 (e^{eps H})_ij dx^i dx^j,  H = [[phi, chi], [chi, -phi]],  f = r^2 - 1/r,
    phi = phi(v, r), chi = chi(v, r);   e^{eps H} = 1 + eps H + eps^2 H^2/2 (exact to eps^2; det = 1 to eps^2).
    Orientation: eps_{vrxy} = + sqrt(-g) = + r^2   (also selectable with ORIENT = -1).

Operators (Weyl tensor C, C~_{abcd} = 1/2 eps_{ab}^{ef} C_{efcd}):
    I = C_abcd C^abcd,  P = C_abcd C~^abcd,  J = C_ab^cd C_cd^ef C_ef^ab,  Q = C_ab^cd C_cd^ef C~_ef^ab,
    and theta(r) P for an arbitrary profile theta (Saremi-Son check).

Extraction.  L2 = eps^2 coefficient of sqrt(-g) * O, a quadratic form in the jets
    (phi, phi_v, phi_r, phi_vv, phi_vr, phi_rr; same for chi).  Cross terms phi-chi come only from P-odd operators.
    The chi equation of motion at O(omega) on the ingoing solution (phi_v = O(w), phi_r = O(w), phi_vr = O(w^2)) is
        E_chi = E_chi^Einstein + lam * sigma(r) * phi_v + O(w^2),
    with sigma read off from the cross terms by the Euler-Lagrange operator (computed symbolically below, no
    hand formula).  The Einstein part is checked to be  N * [ d_r( r^2 (f chi_r + chi_v) ) + d_v( r^2 chi_r ) ],
    i.e. -N d_r Pi_min - ..., with Pi_min = -r^2 (f chi_r + chi_v) the canonical momentum of the shear mode.
    Hence d_r Pi_min = (lam/N) sigma phi_v  =>  Pi_min(oo) = Pi_min(1) + (lam/N) phi_v int_1^oo sigma dr.
    With T^{xy} = Pi_min(oo)/(16 pi G) (calibrated by eta: Pi_min(1) = -chi_v gives eta = 1/16piG, 4 pi eta/s = 1)
    and delta T^{xy} = -p chi - eta chi_v + eta_H phi_v:
        4 pi eta_H / s = (lam/N) int_1^oo sigma(r) dr.
    Odd terms vanish at r -> oo (C -> 0), so they add nothing to the boundary momentum.
"""
import sys
import itertools
import sympy as sp

ORIENT = 1
v, r, x, y = sp.symbols('v r x y')
eps = sp.Symbol('eps')
X = [v, r, x, y]
phi = sp.Function('phi')(v, r)
chi = sp.Function('chi')(v, r)
th = sp.Function('theta')(r)
lamb = sp.Symbol('lam')
f = r ** 2 - 1 / r + lamb * (-8 / r + 24 / r ** 4 - 16 / r ** 7)     # O(lam) J-corrected brane (thermo_check.py)
Nb = 1 - 6 * lamb / r ** 6
Ninv = 1 + 6 * lamb / r ** 6                                          # 1/N to O(lam)
n = 4

# ---------------- truncated-series arithmetic: triples (c0, c1, c2)
Z = (sp.Integer(0),) * 3


def add(a, b):
    return tuple(a[i] + b[i] for i in range(3))


def mul(a, b):
    return (a[0] * b[0], a[0] * b[1] + a[1] * b[0], a[0] * b[2] + a[1] * b[1] + a[2] * b[0])


def scal(c, a):
    return tuple(c * a[i] for i in range(3))


def d(a, var):
    return tuple(sp.diff(a[i], var) for i in range(3))


def _tr(e):
    e = sp.expand(e)
    return e.coeff(lamb, 0) + lamb * e.coeff(lamb, 1)


def simp(a):
    return tuple(_tr(a[i]) for i in range(3))


def is0(a):
    return all(sp.expand(t) == 0 for t in a)


# ---------------- metric and exact-to-eps^2 inverse
H = sp.Matrix([[phi, chi], [chi, -phi]])
H2 = H * H
g = [[Z for _ in range(n)] for _ in range(n)]
gi = [[Z for _ in range(n)] for _ in range(n)]
g[0][0] = (_tr(-Nb ** 2 * f), 0, 0); g[0][1] = g[1][0] = (Nb, 0, 0)
gi[0][1] = gi[1][0] = (Ninv, 0, 0); gi[1][1] = (_tr(f), 0, 0)
for i in range(2):
    for j in range(2):
        g[2 + i][2 + j] = (r ** 2 * (1 if i == j else 0), r ** 2 * H[i, j], r ** 2 * H2[i, j] / 2)
        gi[2 + i][2 + j] = ((1 if i == j else 0) / r ** 2, -H[i, j] / r ** 2, H2[i, j] / (2 * r ** 2))
g = [[tuple(sp.sympify(t) for t in g[a][b]) for b in range(n)] for a in range(n)]
gi = [[tuple(sp.sympify(t) for t in gi[a][b]) for b in range(n)] for a in range(n)]
# check g gi = 1 to eps^2
for a in range(n):
    for b in range(n):
        s = Z
        for c in range(n):
            s = add(s, mul(g[a][c], gi[c][b]))
        assert is0(simp(add(s, ((-1 if a == b else 0), 0, 0)))), (a, b)

# ---------------- Christoffels, Riemann
dg = [[[d(g[a][b], X[c]) for c in range(n)] for b in range(n)] for a in range(n)]
Gam = [[[None] * n for _ in range(n)] for _ in range(n)]
for a in range(n):
    for b in range(n):
        for c in range(b, n):
            s = Z
            for e in range(n):
                if is0(gi[a][e]):
                    continue
                s = add(s, mul(gi[a][e], scal(sp.Rational(1, 2), add(add(dg[e][b][c], dg[e][c][b]), scal(-1, dg[b][c][e])))))
            s = simp(s)
            Gam[a][b][c] = Gam[a][c][b] = s
Rud = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            for e in range(c + 1, n):
                s = add(d(Gam[a][e][b], X[c]), scal(-1, d(Gam[a][c][b], X[e])))
                for h in range(n):
                    s = add(s, add(mul(Gam[a][c][h], Gam[h][e][b]), scal(-1, mul(Gam[a][e][h], Gam[h][c][b]))))
                s = simp(s)
                Rud[a, b, c, e] = s
                Rud[a, b, e, c] = scal(-1, s)
Rl = {}
for a in range(n):
    for b in range(n):
        for c in range(n):
            for e in range(n):
                if c == e:
                    Rl[a, b, c, e] = Z
                    continue
                s = Z
                for h in range(n):
                    if is0(g[a][h]):
                        continue
                    s = add(s, mul(g[a][h], Rud[h, b, c, e]))
                Rl[a, b, c, e] = simp(s)

# ---------------- jets: replace derivatives by symbols (everything below is algebraic)
J = {}
for nm, F in (('p', phi), ('q', chi)):
    J[nm] = sp.Symbol(nm)
    J[nm + 'v'] = sp.Symbol(nm + 'v'); J[nm + 'r'] = sp.Symbol(nm + 'r')
    J[nm + 'vv'] = sp.Symbol(nm + 'vv'); J[nm + 'vr'] = sp.Symbol(nm + 'vr'); J[nm + 'rr'] = sp.Symbol(nm + 'rr')
jet_sub = []
for nm, F in (('p', phi), ('q', chi)):
    jet_sub += [(sp.Derivative(F, v, v), J[nm + 'vv']), (sp.Derivative(F, v, r), J[nm + 'vr']),
                (sp.Derivative(F, r, r), J[nm + 'rr']), (sp.Derivative(F, v), J[nm + 'v']),
                (sp.Derivative(F, r), J[nm + 'r']), (F, J[nm])]


def jet(a):
    out = []
    for t in a:
        t = t.subs(jet_sub[:5] + jet_sub[6:11])   # derivatives first
        t = t.subs({phi: J['p'], chi: J['q']})
        out.append(sp.expand(t))
    return tuple(out)


Rl = {k: jet(vv) for k, vv in Rl.items()}
g = [[jet(t) for t in row] for row in g]
gi = [[jet(t) for t in row] for row in gi]

# ---------------- Ricci, scalar, Weyl (all lower)
Ric = {}
for b in range(n):
    for e in range(n):
        s = Z
        for a in range(n):
            for c in range(n):
                if is0(gi[a][c]):
                    continue
                s = add(s, mul(gi[a][c], Rl[a, b, c, e]))
        Ric[b, e] = simp(s)
Rs = Z
for b in range(n):
    for e in range(n):
        if not is0(gi[b][e]):
            Rs = add(Rs, mul(gi[b][e], Ric[b, e]))
Rs = simp(Rs)
W = {}
for a, b, c, e in itertools.product(range(n), repeat=4):
    s = Rl[a, b, c, e]
    s = add(s, scal(-sp.Rational(1, 2), add(add(mul(g[a][c], Ric[b, e]), scal(-1, mul(g[a][e], Ric[b, c]))),
                                           add(scal(-1, mul(g[b][c], Ric[a, e])), mul(g[b][e], Ric[a, c])))))
    s = add(s, mul(scal(sp.Rational(1, 6), Rs), add(mul(g[a][c], g[b][e]), scal(-1, mul(g[a][e], g[b][c])))))
    W[a, b, c, e] = simp(s)

# mixed W_ab^cd = W_ab ef g^ec g^fd
Wm = {}
for a, b in itertools.product(range(n), repeat=2):
    for c, e in itertools.product(range(n), repeat=2):
        s = Z
        for p_, q_ in itertools.product(range(n), repeat=2):
            if is0(gi[p_][c]) or is0(gi[q_][e]):
                continue
            s = add(s, mul(W[a, b, p_, q_], mul(gi[p_][c], gi[q_][e])))
        Wm[a, b, c, e] = simp(s)
# Levi-Civita eps_{abcd} = ORIENT * sqrt(-g) [abcd], sqrt(-g) = r^2 (exact to eps^2: det e^{eps H} = 1 + O(eps^3))
sqrtg = (_tr(Nb * r ** 2), sp.Integer(0), sp.Integer(0))
# eps_ab^cd = eps_{ab pq} g^pc g^qd
epsm = {}
for a, b, c, e in itertools.product(range(n), repeat=4):
    s = Z
    for p_, q_ in itertools.product(range(n), repeat=2):
        lc = sp.LeviCivita(a, b, p_, q_)
        if lc == 0 or is0(gi[p_][c]) or is0(gi[q_][e]):
            continue
        s = add(s, scal(ORIENT * lc, mul(sqrtg, mul(gi[p_][c], gi[q_][e]))))
    epsm[a, b, c, e] = simp(s)
# dual: Wt_ab^cd = 1/2 eps_ab^ef W_ef^cd
Wtm = {}
for a, b, c, e in itertools.product(range(n), repeat=4):
    s = Z
    for p_, q_ in itertools.product(range(n), repeat=2):
        if is0(epsm[a, b, p_, q_]):
            continue
        s = add(s, scal(sp.Rational(1, 2), mul(epsm[a, b, p_, q_], Wm[p_, q_, c, e])))
    Wtm[a, b, c, e] = simp(s)

nz = lambda D: {k: vv for k, vv in D.items() if not is0(vv)}
Wm_nz, Wtm_nz = nz(Wm), nz(Wtm)


def contract2(A, B):
    """A_ab^cd B_cd^ab (sum over all index values)."""
    s = Z
    for (a, b, c, e), va in A.items():
        vb = B.get((c, e, a, b))
        if vb is not None:
            s = add(s, mul(va, vb))
    return simp(s)


def contract3(A, B, C):
    """A_ab^cd B_cd^ef C_ef^ab."""
    # AB_ab^ef
    AB = {}
    for (a, b, c, e), va in A.items():
        for (c2, e2, p_, q_), vb in B.items():
            if (c2, e2) != (c, e):
                continue
            AB[a, b, p_, q_] = add(AB.get((a, b, p_, q_), Z), mul(va, vb))
    return contract2({k: simp(vv) for k, vv in AB.items()}, C)


I_ = contract2(Wm_nz, Wm_nz)
P_ = contract2(Wm_nz, Wtm_nz)
J_ = contract3(Wm_nz, Wm_nz, Wm_nz)
Q_ = contract3(Wm_nz, Wm_nz, Wtm_nz)
print("background (eps^0):  I =", sp.simplify(I_[0]), " P =", sp.simplify(P_[0]), " J =", sp.simplify(J_[0]),
      " Q =", sp.simplify(Q_[0]))
print("eps^1 pieces (expect 0):", [sp.simplify(t[1]) for t in (I_, P_, J_, Q_)])
jets_p = [J[k] for k in ('p', 'pv', 'pr', 'pvv', 'pvr', 'prr')]
jets_q = [J[k] for k in ('q', 'qv', 'qr', 'qvv', 'qvr', 'qrr')]


def cross_part(expr):
    e = sp.expand(expr)
    out = 0
    for a_ in jets_p:
        for b_ in jets_q:
            c = e.coeff(a_, 1).coeff(b_, 1)
            if c != 0:
                out += sp.factor(c) * a_ * b_
    return out


print("P^(2) cross terms :", cross_part(P_[2]))
print("Q^(2) cross terms :", cross_part(Q_[2]))

# ---------------- Einstein quadratic Lagrangian (for N) : sqrt(-g)(R + 6), eps^2 part
LE2 = sp.expand(r ** 2 * Rs[2])


def EL_from_jets(L2, field):
    """Euler-Lagrange expression E_field[ phi(v,r), chi(v,r) ] (functions restored)."""
    rev = {J['p']: phi, J['pv']: sp.diff(phi, v), J['pr']: sp.diff(phi, r), J['pvv']: sp.diff(phi, v, 2),
           J['pvr']: sp.diff(phi, v, r), J['prr']: sp.diff(phi, r, 2),
           J['q']: chi, J['qv']: sp.diff(chi, v), J['qr']: sp.diff(chi, r), J['qvv']: sp.diff(chi, v, 2),
           J['qvr']: sp.diff(chi, v, r), J['qrr']: sp.diff(chi, r, 2)}
    nm = 'p' if field == 'phi' else 'q'
    out = 0
    for k, ops in ((nm, ()), (nm + 'v', (v,)), (nm + 'r', (r,)), (nm + 'vv', (v, v)), (nm + 'vr', (v, r)), (nm + 'rr', (r, r))):
        dl = sp.diff(L2, J[k]).subs(rev)
        for var in ops:
            dl = -sp.diff(dl, var)
        out += dl
    return sp.expand(out)


Nfac = 1


def hall(L2odd, label, integrate=True):
    Eodd = EL_from_jets(L2odd, 'chi')
    PV = sp.Derivative(phi, v)
    sig, w1, w0, w2 = 0, [], [], []
    for t in sp.Add.make_args(Eodd):
        if not t.has(phi):
            w0.append(t); continue
        ders = [a_ for a_ in t.atoms(sp.Derivative) if a_.expr == phi]
        if ders == [PV] and not (t / PV).has(phi):
            sig += t / PV                       # O(w): phi_v times a function of r
        elif not ders:
            w0.append(t)                        # undifferentiated phi: would be O(1) -- must be absent
        elif all(all(vv == r for vv in a_.variables) for a_ in ders):
            w1.append(t)                        # only r-derivatives: O(w) on the ingoing solution -- must be absent
        else:
            w2.append(t)                        # >= 2 derivatives incl. one v and one more: O(w^2)
    sig = sp.simplify(sig)
    print(f"[{label}] sigma(r) = {sig}")
    print(f"[{label}]   absent-check: undifferentiated-phi terms {w0};  pure-r-derivative terms {w1}")
    print(f"[{label}]   O(w^2) terms dropped: {len(w2)}")
    if not integrate:
        return sig
    integral = sp.integrate(sig, (r, 1, sp.oo))
    out = sp.simplify(integral / Nfac)
    print(f"[{label}]   int_1^oo sigma dr = {sp.simplify(integral)}   ->  4 pi eta_H/s = {out} * coupling")
    return out


