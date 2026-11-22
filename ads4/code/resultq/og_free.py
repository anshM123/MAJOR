"""
og_free (RESULT_Q_OG_CHECK.md): exact free-field <TTTT> in d = 3 (conformally coupled free scalar, connected = O(1/N) part) in the CLS-D
complete gauge of v4c_core, and the SWAPPED (double-commutator) form of Psi+_{k} at fixed celestial angles as a double contour integral.

Gauge (v4c_core): x0 = 0, x2 = -e (patch point of 2^+), x4 = +e, x3 = t3 z3, x1 = -t1 z1 (patch point of 1^+), z(phi) = (1, cos, sin),
light-ray polarisations z1, z3 at 1, 3; measure t1^3 t3^3; f_k(u', v') with u' = s12 s34/(s13 s24), v' = s14 s23/(s13 s24) (patch s_ij).
Correlator: <T4(x4) T3(x3) T1(x1^+) T2(x2^+)>.

Branch / ^+ conventions (derived in RESULT_Q_OG_CHECK.md §2.2).  The only pairs that can be timelike are (12) (both ^+) and (34) (both
unshifted); (13), (14), (23), (24) link a ^+ point to an unshifted one and are TRULY spacelike on the whole domain.  With the (position-
independent) phase convention phi_flat(x^+) := |Omega|^{1/2} phi_cyl, every free-field link is
    same type (12), (34):   g = (x_ij^2)^{-1/2}      (Wightman i eps: "1 left of 2" <-> t1 + i0, "4 left of 3" <-> t3 + i0)
    mixed (13),(14),(23),(24): g = (-x_ij^2)^{-1/2} > 0  on the real domain,
derivatives w.r.t. patch coordinates.  (For MFT this reproduces "naive patch formulas, sign +1 per Wick term" of v4c_core.)

Swapped form.  For t1, t3 > 1/2 all four orderings of [T4,T3][T1,T2] coincide; the four orderings are the boundary values t1 +- i0,
t3 +- i0 of ONE analytic function h(t1, t3) = t1^3 t3^3 f_k G (Lemma C: f's poles carry the same i eps).  Hence
    E+(a,b) := int dt1 dt3 t1^3 t3^3 f_k <[T4,T3][T1,T2]>  =  oint_C dt1 oint_C dt3 h(t1, t3),
C = the circle |t - 1/2| = 1/2 traversed counter-clockwise from t = 0 (i.e. -(hairpin) per variable; for rational h this is
(-2 pi i)^2 Res_{1/2} Res_{1/2}, the v4c_core residue formula).  t = 0 (x1^+ -> x0) is the only endpoint; the integrand is O(t^{-1/2})
there and the substitution theta = pi (1 - cos psi) makes it smooth.  Branches on C: (2t-1)^{-1/2} principal (cut t < 1/2 = the
double-commutator support), (13): (2(1-c))^{-1/2} t1^{-1/2} t3^{-1/2} principal, the rest have Re(base) > 0.

Free scalar, conformal coupling xi = 1/8:  T(x, z) = (3/4)(z.d phi)^2 - (1/4) phi (z.d)^2 phi   (z null).
Connected 4-pt function = sum over perfect matchings of the 8 legs that form one 4-cycle (three classes: 1234, 1243, 1324).
Disconnected (MFT) part = matchings forming two 2-cycles (s: (12)(34), t: (14)(23), u: (13)(24)).
"""
import itertools
import numpy as np
import sympy as sp

E0 = np.array([1.0, 0.0, 0.0])
ETA = np.diag([-1.0, 1.0, 1.0])


def dot(a, b):
    return -a[..., 0]*b[..., 0] + a[..., 1]*b[..., 1] + a[..., 2]*b[..., 2]


def zvec(phi):
    phi = np.asarray(phi, float)
    return np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], axis=-1)


# ------------------------------------------------------------------ directional derivatives of (sigma X^2)^{-1/2}
def _derive(MMAX=2):
    """D(m, n) = (a.d_X)^m (-b.d_X)^n (sigma X.X)^{-1/2} for null a, b, as a function of Y = sigma X.X, r = Y^{-1/2}, A = a.X, B = b.X,
    ab = a.b, sigma.  Returns {(m, n): lambdified f(r, Y, A, B, ab, sg)}."""
    r, Y, A, B, ab, sg = sp.symbols('r Y A B ab sg')

    def dY(e):
        return sp.diff(e, Y) + sp.diff(e, r)*(-r/(2*Y))

    def da(e):
        return dY(e)*2*sg*A + sp.diff(e, B)*ab        # a.dA = a.a = 0

    def db(e):
        return dY(e)*2*sg*B + sp.diff(e, A)*ab
    out = {}
    for n in range(MMAX + 1):
        e = r
        for _ in range(n):
            e = -db(e)
        for m in range(MMAX + 1):
            em = e
            for _ in range(m):
                em = da(em)
            out[(m, n)] = sp.lambdify((r, Y, A, B, ab, sg), sp.simplify(em), 'numpy')
    return out


DFUN = _derive()

# vertex: list of (coef, (mA, mB)) -- number of z-derivatives on each of the two legs
VERTEX = [(0.75, (1, 1)), (-0.25, (0, 2))]

PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
MIXED = {(1, 3), (1, 4), (2, 3), (2, 4)}


def _matchings(legs):
    if not legs:
        yield []
        return
    a = legs[0]
    for i in range(1, len(legs)):
        b = legs[i]
        rest = legs[1:i] + legs[i + 1:]
        for m in _matchings(rest):
            yield [(a, b)] + m


def _classify(match):
    """match: list of ((v,s),(w,s')) ; returns 'X' (invalid), cycle class '1234' / '1243' / '1324', or 's'/'t'/'u' (two 2-cycles)."""
    edges = []
    for (v, _), (w, _) in match:
        if v == w:
            return 'X'
        edges.append(tuple(sorted((v, w))))
    es = sorted(edges)
    if len(set(es)) == 2:                      # two double edges
        p = set(es)
        if (1, 2) in p:
            return 's'
        if (1, 4) in p:
            return 't'
        return 'u'
    nb1 = sorted([w for (v, w) in edges if v == 1] + [v for (v, w) in edges if w == 1])
    return {(2, 4): '1234', (2, 3): '1243', (3, 4): '1324'}[tuple(nb1)]


LEGS = [(v, s) for v in (1, 2, 3, 4) for s in (0, 1)]
MATCHES = {}
for _m in _matchings(LEGS):
    c = _classify(_m)
    if c != 'X':
        MATCHES.setdefault(c, []).append(_m)


# ------------------------------------------------------------------ links
def links(t1, t3, z1, z3, pol2, pol4, y13=None, tube=0.0):
    """Link values D[(i,j)][(m,n)] for all pairs, broadcast over t1, t3 (complex arrays of a common shape).
    z1, z3: (3,) null; pol2, pol4: (3,) null polarisations of T2, T4.  y13 = 1 - cos(phi1 - phi3) = 2 sin^2(dphi/2) (exact (13) link).
    tube: optional timelike-tube shift x1 -> x1 + i*tube*e (control only; then the (13) base is computed directly, principal branch)."""
    t1 = np.asarray(t1, complex); t3 = np.asarray(t3, complex)
    sh = np.broadcast(t1, t3).shape
    t1 = np.broadcast_to(t1, sh); t3 = np.broadcast_to(t3, sh)
    X = {1: -t1[..., None]*z1 + 1j*tube*E0, 2: np.broadcast_to(-E0, sh + (3,)).astype(complex),
         3: t3[..., None]*z3, 4: np.broadcast_to(E0, sh + (3,)).astype(complex)}
    P = {1: z1, 2: pol2, 3: z3, 4: pol4}
    out = {}
    for (i, j) in PAIRS:
        Xij = X[i] - X[j]
        s = dot(Xij, Xij)
        sg = -1.0 if (i, j) in MIXED else 1.0
        Y = sg*s
        if (i, j) == (1, 2):
            Y = Y if tube else 2*t1 - 1 + 0j
            r = Y**-0.5
        elif (i, j) == (3, 4):
            Y = 2*t3 - 1 + 0j
            r = Y**-0.5
        elif (i, j) == (1, 3) and not tube:
            # exact invariants (no cancellation at small dphi): y = 1 - c13 = -z1.z3 = 2 sin^2(dphi/2)
            y = y13
            Y = 2*y*t1*t3
            r = (2*y)**-0.5*t1**-0.5*t3**-0.5
            A = t3*y; B = t1*y; ab = -y
            out[(i, j)] = {mn: DFUN[mn](r, Y, A, B, ab, sg) for mn in DFUN}
            continue
        else:
            r = Y**-0.5
        A = dot(P[i], Xij); B = dot(P[j], Xij); ab = float(dot(P[i], P[j]))
        out[(i, j)] = {mn: DFUN[mn](r, Y, A, B, ab, sg) for mn in DFUN}
    return out


def _link(L, a, b):
    """value of the contraction of leg a = (v, m) with leg b = (w, n)."""
    (v, m), (w, n) = a, b
    if v < w:
        return L[(v, w)][(m, n)]
    return L[(w, v)][(n, m)]


def correlator(L, classes=('1234', '1243', '1324')):
    """sum over vertex terms and matchings; returns {class: array}."""
    out = {}
    for cl in classes:
        tot = 0
        for terms in itertools.product(range(len(VERTEX)), repeat=4):
            coef = 1.0
            der = {}
            for v, ti in zip((1, 2, 3, 4), terms):
                c, (mA, mB) = VERTEX[ti]
                coef *= c
                der[(v, 0)] = mA; der[(v, 1)] = mB
            acc = 0
            for match in MATCHES[cl]:
                p = 1
                for a, b in match:
                    p = p*_link(L, (a[0], der[a]), (b[0], der[b]))
                acc = acc + p
            tot = tot + coef*acc
        out[cl] = tot
    return out


# ------------------------------------------------------------------ subtraction factors and contour
def fk(k, t1, t3, z1, z3, y13=None):
    x1 = -t1[..., None]*z1; x3 = t3[..., None]*z3
    x2 = -E0; x4 = E0
    s12, s34, s13, s14, s23, s24 = (dot(a, a) for a in (x1 - x2, x3 - x4, x1 - x3, x1 - x4, x2 - x3, x2 - x4 + 0*x1))
    if y13 is not None:
        s13 = -2*t1*t3*y13                      # exact
    up = s12*s34/(s13*s24); vp = s14*s23/(s13*s24)
    if k == 2:
        return vp - up
    if k == 3:
        return -np.ones_like(up)
    if k == 4:
        return ((vp - up)**3 - (vp - up))/(4*up*vp)
    if k == 5:
        return -((up - vp)**2 - 1)/(4*up*vp)
    if k == 0:
        return np.ones_like(up)
    raise ValueError(k)


def contour(n, rad=0.5):
    """nodes t, weights w (dt included) on C: t = 1/2 - rad e^{i theta} shifted so that t(0) = 0, ccw; theta = pi(1 - cos psi).
    rad: radius; the circle is centred at rad (passes through 0 and 2 rad)."""
    x, wx = np.polynomial.legendre.leggauss(n)
    psi = (x + 1)*np.pi/2; wpsi = wx*np.pi/2
    th = np.pi*(1 - np.cos(psi)); dth = np.pi*np.sin(psi)
    t = rad - rad*np.exp(1j*th)
    dt = -1j*rad*np.exp(1j*th)*dth
    return t, wpsi*dt


def contour_graded(s0, m=16, q=4.0, rad=0.5):
    """Composite Gauss on C with geometric panels towards both ends theta = 0, 2 pi (the endpoint t = 0, where the integrand has
    structure on the scale |t| ~ dphi^2): panels [0, s0] (theta = s0 u^2), [s0, s0 q], ..., up to pi; mirrored on (pi, 2 pi)."""
    x, wx = np.polynomial.legendre.leggauss(m)
    th, wt = [], []
    u = (x + 1)/2; wu = wx/2                         # first panel, theta = s0 u^2
    th.append(s0*u**2); wt.append(wu*2*s0*u)
    a = s0
    while a < np.pi:
        b = min(a*q, np.pi)
        th.append(a + (b - a)*(x + 1)/2); wt.append(wx*(b - a)/2)
        a = b
    th = np.concatenate(th); wt = np.concatenate(wt)
    th = np.concatenate([th, 2*np.pi - th[::-1]]); wt = np.concatenate([wt, wt[::-1]])
    t = rad - rad*np.exp(1j*th)
    dt = -1j*rad*np.exp(1j*th)
    return t, wt*dt


def Eplus(k, phi1, phi3, pol2, pol4, n=80, rad=0.5, classes=('1234', '1243', '1324'), tube=0.0, graded=None):
    """swapped-form E+(phi1, phi3) = oint oint t1^3 t3^3 f_k G, per class.  Returns dict class -> complex.
    graded: None (cosine-substituted Gauss with n nodes) or (s0, m, q) for contour_graded."""
    z1 = zvec(phi1); z3 = zvec(phi3)
    y13 = 2*np.sin((phi1 - phi3)/2)**2
    t, w = contour(n, rad) if graded is None else contour_graded(*graded, rad=rad)
    T1 = t[:, None]; T3 = t[None, :]
    L = links(T1, T3, z1, z3, pol2, pol4, y13=y13, tube=tube)
    G = correlator(L, classes)
    T1b, T3b = np.broadcast_arrays(T1, T3)
    h0 = T1b**3*T3b**3*fk(k, T1b, T3b, z1, z3, y13=y13)
    W = w[:, None]*w[None, :]
    return {cl: complex(np.sum(W*h0*G[cl])) for cl in classes}


# ------------------------------------------------------------------ EE contraction of the (2,4) structure
def null_basis():
    """5 null vectors whose n n^T span the symmetric traceless 3x3 tensors; Gram (n_a.n_b)^2 and its inverse."""
    ang = [0.3, 1.4, 2.5, 3.9, 5.2]
    ns = []
    for i, a in enumerate(ang):
        # null vectors with varying spatial direction AND a boost so that the 5 rank-1 tensors are independent
        v = np.array([1.0, np.cos(a), np.sin(a)])
        ns.append(v)
    ns = np.array(ns)
    Gm = np.array([[dot(a, b)**2 for b in ns] for a in ns])
    return ns, Gm, np.linalg.inv(Gm)


def EE_combine(fun):
    """fun(pol2, pol4) -> dict; returns sum_ab (G^-1)_ab fun(n_b at 2, n_a at 4) = EE contraction T4^{mu nu} T2_{mu nu}."""
    ns, Gm, Gi = null_basis()
    acc = None
    for a in range(5):
        for b in range(5):
            r = fun(ns[b], ns[a])
            if acc is None:
                acc = {k: Gi[a, b]*v for k, v in r.items()}
            else:
                for k in acc:
                    acc[k] = acc[k] + Gi[a, b]*r[k]
    return acc
