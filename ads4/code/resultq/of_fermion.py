"""
of_fermion (RESULT_Q_FERMION_HEFT.md, task 1): exact connected <TTTT> of the free massless Dirac fermion in d = 3 (and the generalised
free fermion [psibar psi]_{0,2} family) in the CLS-D complete gauge of v4c_core / og_free, and the SWAPPED (double-commutator) form of
Psi+_k at fixed celestial angles as a double contour integral (same contour machinery as og_free; og_free is imported read-only).

Gamma matrices (mostly-plus metric eta = diag(-1, 1, 1), {g^mu, g^nu} = 2 eta^{mu nu}, real Majorana rep):
    g0 = i sigma2, g1 = sigma1, g2 = sigma3;  slash(x) = g^mu eta_{mu nu} x^nu,  slash(x)^2 = x^2.
Propagator (analytic function; the Wightman functions are its boundary values):  S(x) = <psi(x) psibar(0)> = slash(x) (x^2)^{-p},
    p = Delta_psi + 1/2 (= 3/2 for the free fermion; the overall constant is irrelevant: all results are normalised by <TT>).
Vertex (free Dirac T, null z; no improvement term needed):  T(x, z) = psibar zslash (z.d psi) - (z.d psibar) zslash psi.
    GFF: O(x, z) = psibar zslash (z.d psi) + a (z.d psibar) zslash psi, a fixed by conformality of <OO> (a = -1 found for every Delta_psi).
Wick: a closed fermion loop i1 -> i2 -> ... -> i1 (psi at i contracted with psibar at sigma(i)) contributes
    - tr[Gamma_{i1} S_{i1 i2} Gamma_{i2} S_{i2 i3} ...];  connected <TTTT> = sum over the six oriented 4-cycles (classes 1234, 1243, 1324,
    two orientations each); disconnected (MFT) = products of two 2-loops (s: (12)(34), t: (14)(23), u: (13)(24)).

T-shift branch rule for spinor links (derived in RESULT_Q_FERMION_HEFT.md §2.2 from CLS-D's embedding definition O(x^+) = O_emb(-X) for
the TENSOR T, via the 2-pt and 3-pt functions; NOT tuned):
    same-type pairs (12) (both ^+) and (34) (both unshifted): S = slash(x_ij) (x_ij^2)^{-p} with the Wightman i eps (contour C);
    mixed pairs (13), (14), (23), (24):  S = slash(x_ij) (-x_ij^2)^{-p}, base > 0 on the real domain, coefficient +1;
derivatives w.r.t. patch coordinates.  This is the scalar rule of og_free with the spinor numerator slash(x_ij) unchanged.
Validated (of_run.py sec 'val'): <TT> primary; <T3 T1^+> = - flat formula (CLS-D: H^2/|X13|^5); <T1^+ T2^+ T3> and <T1^+ T3 T4> equal the
CLS-D continuation N_flat / prod |X_ij|^{5/2} of the conformal 3-pt function; MFT s, t, u terms equal v4c_core's exact edge terms.
"""
import itertools
import numpy as np
import og_free as F

dot, zvec, E0 = F.dot, F.zvec, F.E0

G0 = np.array([[0, 1], [-1, 0]], complex)
G1 = np.array([[0, 1], [1, 0]], complex)
G2 = np.array([[1, 0], [0, -1]], complex)
GAM = np.stack([G0, G1, G2])                     # gamma^mu (upper index)


def slash(v):
    """slash(v) = gamma^mu eta_{mu nu} v^nu for v[..., 3] -> [..., 2, 2]."""
    v = np.asarray(v)
    return (-v[..., 0, None, None]*G0 + v[..., 1, None, None]*G1 + v[..., 2, None, None]*G2)


def mm(A, B):
    """batched 2x2 matrix product (explicit, fast)."""
    C = np.empty(np.broadcast_shapes(A.shape, B.shape), complex)
    C[..., 0, 0] = A[..., 0, 0]*B[..., 0, 0] + A[..., 0, 1]*B[..., 1, 0]
    C[..., 0, 1] = A[..., 0, 0]*B[..., 0, 1] + A[..., 0, 1]*B[..., 1, 1]
    C[..., 1, 0] = A[..., 1, 0]*B[..., 0, 0] + A[..., 1, 1]*B[..., 1, 0]
    C[..., 1, 1] = A[..., 1, 0]*B[..., 0, 1] + A[..., 1, 1]*B[..., 1, 1]
    return C


def tr(A):
    return A[..., 0, 0] + A[..., 1, 1]


# ------------------------------------------------------------------ one link and its derivatives
def link_mats(xs, w, r, A, B, ab, zi, zj, sg, p):
    """d_i^m d_j^n [slash(x) w^{-p}], w = sg x^2 (x = x_i - x_j), d_i = z_i.d_{x_i}, d_j = z_j.d_{x_j}; m, n in {0, 1}.
    xs = slash(x) [..., 2, 2]; r = w^{-p} (branch fixed by the caller); A = z_i.x, B = z_j.x, ab = z_i.z_j; zi, zj = slash(z_i), slash(z_j).
    f00 = xs r;  f10 = zi r - 2 p sg A xs r/w;  f01 = -zj r + 2 p sg B xs r/w;
    f11 = 2 p sg (B zi + ab xs + A zj) r/w - 4 p (p+1) A B xs r/w^2."""
    r = r[..., None, None]; w = w[..., None, None]
    A = A[..., None, None]; B = B[..., None, None]
    ab = np.asarray(ab)[..., None, None] if np.ndim(ab) else ab
    rw = r/w
    return {(0, 0): xs*r,
            (1, 0): zi*r - 2*p*sg*A*xs*rw,
            (0, 1): -zj*r + 2*p*sg*B*xs*rw,
            (1, 1): 2*p*sg*(B*zi + ab*xs + A*zj)*rw - 4*p*(p + 1)*A*B*xs*rw/w}


def links_generic(X, P, pairs, sigma, p, wover=None):
    """Links for arbitrary points X[i] (..., 3) and null polarisations P[i] (3,), pairs (i<j) -> {(m, n): mats}.
    sigma[(i,j)] = +1 (base x^2, principal branch of the complex power) or -1 (base -x^2).  wover[(i,j)] = (w, r, A, B, ab, xs) override."""
    out = {}
    for (i, j) in pairs:
        if wover and (i, j) in wover:
            w, r, A, B, ab, xs = wover[(i, j)]
            sg = sigma[(i, j)]
        else:
            x = X[i] - X[j]
            sg = sigma[(i, j)]
            w = sg*dot(x, x) + 0j
            r = w**-p
            A = dot(P[i], x); B = dot(P[j], x); ab = float(dot(P[i], P[j])); xs = slash(x)
        out[(i, j)] = link_mats(xs, w, r, A, B, ab, slash(P[i]), slash(P[j]), sg, p)
    return out


def get(L, i, j, m, n):
    """<psi_i (d_i^m) psibar_j (d_j^n)> from the stored unordered pair: S(x_j - x_i) = -S(x_i - x_j)."""
    if i < j:
        return L[(i, j)][(m, n)]
    return -L[(j, i)][(n, m)]


# vertex terms: (coef, derivative on psi, derivative on psibar)
VERT_T = [(1.0, 1, 0), (-1.0, 0, 1)]


def loop(L, cyc, P, vert):
    """sum over vertex terms of  - tr[ Z_{c0} S_{c0 c1} Z_{c1} S_{c1 c2} ... ]  for the oriented cycle cyc = (c0, c1, ..., c_{n-1})."""
    Zs = {i: slash(P[i]) for i in cyc}
    nloop = len(cyc)
    tot = 0
    for terms in itertools.product(range(len(vert)), repeat=nloop):
        coef = 1.0
        dpsi, dbar = {}, {}
        for v, ti in zip(cyc, terms):
            c, mp, mb = vert[ti]
            coef *= c; dpsi[v] = mp; dbar[v] = mb
        M = None
        for a in range(nloop):
            i, j = cyc[a], cyc[(a + 1) % nloop]
            S = get(L, i, j, dpsi[i], dbar[j])
            Mi = mm(Zs[i], S) if M is None else mm(M, mm(Zs[i], S))
            M = Mi
        tot = tot + coef*(-tr(M))
    return tot


CYCLES = {'1234': [(1, 2, 3, 4), (1, 4, 3, 2)], '1243': [(1, 2, 4, 3), (1, 3, 4, 2)], '1324': [(1, 3, 2, 4), (1, 4, 2, 3)]}
MFT = {'s': ((1, 2), (3, 4)), 't': ((1, 4), (2, 3)), 'u': ((1, 3), (2, 4))}
PAIRS = [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
MIXED = {(1, 3), (1, 4), (2, 3), (2, 4)}


def correlator(L, P, classes=('1234', '1243', '1324'), vert=VERT_T):
    out = {}
    for cl in classes:
        if cl in CYCLES:
            out[cl] = sum(loop(L, c, P, vert) for c in CYCLES[cl])
        else:
            (a, b), (c, d) = MFT[cl]
            out[cl] = loop(L, (a, b), P, vert)*loop(L, (c, d), P, vert)
    return out


# ------------------------------------------------------------------ gauge links (branch rules of og_free / §2.2)
def links_gauge(t1, t3, z1, z3, pol2, pol4, y13, p=1.5, shift=True):
    """All six links in the v4c gauge, broadcast over complex t1, t3.  shift=False: control 'no T-shift' (operator 1 at the patch point x1,
    timelike to 3): every mixed link is the naive (x^2)^{-p} with the Wightman branch -> constant phase relative to shift=True."""
    t1 = np.asarray(t1, complex); t3 = np.asarray(t3, complex)
    sh = np.broadcast(t1, t3).shape
    t1 = np.broadcast_to(t1, sh); t3 = np.broadcast_to(t3, sh)
    X = {1: -t1[..., None]*z1, 2: np.broadcast_to(-E0, sh + (3,)).astype(complex),
         3: t3[..., None]*z3, 4: np.broadcast_to(E0, sh + (3,)).astype(complex)}
    P = {1: z1, 2: pol2, 3: z3, 4: pol4}
    sigma = {pr: (-1 if pr in MIXED else 1) for pr in PAIRS}
    wover = {}
    # (12), (34): exact bases, principal branch on C
    for (i, j), tt in (((1, 2), t1), ((3, 4), t3)):
        x = X[i] - X[j]
        w = 2*tt - 1 + 0j
        wover[(i, j)] = (w, w**-p, dot(P[i], x), dot(P[j], x), float(dot(P[i], P[j])), slash(x))
    # (13): exact invariants; slash(x13) = -t1 z1slash - t3 z3slash
    w = 2*y13*t1*t3
    r = (2*y13)**-p*t1**-p*t3**-p
    xs = -t1[..., None, None]*slash(z1) - t3[..., None, None]*slash(z3)
    wover[(1, 3)] = (w, r, t3*y13, t1*y13, -y13, xs)
    L = links_generic(X, P, PAIRS, sigma, p, wover)
    return L, P


def Eplus(k, phi1, phi3, pol2, pol4, n=80, rad=0.5, classes=('1234', '1243', '1324'), p=1.5, vert=VERT_T, q=0, graded=None):
    """swapped-form E+(phi1, phi3) = oint oint t1^w t3^w f_k G per class, w = 2p + q (= Delta + J - 2 = 3 for T at q = 0)."""
    z1 = zvec(phi1); z3 = zvec(phi3)
    y13 = 2*np.sin((phi1 - phi3)/2)**2
    t, wt = F.contour(n, rad) if graded is None else F.contour_graded(*graded, rad=rad)
    T1 = t[:, None]; T3 = t[None, :]
    L, P = links_gauge(T1, T3, z1, z3, pol2, pol4, y13, p=p)
    G = correlator(L, P, classes, vert)
    T1b, T3b = np.broadcast_arrays(T1, T3)
    wgt = 2*p + q
    h0 = T1b**wgt*T3b**wgt*F.fk(k, T1b, T3b, z1, z3, y13=y13)
    W = wt[:, None]*wt[None, :]
    return {cl: complex(np.sum(W*h0*G[cl])) for cl in classes}


# ------------------------------------------------------------------ two-point normalisation
def twopt(x, z, w, p=1.5, vert=VERT_T, sg=1):
    """<O(x, z) O(0, w)> from one 2-loop; sg = -1: base -x^2 (mixed-link rule)."""
    X = {1: np.asarray(x, complex), 2: np.zeros(3, complex)}
    P = {1: z, 2: w}
    L = links_generic(X, P, [(1, 2)], {(1, 2): sg}, p)
    return loop(L, (1, 2), P, vert)


def flat_TT(x, z, w, Delta=3.0):
    """conformal structure ((I z).w)^2/(x^2)^Delta (v4c C_T = 1 normalisation for Delta = 3)."""
    X2 = dot(x, x)
    Iz = z - 2*x*dot(x, z)/X2
    return dot(Iz, w)**2/X2**Delta
