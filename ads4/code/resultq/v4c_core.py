"""
v4c_core (RESULT_Q_V4PRIME.md): exact light-family edge terms of the CLS-D spinning functional Psi+_{k,nu}[f G] on <TTTT> (d = 3),
in the COMPLETE gauge of EE_FOUNDATIONS §1.4 (x0, x2, x4, z0 fixed), where x12^2 and x34^2 are independent variables.

Gauge (final CLS-D configuration 4 > 3 > 0 > 1 > 2, eq. spinningfunctional_definition / scalarfunctional_plus_final):
    x0 = 0,  x2 = -e,  x4 = +e,  e = (1,0,0) (e^2 = -1, mostly-plus metric),  z(phi) = (1, cos phi, sin phi),
    x3 = t3 z3 (future light ray of x0, first half of L[T3]),  x1 = -t1 z1 (patch coordinate of x1^+, second half of L[T1]).
    Light-ray polarisations z1, z3 (CLS-D: O_1(x1^+, x_{1+0}), O_3(x3, x_30)); light-transform weight t^(Delta+J-2) = t^3 for T.
    Exact invariants:  x12^2 = 2t1-1, x34^2 = 2t3-1, x14^2 = -(1+2t1), x23^2 = -(1+2t3), x13^2 = 2 t1 t3 (z1.z3), x24^2 = -4
    u' = (1-2t1)(1-2t3)/(8 t1 t3 (1-c)),  v' = (1+2t1)(1+2t3)/(8 t1 t3 (1-c)),  c = cos(phi1 - phi3).
    The ^+ images change only signs (-1)^(Delta+J) = -1 per 2-pt factor that links a ^+ point to a non-^+ point: every MFT Wick term
    has an even number of them, so all three Wick terms enter with their naive patch formulas (sign +1).
    Edge: the double commutator Delta_12 Delta_34 flips only the (12), (34) i-eps (Lemma C); along each light ray x12^2, x34^2 are
    linear in t, so int dt1 dt3 Delta_12 Delta_34 [g] = (-2 pi i)^2 Res_{t1=1/2} Res_{t3=1/2} g  (exact, rational g).
    The light-cone edge u' = 0 is t1 = t3 = 1/2 for every (phi1, phi3); there v' = 2/(1-c): the CORNER r -> 0, rho -> 1 is phi1 -> phi3.
    (Away from the edge, phi1 -> phi3 at fixed t is the large-eta Regge corner at fixed rho = 2t: ((1-rho)/(1+rho))^2 = u'/v'.)

T two-point function (C_T = 1):  <T_{mu nu}(x) T(y, z)> = v_mu v_nu / (X^2)^3,  v = I(X) z,  I(X) = 1 - 2 X X / X^2,  X = x - y (z null).
Dual (2,4) structure: in this gauge the EE structure is E ~ (I(x02) z2).(I(x04) z4) = z2.z4 (verified numerically below), so
EE <-> full contraction T4^{mu nu} T2_{mu nu}. Since the gauge is complete, the sum rule must hold for EVERY constant (2,4) tensor S
(EE_FOUNDATIONS §1.4): the full 5x5 tensor of the light side is also returned.
Celestial j' = 0 structure (delta1 = delta3 = Delta_T - 1 = 2, delta = 1/2 + i nu):
    C = (-2 z1.z3)^{(5/2 + i nu)/2} / [(-2 z1.z0)(-2 z3.z0)]^{(1/2 + i nu)/2}.
"""
import numpy as np

ETA = np.diag([-1.0, 1.0, 1.0])
E0 = np.array([1.0, 0.0, 0.0])


def dot(a, b):
    """Minkowski dot over the last axis."""
    return -a[..., 0]*b[..., 0] + a[..., 1]*b[..., 1] + a[..., 2]*b[..., 2]


def zvec(phi):
    phi = np.asarray(phi, float)
    return np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], axis=-1)


def inv_vec(X, z):
    """I(X) z = z - 2 X (X.z)/X^2 (complex ok)."""
    X2 = dot(X, X)
    return z - 2.0*X*(dot(X, z)/X2)[..., None]


def lower(v):
    return v*np.array([-1.0, 1.0, 1.0])


# symmetric traceless basis (5 elements, orthonormal under the trace pairing with eta) for the (2,4) tensor components
def st_basis():
    B = []
    # in lower-index components; trace with eta^{mu nu}
    def sym(i, j):
        M = np.zeros((3, 3)); M[i, j] += 1; M[j, i] += 1; return M/2
    cand = [sym(0, 1), sym(0, 2), sym(1, 2), np.diag([0., 1., -1.]), np.diag([2., 1., 1.])]   # traceless w.r.t. eta^{-1}=diag(-1,1,1)
    for M in cand:
        tr = -M[0, 0] + M[1, 1] + M[2, 2]
        assert abs(tr) < 1e-14, tr
        B.append(M)
    return B


def wick_terms(t1, t3, z1, z3, k, which=('s', 't', 'u'), eps13=None, full=False):
    """t1, t3: complex arrays (broadcast);  z1, z3: (..., 3).  Returns dict of the integrand g = t1^3 t3^3 f_k W (EE-contracted, and
    optionally the (v v)(w w) tensor factors for full components).  eps13: optional complex shift of x1 (x1 -> x1 + i eps13 e) used only
    to regularise the collinear (13) singularity of W_u."""
    t1 = np.asarray(t1, complex); t3 = np.asarray(t3, complex)
    x1 = -t1[..., None]*z1; x3 = t3[..., None]*z3
    x2 = -E0; x4 = E0
    if eps13 is not None:
        x1 = x1 + 1j*eps13*E0
    x12 = x1 - x2; x34 = x3 - x4; x13 = x1 - x3; x14 = x1 - x4; x23 = x2 - x3; x24 = x2 - x4
    s12, s34, s13, s14, s23, s24 = (dot(a, a) for a in (x12, x34, x13, x14, x23, x24))
    up = s12*s34/(s13*s24); vp = s14*s23/(s13*s24)
    if k == 2:
        f = vp - up
    elif k == 4:
        f = ((vp - up)**3 - (vp - up))/(4*up*vp)
    elif k == 0:
        f = np.ones_like(up)
    else:
        raise ValueError(k)
    meas = t1**3*t3**3
    out = {}
    if 's' in which:      # <T4 T3><T1 T2>
        v = inv_vec(x34, z3); w = inv_vec(x12, z1)
        den = s34**3*s12**3
        out['s'] = meas*f*dot(v, w)**2/den
        if full:
            out['s_vec'] = (v, w, meas*f/den)
    if 't' in which:      # <T4 T1><T3 T2>
        a = inv_vec(x14, z1); b = inv_vec(x23, z3)
        den = s14**3*s23**3
        out['t'] = meas*f*dot(a, b)**2/den
        if full:
            out['t_vec'] = (a, b, meas*f/den)
    if 'u' in which:      # <T4 T2><T3 T1>, T4.T2 fully contracted (EE); numerator N vanishes identically for light-ray polarisations
        N = dot(z3, z1)*s13 - 2*dot(z3, x13)*dot(z1, x13)
        Il = ETA - 2*np.outer(lower(x24), lower(x24))/s24
        Iu = np.linalg.inv(ETA)
        A = np.einsum('mr,ns,mr,ns->', Il, Il, Iu, Iu); B = np.einsum('ms,nr,mr,ns->', Il, Il, Iu, Iu)
        Cq = np.einsum('mn,rs,mr,ns->', ETA, ETA, Iu, Iu)
        T42 = (0.5*(A + B) - Cq/3.0)/s24**3
        out['u'] = meas*f*T42*N**2/s13**5
    return out


def taylor2(fun, n1, n3, r=0.2, M=32):
    """Taylor coefficients a[j,k] (j < n1, k < n3) of fun(t1, t3) at (1/2, 1/2) by a bivariate Cauchy/FFT on |t-1/2| = r.
    fun returns arrays of shape (M, M, ...)."""
    th = 2*np.pi*np.arange(M)/M
    w = r*np.exp(1j*th)
    T1 = 0.5 + w[:, None]; T3 = 0.5 + w[None, :]
    vals = fun(T1, T3)
    co = np.fft.fft2(vals, axes=(0, 1))/(M*M)
    jj = np.arange(n1); kk = np.arange(n3)
    sc = (r**(-jj))[:, None]*(r**(-kk))[None, :]
    return co[:n1, :n3]*sc.reshape(sc.shape + (1,)*(co.ndim - 2))


def edge_residue(phi1, phi3, k, term, r=0.2, M=32):
    """(-2 pi i)^2 Res_{t1=1/2} Res_{t3=1/2} of the EE-contracted integrand of Wick term 's' or 't' (u is identically 0 off phi1 = phi3).
    phi1, phi3: 1-d arrays (same length).  Pole orders: s: 6 in each variable (5 from the 2-pt function, 1 from f, k = 4), t: 1 (k = 4)."""
    z1 = zvec(phi1); z3 = zvec(phi3)
    if term == 's':
        n = 6 if k == 4 else 5
    else:
        n = 1 if k == 4 else 0
    if n == 0:
        return np.zeros(len(phi1), complex)

    def fun(T1, T3):
        T1b = T1[..., None]; T3b = T3[..., None]
        g = wick_terms(T1b, T3b, z1[None, None], z3[None, None], k, which=(term,))[term]
        return g*(2*T1b - 1)**n*(2*T3b - 1)**n
    a = taylor2(fun, n, n, r, M)
    res = a[n - 1, n - 1]/2.0**(2*n)          # Res_t [P(t)/(2t-1)^n] = P^{(n-1)}(1/2)/((n-1)! 2^n)
    return (-2j*np.pi)**2*res


def edge_residue_full(phi1, phi3, k, term, r=0.2, M=32):
    """Full (2,4) tensor of the edge term: returns array (N, 3, 3, 3, 3) of W^{mu nu}_{(4)} {}^{rho sigma}_{(2)} (upper indices),
    before projection.  Implemented through the vectors: s: v^mu v^nu w^rho w^sigma, t: a^mu a^nu b^rho b^sigma."""
    z1 = zvec(phi1); z3 = zvec(phi3)
    n = (6 if k == 4 else 5) if term == 's' else (1 if k == 4 else 0)
    if n == 0:
        return np.zeros((len(phi1), 3, 3, 3, 3), complex)

    def fun(T1, T3):
        T1b = T1[..., None]; T3b = T3[..., None]
        d = wick_terms(T1b, T3b, z1[None, None], z3[None, None], k, which=(term,), full=True)
        p, q, sc = d[term + '_vec']
        ten = np.einsum('...a,...b,...c,...d->...abcd', p, p, q, q)*sc[..., None, None, None, None]
        return ten*((2*T1b - 1)**n*(2*T3b - 1)**n)[..., None, None, None, None]
    a = taylor2(fun, n, n, r, M)
    return (-2j*np.pi)**2*a[n - 1, n - 1]/2.0**(2*n)


def celestial(phi1, phi3, phi0, nu):
    z1, z3, z0 = zvec(phi1), zvec(phi3), zvec(phi0)
    d = 0.5 + 1j*nu
    A13 = -2*dot(z1, z3); A10 = -2*dot(z1, z0); A30 = -2*dot(z3, z0)
    return A13**((2 + d)/2)/(A10*A30)**(d/2)
