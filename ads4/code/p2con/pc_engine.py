"""
p2con/pc_engine.py (P2_CONTRACT.md §1): EXACT leading-Regge EE plane-wave image engine for ARBITRARY constant polarisation tensors
and ARBITRARY polynomial insertions -- i.e. the full rank-3 (three-channel) leading-Regge EE kernel, no helicity model.

Kernel: the split representation of H7 §18.1 (h7i14/v3_regge_kernel.py; read-only import of its monomial lists).  The raw (eta-kernel,
f-dressed) image of a plane wave e^{p.x + pb.y} times polarisation tensors and a polynomial Q(x) R(y) is
      raw = - eta^{mu nu} d/dp^mu d/dpb^nu  int Dz0  Ahat_6(p, z0)[z2', z4'; Q]  Bhat_6(pb, z0)[z1, z3; R]          (b = a - 1 = 6)
with Ahat = sum_monomials coef(z0.e, z0.f, e.f) K(c, lam) Q(d_p) (e.d)^i (f.d)^j [S^-lam R^kap](p) (v3 docstring), similarly Bhat.
The no-gradient version with b = 7 is the 'f-stripped' image (Dp = 1 kernel; used only in checks).

Implementation (no quadrature anywhere except the final celestial integral, which is done exactly below):
  * each side is evaluated in the REST FRAME of its momentum (p0 for side A = pair (2,4), pb0 for side B = pair (1,3)); there
    S = -z0.p0 = 1 on the parametrised celestial circle z0 = (1, cos phi, sin phi), so the side function is an exact trigonometric
    polynomial in phi; its harmonics A_k are obtained by an exact DFT (32 points, |k| <= 12 needed);
  * derivatives: the full 3-variable Taylor jet of S^-lam R^kap to total order 8 (exact truncated power series), contracted with the
    operator polynomial (e.d)^i (f.d)^j Q(d) d_a;
  * polarisation tensors on each leg: null-square basis n(t_j) (x) n(t_j), t_j = 2 pi j/5 (the five-point basis of traceless symmetric
    rank-2 tensors in 3d); a mode e^{ikt} (P2_SPINNING's helicity-frame modes) has weights w_j = e^{ikt_j}/(5 c_|k|), c = (3/2, -1, 1/4);
  * the celestial pairing of harmonic k (frame p0) with harmonic k' (frame pb0) is the SO(2,1) matrix element
        M_kk'(nu, beta) = int dpsi S_A^{-(1-delta)} S_B^{-delta} e^{i k phi_A} e^{i k' phi_B}       (delta = 1/2 + i nu)
    computed (i) by the trapezoid rule on the global circle (spectral; used for nu beta <~ 40) and (ii) in closed form
        M_kk' = 2pi sqrt2^{n+n'} / ((1-delta)_n (delta)_n') sum_j C(n,j) C(n',j) j! eta_AB^j eta_A^{n-j} eta_B^{n'-j} P_l^{(n+n'-j)}(eta),
    n = |k|, eta_A = -eps_A.pb0, eta_B = -eps_B.p0, eta_AB = -eps_A.eps_B  (used for nu beta >~ 40, where the P, P' reduction is stable).
All nu-only normalisations (K^ref of each side) cancel in ratios to the leading image computed with the same engine.
"""
import os, sys, functools, math
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14')); sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
import numpy as np
import sympy as sp
from scipy.special import loggamma
import v3_regge_kernel as V          # read-only (monomial lists)
import v4_hclosed as H4               # read-only (sympy monomials)
import s1_offdiag_kernel as S1        # read-only (the EE (2,4) structure V0^4 (z2'.z4')^2, A_monos_struct('EE'))

NORD = 8
IDX = [(a, b, d - a - b) for d in range(NORD + 1) for a in range(d, -1, -1) for b in range(d - a, -1, -1)]
POS = {t: i for i, t in enumerate(IDX)}
NM = len(IDX)
DEG = np.array([sum(t) for t in IDX])
FACT = np.array([math.factorial(a)*math.factorial(b)*math.factorial(c) for a, b, c in IDX], float)
_PROD = []
for i, ta in enumerate(IDX):
    js, ks = [], []
    for j, tb in enumerate(IDX):
        s = (ta[0] + tb[0], ta[1] + tb[1], ta[2] + tb[2])
        if sum(s) <= NORD:
            js.append(j); ks.append(POS[s])
    _PROD.append((np.array(js), np.array(ks)))
MET = np.array([-1.0, 1.0, 1.0])
NPHI = 32
PHI = 2*np.pi*np.arange(NPHI)/NPHI
ZREST = np.stack([np.ones(NPHI), np.cos(PHI), np.sin(PHI)], -1)      # z0 in a rest frame (S = 1)
KH = 12                                                             # harmonics |k| <= KH kept (checked: higher ones vanish)
TJ = 2*np.pi*np.arange(5)/5
CK = {0: 1.5, 1: -1.0, 2: 0.25}


def mdot(u, v):
    return -u[..., 0]*v[..., 0] + u[..., 1]*v[..., 1] + u[..., 2]*v[..., 2]


def jmul(U, V_):
    """product of jets (..., NM) (broadcasting)."""
    sh = np.broadcast_shapes(U.shape, V_.shape)
    out = np.zeros(sh, complex)
    for i, (js, ks) in enumerate(_PROD):
        u = U[..., i]
        if not np.any(u):
            continue
        out[..., ks] += u[..., None]*V_[..., js]
    return out


def poly_lin(vec):
    """jet of the linear form sum_a vec[a] t^a (vec (..., 3))."""
    out = np.zeros(vec.shape[:-1] + (NM,), complex)
    out[..., POS[(1, 0, 0)]] = vec[..., 0]; out[..., POS[(0, 1, 0)]] = vec[..., 1]; out[..., POS[(0, 0, 1)]] = vec[..., 2]
    return out


def one():
    o = np.zeros(NM, complex); o[0] = 1; return o


def binom_arr(x, n):
    r = np.ones_like(x, dtype=complex)
    for k in range(n):
        r = r*(x - k)/(k + 1)
    return r


@functools.lru_cache(maxsize=None)
def _hpow():
    """powers of h_S (per phi) and h_R: S = 1 + h_S, R = 1 + h_R in rest coordinates t around p0 = (1,0,0)."""
    hS = poly_lin(np.stack([np.ones(NPHI), -np.cos(PHI), -np.sin(PHI)], -1))      # S = -z0.(p0+t) = 1 + t0 - cos t1 - sin t2
    hR = np.zeros(NM, complex)                                                      # R = -(p0+t)^2 = 1 + 2t0 + t0^2 - t1^2 - t2^2
    hR[POS[(1, 0, 0)]] = 2; hR[POS[(2, 0, 0)]] = 1; hR[POS[(0, 2, 0)]] = -1; hR[POS[(0, 0, 2)]] = -1
    pS = [np.broadcast_to(one(), (NPHI, NM)).copy()]; pR = [one()]
    for n in range(NORD):
        pS.append(jmul(pS[-1], hS)); pR.append(jmul(pR[-1], hR))
    return pS, pR


def jet_SR(lam, kap):
    """derivative tensor D[nu, phi, gamma] = gamma! * Taylor coefficient of S^-lam R^kap (lam, kap arrays over nu)."""
    pS, pR = _hpow()
    lam = np.atleast_1d(lam); kap = np.atleast_1d(kap)
    JS = sum(binom_arr(-lam, n)[:, None, None]*pS[n][None] for n in range(NORD + 1))
    JR = sum(binom_arr(kap, n)[:, None]*pR[n][None] for n in range(NORD + 1))
    return jmul(JS, JR[:, None, :])*FACT


def opoly_dir(v):
    """jet-index polynomial of the directional derivative v.d = sum_a v^a d_a (v: 3 upper rest components)."""
    return poly_lin(np.asarray(v, complex))


def opoly_pow(P, n):
    r = one()
    for _ in range(n):
        r = jmul(r, P)
    return r


def opoly_ins(alpha):
    """insertion x^alpha (rest contravariant components) -> operator prod (eta^{bb} d_b)^alpha_b."""
    r = one()
    for b in range(3):
        for _ in range(alpha[b]):
            v = np.zeros(3); v[b] = MET[b]
            r = jmul(r, poly_lin(v))
    return r


def nvec(t):
    return np.array([1.0, np.cos(t), np.sin(t)])


# --------------------------------------------------------------------------------------------- monomial groups of each side
def _lk(c, lam):
    mu = c + 3 - lam
    return (mu - 1)*np.log(2.0) + 0.5*np.log(np.pi) + loggamma(1 + c/2) + loggamma((c + 3)/2 - lam)


ASTRUCT = 'EE'          # the EE functional's (2,4) structure (EE_FUNCTIONAL); 'q2' = v3's q^2 (checks only)


@functools.lru_cache(maxsize=None)
def side_groups(side, spin, astruct=None):
    """[(k(delta), m(delta), [(coef_fun(arg1, arg2, arg3, delta), i, j)])] grouped by (k, m) (v3 / s1 monomial lists)."""
    ast = ASTRUCT if astruct is None else astruct
    if side == 'A':
        monos = (S1.A_monos_struct(ast) if ast != 'q2' else V.A_monomials(True)) if spin else V.A_monomials(False)
    else:
        monos = V.B_monomials(spin)
    dt = 0.37 + 1.91j
    groups = {}
    for (cf, kf, mf, i, j) in monos:
        key = (round(complex(kf(dt)).real, 9), round(complex(kf(dt)).imag, 9), round(complex(mf(dt)).real, 9), round(complex(mf(dt)).imag, 9))
        groups.setdefault(key, [kf, mf, []])[2].append((cf, i, j))
    return [tuple(v) for v in groups.values()]


def side_table(side, spin, nus, pol_pairs, alphas, grad, b=6, astruct=None):
    """Harmonic table of one side.
    pol_pairs: list of (e, f) rest-frame vectors (A: (z2', z4'); B: (z1, z3));  alphas: list of insertion multi-indices (rest frame);
    grad: list of gradient components (None = no gradient, or 0, 1, 2).
    Returns T[nu, ipol, ialpha, igrad, k] for k = -KH..KH (DFT harmonics in the rest-frame angle), with the side's nu-only factor
    K_ref removed (K_ref = K of the scalar group: A: (c, lam) = (b-3 + 1-delta, 1-delta); B: (b-3 + delta, delta))."""
    nus = np.atleast_1d(np.asarray(nus, complex))
    delta = 0.5 + 1j*nus
    if side == 'A':
        lref = _lk((b - 3) + (1 - delta), 1 - delta)
    else:
        lref = _lk((b - 3) + delta, delta)
    groups = side_groups(side, spin, astruct)
    npol, nal, ngr = len(pol_pairs), len(alphas), len(grad)
    out = np.zeros((len(nus), npol, nal*ngr, NPHI), complex)
    insP = [opoly_ins(a) for a in alphas]
    grP = [one() if g is None else opoly_dir(np.eye(3)[g]) for g in grad]
    IG = np.array([jmul(insP[ia], grP[ig]) for ia in range(nal) for ig in range(ngr)])        # (nal*ngr, NM)
    for (kf, mf, terms) in groups:
        k = np.asarray(kf(delta), complex)*np.ones_like(delta); m = np.asarray(mf(delta), complex)*np.ones_like(delta)
        c = (b - 3) + k; lam = -m; kap = (2*lam - c - 3)/2
        Kr = np.exp(_lk(c, lam) - lref)                                     # (nu,)
        D = jet_SR(lam, kap)                                                # (nu, phi, NM)
        ij = sorted(set((i, j) for (_, i, j) in terms))
        for ip, (e, f) in enumerate(pol_pairs):
            z = ZREST
            a1, a2, a3 = mdot(z, e[None, :]*np.ones((NPHI, 1))), mdot(z, f[None, :]*np.ones((NPHI, 1))), np.full(NPHI, mdot(e, f))
            eP = opoly_dir(e); fP = opoly_dir(f)
            for (i, j) in ij:
                W = np.zeros((len(nus), NPHI), complex)
                for (cf, ii, jj) in terms:
                    if (ii, jj) != (i, j):
                        continue
                    W += np.asarray(cf(a1[None, :], a2[None, :], a3[None, :], delta[:, None]), complex)*np.ones((len(nus), NPHI))
                if not np.any(W):
                    continue
                base = jmul(opoly_pow(eP, i), opoly_pow(fP, j))
                Qs = jmul(base[None, :], IG)                                # (nal*ngr, NM)
                val = np.einsum('npg,qg->npq', D, Qs)                       # (nu, phi, q)
                out[:, ip] += np.transpose(Kr[:, None, None]*W[:, :, None]*val, (0, 2, 1))
    out = out.reshape(len(nus), npol, nal, ngr, NPHI)
    H = np.fft.fft(out, axis=-1)/NPHI                                       # harmonic k at index k mod NPHI
    ks = np.arange(-KH, KH + 1)
    return H[..., ks % NPHI], np.abs(H[..., (np.arange(KH + 1, NPHI - KH)) % NPHI]).max() if NPHI > 2*KH + 1 else 0.0


# --------------------------------------------------------------------------------------------- frames, pairing
def boost(r):
    ch, sh = np.cosh(r), np.sinh(r)
    return np.array([[ch, sh, 0.0], [sh, ch, 0.0], [0.0, 0.0, 1.0]])


def frames(beta):
    """rest bases (columns) of side A (p0 = B(beta/2) e0) and side B (pb0 = B(-beta/2) e0) in global Cartesian coordinates."""
    return boost(beta/2), boost(-beta/2)


def gradient_metric(beta):
    """G^{ab} = - eta^{aa} eta^{bb} (E_a . F_b): raw = sum_ab G^{ab} (d_a Ahat)(d_b Bhat)."""
    EA, EB = frames(beta)
    G = np.zeros((3, 3))
    for a in range(3):
        for c in range(3):
            G[a, c] = -MET[a]*MET[c]*mdot(EA[:, a], EB[:, c])
    return G


def M_quad(nus, beta, npsi=None):
    """M[nu, k, k'] by the trapezoid rule on the global celestial circle (k, k' = -KH..KH)."""
    nus = np.atleast_1d(np.asarray(nus, complex))
    if npsi is None:
        npsi = int(max(64, 2**np.ceil(np.log2(8*np.max(np.abs(nus))*beta + 64))))
    EA, EB = frames(beta)
    psi = 2*np.pi*np.arange(npsi)/npsi
    z = np.stack([np.ones(npsi), np.cos(psi), np.sin(psi)], -1)
    SA = -mdot(z, EA[:, 0][None, :]); SB = -mdot(z, EB[:, 0][None, :])
    eA = (mdot(z, EA[:, 1][None, :]) + 1j*mdot(z, EA[:, 2][None, :]))/SA    # e^{i phi_A}
    eB = (mdot(z, EB[:, 1][None, :]) + 1j*mdot(z, EB[:, 2][None, :]))/SB
    ks = np.arange(-KH, KH + 1)
    YA = eA[None, :]**ks[:, None]; YB = eB[None, :]**ks[:, None]           # (k, psi)
    delta = 0.5 + 1j*nus
    W = np.exp(-(1 - delta)[:, None]*np.log(SA)[None, :] - delta[:, None]*np.log(SB)[None, :])*(2*np.pi/npsi)
    return np.einsum('np,kp,lp->nkl', W, YA, YB)


def legendre_PdP(nus, beta):
    """P_l(cosh beta) and dP_l/deta for l = -1/2 + i nu (complex nu ok): Mehler-Dirichlet (n2_common.Pnu) and
    (eta^2 - 1) P_l' = l (eta P_l - P_{l-1})."""
    sys.path.insert(0, os.path.join(HERE, '..', 'n2sub'))
    import n2_common as N
    nus = np.atleast_1d(np.asarray(nus, complex))
    eta = np.cosh(beta)
    nph = int(600 + 3*np.max(np.abs(nus))*beta)
    P = N.Pnu(beta, nus, nphi=nph); Pm = N.Pnu(beta, nus + 1j, nphi=nph)
    l = -0.5 + 1j*nus
    return P, l*(eta*P - Pm)/(eta*eta - 1)


def M_closed(nus, beta):
    """closed form of M[nu, k, k'] reduced to (coefP, coefdP): M = cP P_l + cD P_l'  (stable for nu beta >> |k| + |k'|)."""
    nus = np.atleast_1d(np.asarray(nus, complex))
    delta = 0.5 + 1j*nus; l = delta - 1
    eta = np.cosh(beta)
    EA, EB = frames(beta)
    epsA = {s: (EA[:, 1] + s*1j*EA[:, 2])/np.sqrt(2) for s in (1, -1)}
    epsB = {s: (EB[:, 1] + s*1j*EB[:, 2])/np.sqrt(2) for s in (1, -1)}
    p0, pb0 = EA[:, 0], EB[:, 0]
    mmax = 2*KH
    c = [np.ones_like(nus), np.zeros_like(nus)]; d = [np.zeros_like(nus), np.ones_like(nus)]
    for k in range(mmax - 1):
        fac = (k*(k + 1) - l*(l + 1))
        c.append(-(2*(k + 1)*eta*c[k + 1] + fac*c[k])/(eta*eta - 1))
        d.append(-(2*(k + 1)*eta*d[k + 1] + fac*d[k])/(eta*eta - 1))
    ks = np.arange(-KH, KH + 1)
    cP = np.zeros((len(nus), len(ks), len(ks)), complex); cD = np.zeros_like(cP)
    for ia, k in enumerate(ks):
        n = abs(k); sA = 1 if k >= 0 else -1
        for ib, kk in enumerate(ks):
            nn = abs(kk); sB = 1 if kk >= 0 else -1
            eA, eB = epsA[sA], epsB[sB]
            hA = -mdot(eA, pb0); hB = -mdot(eB, p0); hAB = -mdot(eA, eB)
            pre = 2*np.pi*np.sqrt(2)**(n + nn)/(poch(1 - delta, n)*poch(delta, nn))
            for j in range(min(n, nn) + 1):
                w = math.comb(n, j)*math.comb(nn, j)*math.factorial(j)*hAB**j*hA**(n - j)*hB**(nn - j)
                mm = n + nn - j
                cP[:, ia, ib] += pre*w*c[mm]; cD[:, ia, ib] += pre*w*d[mm]
    return cP, cD


def poch(x, n):
    r = np.ones_like(x)
    for k in range(n):
        r = r*(x + k)
    return r


# --------------------------------------------------------------------------------------------- polarisation-mode bases
def mode_weights(k):
    """weights w_j of the null-square basis n(t_j) (x) n(t_j) for the leg tensor whose null-circle function is e^{ikt}."""
    return np.exp(1j*k*TJ)/(5*CK[abs(k)])


POLS5 = [nvec(t) for t in TJ]
POLPAIRS = [(POLS5[a], POLS5[b]) for a in range(5) for b in range(5)]         # (first leg j, second leg j') -> index 5a + b
MODES = [-2, -1, 0, 1, 2]


def to_modes(T):
    """T[..., ipolpair(25), rest...] (pol-grid) -> T[..., kfirst(5), ksecond(5), rest...] (modes)."""
    W = np.array([mode_weights(k) for k in MODES])                          # (5 modes, 5 grid)
    sh = T.shape
    T5 = T.reshape(sh[:-4] + (5, 5) + sh[-3:]) if False else None
    return W


ALPHAS3 = [t for t in IDX if sum(t) <= 3]
