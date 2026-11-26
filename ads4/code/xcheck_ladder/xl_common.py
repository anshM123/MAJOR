"""
xcheck_ladder/xl_common.py -- independent re-implementation of the exact finite-mass ladder action H_ang,b of Omega^EE
(prl/sec_fullpos.tex, Eqs. fp-U, fp-zy, fp-moments, fp-H, fp-single), written from the TeX, not from code/fullpos.

Two sources of the moments Z_k, Y_k (k = 0..K):
  'cert'  : the certified Arb moments of code/fullpos/cache/moments_<Dg>_K1000.json (read-only), midpoints;
  'quad'  : moments of a DISCRETE measure  sum_i w_i g(nu_i) delta(nu - nu_i)  (Gauss-Legendre on the continuum weight, or the
            published midpoint grid), accumulated in Arb at high precision.  Because the action is linear in the measure,
            H[quad] = sum_i w_i g_i K_T(nu_i) EXACTLY, where K_T(nu) is the exact per-nu finite-T ladder kernel; K_T is bounded on
            nu < m, so the quadrature error of H is that of an ordinary smooth integral, NOT amplified by the 10^k cancellation of
            the ladder.  (The cancellation is handled by the Arb precision of the moment accumulation; checked by varying prec.)
Nothing outside code/xcheck_ladder/ is written.
"""
import os, sys, json, math
sys.dont_write_bytecode = True
import numpy as np
from fractions import Fraction as Fr
from flint import arb, fmpq, fmpq_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
OUT = os.path.join(HERE, 'out')
K0 = Fr('27.48509817')


def lam_int(Dg):
    return (Dg - 2)*(Dg - 1)


# ------------------------------------------------------------------ NF02 weights (exact)
def nf02_exact():
    nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    p2 = [Fr(c)/K0 for c in nf['psi2_exact']]
    p3 = [Fr(c)/K0 for c in nf['psi3_exact']]
    return p2, p3


# ------------------------------------------------------------------ CDH polynomials (Eq. fp-zy and the d=1 recurrence), exact in x = nu^2
def cdh_f(M):
    """f_0..f_M: f_0 = 0, f_1 = 1, (m+1)^2 f_{m+1} = -(s - 2 m^2) f_m - (m-1)^2 f_{m-1}, s = x + 1/4."""
    X = fmpq_poly([0, 1])
    f = [fmpq_poly([0]), fmpq_poly([1])]
    for m in range(1, M):
        f.append((-(X + fmpq(1, 4) - 2*m*m)*f[m] - (m - 1)**2*f[m - 1])*fmpq(1, (m + 1)**2))
    return f


def cdh_q(N):
    """q_0..q_N: q_0 = 1, (n+4)^2 q_{n+1} = -(s + 4 - 2(n+2)^2) q_n - n^2 q_{n-1}."""
    X = fmpq_poly([0, 1])
    q = [fmpq_poly([1])]
    for n in range(N):
        prev = q[n - 1] if n >= 1 else fmpq_poly([0])
        q.append((-(X + fmpq(1, 4) + 4 - 2*(n + 2)**2)*q[n] - n*n*prev)*fmpq(1, (n + 4)**2))
    return q


def zy_polys(K):
    """z_{k+2}, y_{k+2} (k = 0..K) as exact fmpq_poly in x."""
    f = cdh_f(K + 6)
    Z, Y = [], []
    for k in range(K + 1):
        m = k + 2
        Z.append(12*f[m] - (m - 1)*f[m - 1] + (m + 1)*f[m + 1])
        Y.append(fmpq(12 - m*m, 2)*f[m] + fmpq(13 - 7*m, 2)*f[m - 1] + fmpq(13 + 7*m, 2)*f[m + 1]
                 + fmpq((m - 2)**2, 4)*f[m - 2] + fmpq((m + 2)**2, 4)*f[m + 2])
    return Z, Y


def Q_exact(K, lam):
    """Q_k = int_0^1 psi3(p) q_k(lam p^2) dp (exact)."""
    _, p3 = nf02_exact()
    p3 = [fmpq(c.numerator, c.denominator) for c in p3]
    mom = [sum((c/(2*i + n + 1) for n, c in enumerate(p3)), fmpq(0)) for i in range(K + 2)]
    q = cdh_q(K)
    out = []
    for k in range(K + 1):
        cs = q[k].coeffs()
        out.append(sum((c*fmpq(lam)**i*mom[i] for i, c in enumerate(cs)), fmpq(0)))
    return out


# ------------------------------------------------------------------ the continuum weight (Dg = 300, 1000), double precision
def load_cont(Dg):
    d = json.load(open(os.path.join(DATA, 'k2', f'k2_continuum_D_{Dg}_rep.json')))
    cs = np.array([[float(c), float(s)] for c, s in d['centres_sigma']])
    x = np.array([float(v) for v in d['x_re']]); y = np.array([float(v) for v in d['y_re']])
    return dict(c=cs[:, 0], s=cs[:, 1], x=x, y=y, Lam=math.sqrt(lam_int(Dg)), Dg=Dg)


def chi_re(W, nu):
    """chi_re(nu) = Q(nu) sum_k x_k [G_k(nu-c_k) + G_k(nu+c_k)] (only G(nu) for c_k = 0) + sum_j y_j nu^{2j} l(nu) e^{-nu^2/2}."""
    nu = np.asarray(nu, float)
    out = np.zeros_like(nu)
    c, s, x = W['c'], W['s'], W['x']
    for i0 in range(0, len(nu), 2000):
        v = nu[i0:i0 + 2000][:, None]
        B = x*(np.exp(-(v - c)**2/(2*s*s)) + np.where(c > 0, np.exp(-(v + c)**2/(2*s*s)), 0.0))
        out[i0:i0 + 2000] = B.sum(1)
    x2 = nu*nu
    Q = x2/(x2 + 0.25)**2
    L = (W['y'][0] + W['y'][1]*x2 + W['y'][2]*x2*x2)*(x2 + 1)*np.log1p(x2)*np.exp(-x2/2)
    return Q*out + L


def theta(nu, Lam):
    from scipy.special import erfc
    nuc = Lam + 100; a = 10.0
    # (1/2)[erf((nu+nuc)/a) - erf((nu-nuc)/a)] = 1 - (1/2)erfc((nu+nuc)/a) - (1/2)erfc((nuc-nu)/a)
    return 1 - 0.5*erfc((nu + nuc)/a) - 0.5*erfc((nuc - nu)/a)


def g_of_chi(chi, nu, Lam):
    x = nu*nu
    return -9216.0/Lam*(4*x + 1)/((4*x + 9)*(4*x + 49))*chi


def cont_nodes(Lam, nup=24, numax_extra=70.0):
    """Gauss-Legendre panels on [0, Lam + numax_extra]: width 0.125 on [0,1], 0.25 on [1,4], 0.5 on [4,40], 1.0 beyond."""
    edges = list(np.arange(0, 1, 0.125)) + list(np.arange(1, 4, 0.25)) + list(np.arange(4, 40, 0.5)) + \
        list(np.arange(40, Lam + numax_extra, 1.0)) + [Lam + numax_extra]
    t, w = np.polynomial.legendre.leggauss(nup)
    N, Wt = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        N.append((a + b)/2 + (b - a)/2*t); Wt.append((b - a)/2*w)
    return np.concatenate(N), np.concatenate(Wt)


def grid_measure(Dg):
    """published grid weight: nu_k, dnu, chi_re = Re X / (R_- + R_+), read as zero beyond Lam + 10 (the *_trunc files)."""
    import gzip
    if Dg in (300, 1000):
        fn = os.path.join(DATA, f'd_abel_X_{Dg}_trunc.txt'); op = open
    else:
        fn = os.path.join(DATA, 'gap', f'd_abel_X_{Dg}_trunc.txt.gz'); op = gzip.open
    rows = []
    with op(fn, 'rt') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            rows.append([float(v) for v in line.split()[:2]])
    a = np.array(rows)
    nu, ReX = a[:, 0], a[:, 1]
    dnu = float(np.round(nu[1] - nu[0], 6))
    # R_- + R_+ = 2 Re R_-  (kernels.json D_k2_EE; R_+ = conj R_- on the real axis)
    z = nu.astype(complex)
    I = 1j
    Rm = -2*(2*z - I)*(2*z + I)**2*(2*z + 3*I)*(2*z + 5*I)*(2*z + 11*I)/(z*(2*z - 7*I)*(2*z - 3*I)*(2*z + 7*I))
    Rp = -2*(2*z - 11*I)*(2*z - 5*I)*(2*z - 3*I)*(2*z - I)**2*(2*z + I)/(z*(2*z - 7*I)*(2*z + 3*I)*(2*z + 7*I))
    chi = ReX/(Rm + Rp).real
    return nu, np.full_like(nu, dnu), chi


# ------------------------------------------------------------------ discrete-measure power moments in Arb
def power_moments(nu, wg, J2max, prec, ranges=None):
    """S_j = sum_i wg_i nu_i^{2j}, j = 0..J2max, in Arb at `prec` bits (nu_i, wg_i exact doubles).
    ranges: list of (lo, hi) nu-intervals -> one S-vector per range."""
    ctx.prec = prec
    if ranges is None:
        ranges = [(-1, 1e30)]
    res = []
    for lo, hi in ranges:
        S = [arb(0)]*(J2max + 1)
        sel = np.nonzero((nu > lo) & (nu <= hi))[0]
        for i in sel:
            x = arb(float(nu[i]))**2
            t = arb(float(wg[i]))
            for j in range(J2max + 1):
                S[j] += t
                t *= x
        res.append(S)
    return res


def moments_from_power(S, K):
    """Z_k = sum_j zcoef_{k,j} S_j etc."""
    Zp, Yp = zy_polys(K)
    Z = [sum((arb(c)*S[j] for j, c in enumerate(p.coeffs())), arb(0)) for p in Zp]
    Y = [sum((arb(c)*S[j] for j, c in enumerate(p.coeffs())), arb(0)) for p in Yp]
    return Z, Y


def load_cert_moments(Dg=1000, K=1000, kmax=None):
    d = json.load(open(os.path.join(ADS4, 'code', 'fullpos', 'cache', f'moments_{Dg}_K{K}.json')))
    n = len(d['Z']) if kmax is None else kmax + 1
    Z = [arb(s) for s in d['Z'][:n]]; Y = [arb(s) for s in d['Y'][:n]]
    return Z, Y


# ------------------------------------------------------------------ the ladder action, Eq. (fp-H) and Eq. (fp-single)
def H_ang(J, b, T, Z, Y, Q, lam, parts=False):
    """H_ang,b = (1-r)^3 sum_{k<=n} D_k^2 U_{b,k} [Z_k(1/(k+3) - r/(k+1)) + (2Y_k + V_{b,k}Z_k + sigma_b Lam^2 Q_k)/T]."""
    T = arb(T) if not isinstance(T, arb) else T
    C = J*(J + 1); r = arb(C)/T; n = J - b
    sig = 1 if b == 0 else -1
    U = arb(1); h2 = arb(0); h3 = arb(0)
    for k in range(n + 1):
        Dk2 = math.comb(k + 3, 3)**2
        V = arb(fmpq(4*k, k + 3)) if b == 0 else 4 + arb(fmpq(12, k + 1))
        h2 += Dk2*U*(Z[k]*(1/arb(k + 3) - r/(k + 1)) + (2*Y[k] + V*Z[k])/T)
        h3 += Dk2*U*(sig*arb(lam)*arb(Q[k]))/T
        if k < n:
            uk = (b + k)*(b + k + 1)
            U = U*(C - uk)/(T - uk)
    w3 = (1 - r)**3
    if parts:
        return w3*h2, w3*h3
    return w3*(h2 + h3)


def H2_single(J, b, T, Z):
    """Eq. (fp-single): H^(2) = (1-r)^3 sum_{k=0}^{n+1} V_{b,k} P_k, P_k = D_k^2 L_k, L_k = (Z_k - k/(k+3) Z_{k-1})/(k+3)."""
    T = arb(T) if not isinstance(T, arb) else T
    C = J*(J + 1); r = arb(C)/T; n = J - b
    Us = [arb(1)]
    for k in range(n + 1):
        uk = (b + k)*(b + k + 1)
        Us.append(Us[-1]*(C - uk)/(T - uk))       # Us[n+1] = 0 exactly (C = u_n)
    tot = arb(0)
    for k in range(n + 2):
        L = Z[0]/3 if k == 0 else (Z[k] - arb(fmpq(k, k + 3))*Z[k - 1])/(k + 3)
        P = math.comb(k + 3, 3)**2*L
        ukm1 = (b + k - 1)*(b + k)                  # u_{k-1}, continued to k = 0: (b-1) b
        V = Us[k]*(1 + (arb(fmpq((k + 3)**2, 2)) - ukm1)/T) + (arb(fmpq(k*k, 2))/T*Us[k - 1] if k >= 1 else 0)
        tot += V*P
    return (1 - r)**3*tot


def T_of(Delta=None, mfl=None, J=None, Dg=1000):
    """T = (Delta-1)(Delta-2) (exact Fraction for integer Delta), or T = m^2 + J(J+1) with m = mfl*Lam (mfl a decimal string)."""
    if Delta is not None:
        D = Fr(Delta)
        return (D - 1)*(D - 2)
    lam = lam_int(Dg)
    return Fr(mfl)**2*lam + J*(J + 1)


def to_arb(fr):
    return arb(fmpq(fr.numerator, fr.denominator))
