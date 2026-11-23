r"""
lrc_blocks (LRC_CLOSE.md §2): exact d = 3 Casimir eigenfunctions with CONTINUOUS spin J' = -1/2 + i nu, by the radial Casimir
recursion of reggemoment/blocks3d.py (read, not imported; same operator, generalised to complex j and to all offsets m in [-n, n]).

Normalisation:  Gh(r, eta; nu) = sum_{n even} r^n sum_{m = -n..n step 2} B[n, m](nu) P_{J'+m}(eta),   B[0, 0] = 1,
so that Gh -> P_{J'}(eta) as r -> 0 at fixed eta (leading Regge).  The block itself is (4r)^{Dp} Gh; Dp is its dimension.
(r, eta) are the radial coordinates of the (13)(42) channel; on the Regge box rho = r e^beta (larger), rhob = r e^-beta, eta = cosh beta.

Recursion (blocks3d.py header, d = 3, equal external dimensions): with s = r^2, K(E, j) = E(E-3) + j(j+1) - C,
  s^0: K/2;  s^1: (K/2 + E)(1 - 4 eta^2) + 4 eta(1-eta^2) d_eta;  s^2: (K/2 + 2E)(4 eta^2 - 1) - 4 eta(1-eta^2) d_eta;  s^3: -(K/2 + 3E),
acting on r^E P_j;  eta^2 P_j and eta(1-eta^2) P_j' re-expanded in P_{j-2}, P_j, P_{j+2} (valid for complex j: no truncation).

Legendre functions of complex degree on eta = cosh beta >= 1:  P_{J'+k}, k >= 0, by the (stable, dominant) forward recurrence from
mpmath.legenp starting values, stored scaled  Pt[k] = P_{J'+k} e^{-k beta};  P_{J'-k} = P_{-1-J'+k} = conj(P_{J'+k}) for real nu.
"""
import numpy as np
import mpmath as mp


def legendre_coeffs(j):
    """complex j (array): eta^2 P_j = e_up P_{j+2} + e_mid P_j + e_dn P_{j-2};  eta(1-eta^2)P_j' = d_up P_{j+2} + d_mid P_j + d_dn P_{j-2}."""
    j = np.asarray(j, dtype=complex)
    e_up = (j + 1)*(j + 2)/((2*j + 1)*(2*j + 3))
    e_mid = (j + 1)**2/((2*j + 1)*(2*j + 3)) + j**2/((2*j + 1)*(2*j - 1))
    e_dn = j*(j - 1)/((2*j + 1)*(2*j - 1))
    c = j*(j + 1)/(2*j + 1)
    d_up = -c*(j + 2)/(2*j + 3)
    d_mid = c*j/(2*j - 1) - c*(j + 1)/(2*j + 3)
    d_dn = c*(j - 1)/(2*j - 1)
    return e_up, e_mid, e_dn, d_up, d_mid, d_dn


def block_coeffs(Dp, nu, N):
    """B[lev, i], n = 2 lev (lev = 0..N/2), m = -N + 2 i (i = 0..N).  Complex spin J' = -1/2 + i nu, dimension Dp."""
    Jp = -0.5 + 1j*nu
    nlev = N//2 + 1
    m = -N + 2*np.arange(N + 1)
    jj = Jp + m
    C = Dp*(Dp - 3) + Jp*(Jp + 1)
    e_up, e_mid, e_dn, d_up, d_mid, d_dn = legendre_coeffs(jj)
    B = np.zeros((nlev, N + 1), dtype=complex)
    i0 = N//2
    B[0, i0] = 1.0

    def apply(vec, a2, aD, a1):
        out = (a2*e_mid + aD*d_mid + a1)*vec
        t = (a2*e_up + aD*d_up)*vec
        out[1:] += t[:-1]
        t = (a2*e_dn + aD*d_dn)*vec
        out[:-1] += t[1:]
        return out

    for lev in range(1, nlev):
        n = 2*lev
        rhs = np.zeros(N + 1, dtype=complex)
        E = Dp + n - 2
        K = E*(E - 3) + jj*(jj + 1) - C
        rhs += apply(B[lev - 1], -4*(K/2 + E), 4.0, (K/2 + E))
        if lev >= 2:
            E = Dp + n - 4
            K = E*(E - 3) + jj*(jj + 1) - C
            rhs += apply(B[lev - 2], 4*(K/2 + 2*E), -4.0, -(K/2 + 2*E))
        if lev >= 3:
            E = Dp + n - 6
            K = E*(E - 3) + jj*(jj + 1) - C
            rhs += -(K/2 + 3*E)*B[lev - 3]
        E = Dp + n
        K = E*(E - 3) + jj*(jj + 1) - C
        lo, hi = i0 - lev, i0 + lev
        B[lev, lo:hi + 1] = -rhs[lo:hi + 1]/(K[lo:hi + 1]/2)
    return B, m


def legendre_scaled(nu, beta, kmax, dps=30):
    """Pt[k] = P_{J'+k}(cosh beta) e^{-k beta}, k = 0..kmax, J' = -1/2 + i nu."""
    with mp.workdps(dps):
        Jp = mp.mpf(-0.5) + 1j*mp.mpf(nu)
        z = mp.cosh(mp.mpf(beta))
        p0 = complex(mp.legenp(Jp, 0, z, type=3))
        p1 = complex(mp.legenp(Jp + 1, 0, z, type=3))*float(mp.e**(-mp.mpf(beta)))
    eta = np.cosh(beta); emb = np.exp(-beta)
    Pt = np.empty(kmax + 1, dtype=complex)
    Pt[0] = p0
    if kmax >= 1:
        Pt[1] = p1
    Jp = -0.5 + 1j*nu
    for k in range(1, kmax):
        mu = Jp + k
        Pt[k + 1] = ((2*mu + 1)*eta*emb*Pt[k] - mu*emb**2*Pt[k - 1])/(mu + 1)
    return Pt


def radial_Y(B, m, Pt, beta):
    """Y[lev] = sum_m B[lev, m] e^{-(n - |m|) beta} Pt_{|m|} (conj for m < 0), so that Gh(rho, beta) = sum_lev Y[lev] rho^n."""
    nlev, nm = B.shape
    am = np.abs(m)
    P = np.where(m >= 0, Pt[am], np.conj(Pt[am]))
    n = 2*np.arange(nlev)[:, None]
    with np.errstate(over='ignore', invalid='ignore'):
        fac = np.where(am[None, :] <= n, np.exp(-(n - am[None, :])*beta), 0.0)
    return (B*fac*P[None, :]).sum(axis=1)


def Gh(Y, rho):
    """sum_lev Y[lev] rho^{2 lev} (rho may be an array)."""
    rho = np.atleast_1d(rho)
    nlev = len(Y)
    pw = rho[:, None]**(2*np.arange(nlev))[None, :]
    return pw @ Y


def edge_F(Dp, nu, rho, pmax=4000):
    """leading light-cone (edge) ratio exact/LR: 2F1(1/2, h; h + 1/2; rho^2), h = (Dp + J')/2 (equal external dims)."""
    h = (Dp - 0.5 + 1j*nu)/2
    x = rho**2
    term = 1.0 + 0j; s = 1.0 + 0j
    for p in range(1, pmax):
        term *= (0.5 + p - 1)*(h + p - 1)/((h + 0.5 + p - 1)*p)*x
        s += term
        if abs(term) < 1e-18*abs(s):
            break
    return s


def c_coef(mu):
    """leading coefficient of P_mu(eta) ~ c_mu (2 eta)^mu:  Gamma(mu + 1/2)/(sqrt(pi) Gamma(mu + 1))."""
    return complex(mp.gamma(mu + mp.mpf(0.5))/(mp.sqrt(mp.pi)*mp.gamma(mu + 1)))
