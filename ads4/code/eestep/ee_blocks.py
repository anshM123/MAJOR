r"""
ee_blocks (EESTEP_CLOSE.md §1.2): d = 3 continuous-spin Casimir eigenfunctions with UNEQUAL external weights on the (24) side,
(a, b) = (0, b), by a radial recursion that generalises lrc/lrc_blocks.py (read, not imported) to b != 0.

Normalisation (as lrc_blocks):  Gh_b(r, eta; nu) = sum_{n >= 0} r^n sum_{m = -n..n, m = n mod 2} B[n, m](nu) P_{J'+m}(eta),  B[0,0] = 1,
J' = l = -1/2 + i nu, dimension Dp (= 1 here).  Block = (4r)^Dp Gh_b.  Gh_b -> P_l(eta) as r -> 0 at fixed eta, for every b.

Casimir.  Dolan-Osborn with a = 0:  D_z = z^2(1-z) d_z^2 - (b+1) z^2 d_z.  Relative to b = 0 the (1/2)-Casimir acquires
    -b (z^2 d_z + zb^2 d_zb) = -(4 b r/den) [eta (1-r^2) r d_r + (1+r^2)(eta^2-1) d_eta],   den = (1+r^2)^2 - 4 r^2 eta^2
(checked symbolically, chk_op in the note).  Multiplying by (1-r^2) den as in blocks3d.py, on r^E P_j it adds
    r^1: -4b [E eta P_j - (1-eta^2) P_j'],   r^3: +8b E eta P_j,   r^5: -4b [E eta P_j + (1-eta^2) P_j']
to the b = 0 operator (s^0..s^3 of blocks3d.py).  eta P_j and (1-eta^2) P_j' are re-expanded in P_{j+-1} (exact for complex j).
Odd levels n now appear; the denominators K(Dp+n, l+m)/2 are b-independent.
"""
import numpy as np
import mpmath as mp


def legendre_coeffs2(j):
    """eta^2 P_j = e_up P_{j+2} + e_mid P_j + e_dn P_{j-2};  eta(1-eta^2)P_j' = d_up P_{j+2} + d_mid P_j + d_dn P_{j-2}."""
    j = np.asarray(j, dtype=complex)
    e_up = (j + 1)*(j + 2)/((2*j + 1)*(2*j + 3))
    e_mid = (j + 1)**2/((2*j + 1)*(2*j + 3)) + j**2/((2*j + 1)*(2*j - 1))
    e_dn = j*(j - 1)/((2*j + 1)*(2*j - 1))
    c = j*(j + 1)/(2*j + 1)
    d_up = -c*(j + 2)/(2*j + 3)
    d_mid = c*j/(2*j - 1) - c*(j + 1)/(2*j + 3)
    d_dn = c*(j - 1)/(2*j - 1)
    return e_up, e_mid, e_dn, d_up, d_mid, d_dn


def legendre_coeffs1(j):
    """eta P_j = f_up P_{j+1} + f_dn P_{j-1};   (1-eta^2) P_j' = g_up P_{j+1} + g_dn P_{j-1}."""
    j = np.asarray(j, dtype=complex)
    f_up = (j + 1)/(2*j + 1); f_dn = j/(2*j + 1)
    c = j*(j + 1)/(2*j + 1)
    return f_up, f_dn, -c, c


def block_coeffs_b(nu, b, N, Dp=1.0):
    """B[n, i], n = 0..N, m = -N + i (i = 0..2N).  Complex spin l = -1/2 + i nu (nu may be complex), dimension Dp, (a, b) = (0, b)."""
    Jp = -0.5 + 1j*nu
    M = 2*N + 1
    m = -N + np.arange(M)
    jj = Jp + m
    C = Dp*(Dp - 3) + Jp*(Jp + 1)
    e_up, e_mid, e_dn, d_up, d_mid, d_dn = legendre_coeffs2(jj)
    f_up, f_dn, g_up, g_dn = legendre_coeffs1(jj)
    B = np.zeros((N + 1, M), dtype=complex)
    i0 = N
    B[0, i0] = 1.0

    def apply2(vec, a2, aD, a1):
        out = (a2*e_mid + aD*d_mid + a1)*vec
        t = (a2*e_up + aD*d_up)*vec
        out[2:] += t[:-2]
        t = (a2*e_dn + aD*d_dn)*vec
        out[:-2] += t[2:]
        return out

    def apply1(vec, aeta, aD):
        """(aeta * eta + aD * (1-eta^2) d_eta) sum vec_j P_j  ->  coefficients (shift +-1)."""
        out = np.zeros_like(vec)
        t = (aeta*f_up + aD*g_up)*vec
        out[1:] += t[:-1]
        t = (aeta*f_dn + aD*g_dn)*vec
        out[:-1] += t[1:]
        return out

    for n in range(1, N + 1):
        rhs = np.zeros(M, dtype=complex)
        if n >= 2:
            E = Dp + n - 2
            K = E*(E - 3) + jj*(jj + 1) - C
            rhs += apply2(B[n - 2], -4*(K/2 + E), 4.0, (K/2 + E))
        if n >= 4:
            E = Dp + n - 4
            K = E*(E - 3) + jj*(jj + 1) - C
            rhs += apply2(B[n - 4], 4*(K/2 + 2*E), -4.0, -(K/2 + 2*E))
        if n >= 6:
            E = Dp + n - 6
            K = E*(E - 3) + jj*(jj + 1) - C
            rhs += -(K/2 + 3*E)*B[n - 6]
        if b != 0:
            E = Dp + n - 1
            rhs += apply1(B[n - 1], -4*b*E, 4*b)
            if n >= 3:
                E = Dp + n - 3
                rhs += apply1(B[n - 3], 8*b*E, 0.0)
            if n >= 5:
                E = Dp + n - 5
                rhs += apply1(B[n - 5], -4*b*E, -4*b)
        E = Dp + n
        K = E*(E - 3) + jj*(jj + 1) - C
        lo, hi = i0 - n, i0 + n
        sl = slice(lo, hi + 1, 2)                    # only m = n (mod 2); the other parity is identically zero
        B[n, sl] = -rhs[sl]/(K[sl]/2)
    return B, m


def legendre_scaled(nu, beta, kmax, dps=30):
    """Pt[k] = P_{J'+k}(cosh beta) e^{-k beta}, k = 0..kmax, J' = -1/2 + i nu (real nu): stable forward recurrence (as lrc_blocks)."""
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
    """Y[n] = sum_m B[n, m] e^{-(n - |m|) beta} Pt_{|m|} (conj for m < 0; real nu), so that Gh_b(rho, beta) = sum_n Y[n] rho^n."""
    N1, M = B.shape
    am = np.abs(m)
    P = np.where(m >= 0, Pt[am], np.conj(Pt[am]))
    n = np.arange(N1)[:, None]
    with np.errstate(over='ignore', invalid='ignore'):
        fac = np.where(am[None, :] <= n, np.exp(-(n - am[None, :])*beta), 0.0)
    return (B*fac*P[None, :]).sum(axis=1)


def Gh(Y, rho):
    rho = np.atleast_1d(np.asarray(rho, dtype=complex))
    pw = rho[:, None]**np.arange(len(Y))[None, :]
    return pw @ Y


def v_cross(rho, beta):
    """v = (1-z)(1-zb) = [(1-rho)(1-rhob)/((1+rho)(1+rhob))]^2, rhob = rho e^{-2 beta}  (= (x12^2/x14^2)^2 in the H7 frame)."""
    rb = rho*np.exp(-2*beta)
    return ((1 - rho)*(1 - rb)/((1 + rho)*(1 + rb)))**2


def Gh_point(nu, b, rho, beta, N=200, Dp=1.0):
    """Gh_b at one point.  b < 0 via the Euler-type identity G^(0,b) = v^{-b} G^(0,-b) (exact; avoids the cancellations of the
    negative-b series, which lose up to 9 digits at rho = 0.85, b = -4)."""
    bb = abs(b)
    B, m = block_coeffs_b(nu, bb, N, Dp)
    Pt = legendre_scaled(nu, beta, N + 1)
    g = Gh(radial_Y(B, m, Pt, beta), rho)[0]
    return g*v_cross(rho, beta)**bb if b < 0 else g


def edge_Fb(nu, b, rho, Dp=1.0, dps=30):
    """Top-shell (light-cone edge) factor for (a, b) = (0, b):  F_b = (z/(4 rho))^h 2F1(h, h+b; 2h; z),  z = 4 rho/(1+rho)^2,
    h = (Dp + l)/2.  b = 0: = 2F1(1/2, h; h+1/2; rho^2) (lrc_blocks.edge_F)."""
    with mp.workdps(dps):
        h = (mp.mpf(Dp) - mp.mpf(0.5) + 1j*mp.mpmathify(nu))/2
        rr = mp.mpf(rho)
        z = 4*rr/(1 + rr)**2
        return complex((z/(4*rr))**h*mp.hyp2f1(h, h + b, 2*h, z))


def c_coef(mu):
    """P_mu(eta) ~ c_mu (2 eta)^mu:  c_mu = Gamma(mu + 1/2)/(sqrt(pi) Gamma(mu + 1))."""
    return complex(mp.gamma(mu + mp.mpf(0.5))/(mp.sqrt(mp.pi)*mp.gamma(mu + 1)))
