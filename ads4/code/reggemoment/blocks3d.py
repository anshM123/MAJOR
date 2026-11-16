"""
Scalar conformal blocks in d = 3 for equal external dimensions, by the Casimir recursion in radial coordinates.

Conventions.  z = 4 rho/(1+rho)^2, zbar likewise;  r = sqrt(rho rhobar), eta = (rho+rhobar)/(2 r).
  g_{Delta,J}(z,zbar) = (4 r)^Delta sum_{n = 0,2,4,...} r^n sum_j B_{n,j} P_j(eta),      B_{0,J} = J!/(1/2)_J,
so that g ~ z^{(Delta-J)/2} zbar^{(Delta+J)/2} for z << zbar << 1 (the normalisation in which the mean-field OPE
coefficients are CHMRS (2.33); checked in chmrs_formulas.py / s0_block_checks.py).

Derivation of the recursion (done by hand, verified numerically in s0_block_checks.py against the Casimir PDE in
(z, zbar)).  With D_z = z^2(1-z) d_z^2 - z^2 d_z, the Casimir equation is
  [2 D_z + 2 D_zbar + 2 (d-2) z zbar/(z-zbar) ((1-z) d_z - (1-zbar) d_zbar)] g = [Delta(Delta-d) + J(J+d-2)] g.
In rho: D_z = rho^2 d_rho^2 - 2 rho^3/(1-rho^2) d_rho, and with th = r d_r, L = -(1-eta^2) d_eta^2 + (d-1) eta d_eta,
  (1/2) Casimir = (1/2)[th^2 - d th + L] - (d-2) r^2/(1-r^2) th - (2 r^2 (2 eta^2-1) - 2 r^4)/den th
                  + 4 r^2 eta (1-eta^2)/den d_eta,          den = (1+r^2)^2 - 4 r^2 eta^2.
Multiplying by (1-r^2) den, with s = r^2, K(E,j) = E(E-3) + j(j+1) - C (d = 3), acting on r^E P_j:
  s^0 : K/2
  s^1 : (K/2 + E)(1 - 4 eta^2) + 4 eta (1-eta^2) d_eta
  s^2 : (K/2 + 2E)(4 eta^2 - 1) - 4 eta (1-eta^2) d_eta
  s^3 : -(K/2 + 3E)
and  eta^2 P_j,  eta (1-eta^2) P_j'  are re-expanded in P_{j-2}, P_j, P_{j+2}.
"""
import numpy as np
import mpmath as mp


def _eta2_coeffs(j):
    """eta^2 P_j = up*P_{j+2} + mid*P_j + dn*P_{j-2}."""
    j = np.asarray(j, dtype=float)
    up = (j+1)*(j+2)/((2*j+1)*(2*j+3))
    mid = (j+1)**2/((2*j+1)*(2*j+3)) + np.where(j > 0, j**2/((2*j+1)*(2*j-1)), 0.0)
    dn = np.where(j > 1, j*(j-1)/((2*j+1)*(2*j-1)), 0.0)
    return up, mid, dn


def _etaD_coeffs(j):
    """eta (1-eta^2) P_j' = up*P_{j+2} + mid*P_j + dn*P_{j-2}."""
    j = np.asarray(j, dtype=float)
    c = j*(j+1)/(2*j+1)
    up = -c*(j+2)/(2*j+3)
    mid = np.where(j > 0, c*j/(2*j-1), 0.0) - c*(j+1)/(2*j+3)
    dn = np.where(j > 1, c*(j-1)/(2*j-1), 0.0)
    return up, mid, dn


class Block3d:
    """Radial coefficients B_{n,j}, n = 0, 2, ..., N.  Stored as a dense array Bt[n/2, i], j = jpar + 2 i."""

    def __init__(self, Delta, J, N):
        self.Delta, self.J, self.N = float(Delta), int(J), int(N)
        Delta, J = self.Delta, self.J
        self.jpar = J % 2
        nlev = N//2 + 1
        self.jmax = J + N
        nj = (self.jmax - self.jpar)//2 + 1
        jj = self.jpar + 2*np.arange(nj)
        self.jj = jj
        C = Delta*(Delta-3) + J*(J+1)
        e_up, e_mid, e_dn = _eta2_coeffs(jj)
        d_up, d_mid, d_dn = _etaD_coeffs(jj)
        B = np.zeros((nlev, nj))
        i0 = (J - self.jpar)//2
        B[0, i0] = float(mp.factorial(J)/mp.rf(mp.mpf(1)/2, J))

        def apply(vec, a_eta2, a_etaD, a_one):
            """coefficients of P_j in  (a_eta2*eta^2 + a_one + a_etaD*eta(1-eta^2)d_eta) sum_j vec_j P_j ;
            a_* are arrays over the source index j."""
            out = (a_eta2*e_mid + a_etaD*d_mid + a_one)*vec
            t = (a_eta2*e_up + a_etaD*d_up)*vec          # goes to j+2
            out[1:] += t[:-1]
            t = (a_eta2*e_dn + a_etaD*d_dn)*vec          # goes to j-2
            out[:-1] += t[1:]
            return out

        for lev in range(1, nlev):
            n = 2*lev
            rhs = np.zeros(nj)
            # s^1 acting on level n-2
            E = Delta + n - 2
            K = E*(E-3) + jj*(jj+1) - C
            rhs += apply(B[lev-1], -4*(K/2 + E), 4.0*np.ones(nj), (K/2 + E))
            if lev >= 2:
                E = Delta + n - 4
                K = E*(E-3) + jj*(jj+1) - C
                rhs += apply(B[lev-2], 4*(K/2 + 2*E), -4.0*np.ones(nj), -(K/2 + 2*E))
            if lev >= 3:
                E = Delta + n - 6
                K = E*(E-3) + jj*(jj+1) - C
                rhs += -(K/2 + 3*E)*B[lev-3]
            E = Delta + n
            K = E*(E-3) + jj*(jj+1) - C
            lo = max(0, (J - n - self.jpar + 1)//2)
            hi = (J + n - self.jpar)//2
            B[lev, lo:hi+1] = -rhs[lo:hi+1]/(K[lo:hi+1]/2)
        self.B = B

    # ------------------------------------------------------------------
    def h_abs(self, rho_big, rho_small, tol=1e-14, nmax=None):
        """sum_{n,j} B_{n,j} r^n P_j(|eta|)  for rho, rhobar real of the same sign with moduli rho_big >= rho_small
        (so |eta| = cosh chi >= 1).  For rho, rhobar < 0 the block is (-1)^J (4 r)^Delta h_abs;  for rho, rhobar > 0 it is
        (4 r)^Delta h_abs.  Vectorised over points."""
        rb = np.atleast_1d(np.asarray(rho_big, dtype=float))
        rs = np.atleast_1d(np.asarray(rho_small, dtype=float))
        npt = rb.size
        N = self.N if nmax is None else min(self.N, nmax)
        J = self.J
        lb, ls = np.log(rb), np.log(rs)
        # scaled Legendre: Pt_j = P_j(cosh chi) e^{-j chi},  e^{chi} = sqrt(rb/rs)
        emchi = np.sqrt(rs/rb)
        eta = 0.5*(1/emchi + emchi)
        jmaxN = J + N
        Pt = np.empty((npt, jmaxN + 1))
        Pt[:, 0] = 1.0
        if jmaxN >= 1:
            Pt[:, 1] = eta*emchi
        for j in range(1, jmaxN):
            Pt[:, j+1] = ((2*j+1)*eta*emchi*Pt[:, j] - j*emchi**2*Pt[:, j-1])/(j+1)
        out = np.zeros(npt)
        active = np.ones(npt, dtype=bool)
        for lev in range(0, N//2 + 1):
            n = 2*lev
            lo = max(0, (J - n - self.jpar + 1)//2)
            hi = (J + n - self.jpar)//2
            j = self.jj[lo:hi+1]
            idx = np.nonzero(active)[0]
            if idx.size == 0:
                break
            # r^n e^{j chi} = rb^{(n+j)/2} rs^{(n-j)/2}
            W = np.exp(np.outer(lb[idx], (n + j)/2.0) + np.outer(ls[idx], (n - j)/2.0))
            term = (W*Pt[idx][:, j]) @ self.B[lev, lo:hi+1]
            out[idx] += term
            if lev > 4:
                small = np.abs(term) < tol*np.abs(out[idx])
                # geometric tail: stop a point only when the term is far below tolerance
                active[idx[small & (rb[idx]**2 < 0.999)]] = False
        return out

    def g(self, rho, rhobar, **kw):
        """Block value for real rho, rhobar of equal sign, |rho|,|rhobar| < 1."""
        rho = np.atleast_1d(np.asarray(rho, dtype=float)); rhobar = np.atleast_1d(np.asarray(rhobar, dtype=float))
        a, b = np.abs(rho), np.abs(rhobar)
        big, small = np.maximum(a, b), np.minimum(a, b)
        sgn = np.where(rho < 0, (-1.0)**self.J, 1.0)
        r = np.sqrt(big*small)
        return sgn*np.exp(self.Delta*np.log(4*r))*self.h_abs(big, small, **kw)


def block_mp(Delta, J, N, rho, rhobar):
    """Slow mpmath reference implementation of the same recursion (independent code path: explicit loops, exact
    Legendre algebra), used only for validation.  Returns g(rho, rhobar)."""
    Delta = mp.mpf(Delta)
    C = Delta*(Delta-3) + J*(J+1)
    def e2(j):
        up = mp.mpf((j+1)*(j+2))/((2*j+1)*(2*j+3))
        mid = mp.mpf((j+1)**2)/((2*j+1)*(2*j+3)) + (mp.mpf(j*j)/((2*j+1)*(2*j-1)) if j > 0 else 0)
        dn = mp.mpf(j*(j-1))/((2*j+1)*(2*j-1)) if j > 1 else mp.mpf(0)
        return up, mid, dn
    def eD(j):
        c = mp.mpf(j*(j+1))/(2*j+1)
        up = -c*(j+2)/(2*j+3)
        mid = (c*j/(2*j-1) if j > 0 else 0) - c*(j+1)/(2*j+3)
        dn = c*(j-1)/(2*j-1) if j > 1 else mp.mpf(0)
        return up, mid, dn
    B = {0: {J: mp.factorial(J)/mp.rf(mp.mpf(1)/2, J)}}
    def K(E, j):
        return E*(E-3) + j*(j+1) - C
    for n in range(2, N+1, 2):
        rhs = {}
        def add(j, v):
            if j >= 0:
                rhs[j] = rhs.get(j, 0) + v
        for (shift, a2f, aDf, a1f) in ((2, lambda Kh, E: -4*(Kh+E), 4, lambda Kh, E: Kh+E),
                                       (4, lambda Kh, E: 4*(Kh+2*E), -4, lambda Kh, E: -(Kh+2*E)),
                                       (6, lambda Kh, E: 0, 0, lambda Kh, E: -(Kh+3*E))):
            if n - shift < 0:
                continue
            E = Delta + n - shift
            for j, v in B[n-shift].items():
                Kh = K(E, j)/2
                u2, m2_, d2 = e2(j); uD, mD, dD = eD(j)
                a2, a1 = a2f(Kh, E), a1f(Kh, E)
                add(j+2, (a2*u2 + aDf*uD)*v); add(j, (a2*m2_ + aDf*mD + a1)*v); add(j-2, (a2*d2 + aDf*dD)*v)
        B[n] = {j: -v/(K(Delta+n, j)/2) for j, v in rhs.items() if abs(j-J) <= n}
    rho, rhobar = mp.mpf(rho), mp.mpf(rhobar)
    r = mp.sqrt(rho*rhobar); eta = (rho+rhobar)/(2*r)
    s = mp.mpf(0)
    for n, lev in B.items():
        s += r**n*sum(v*mp.legendre(j, eta) for j, v in lev.items())
    return (4*r)**Delta*s
