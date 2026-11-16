"""
Stage 0: direct test of the heavy-block Regge-moment formula of arXiv:2106.10274 in d = 3 with exact blocks.

Computed (no large-Delta approximation anywhere):
   Pihat_{k,nu}[G^s_{Delta,J}] / (2 sin^2(pi (tau - 2 Dphi)/2))
       = (1/4) int_0^1 int_0^1 d rho d rhobar |rho-rhobar|^{d-2} (rho rhobar)^{(k-d-1)/2} P_{(2-d)/2+i nu}(eta) G^s_{Delta,J}
   [(C.1) with (2.39), (2.40); the change of variables (r, eta) -> (rho, rhobar) is worked out in REGGE_MOMENT_TTTT.md]
   G^s = (w wbar)^{-Dphi} g_{Delta,J}(w, wbar),  w = -(1-rho)^2/(4 rho) < 0   (rho = u-channel radial coordinate (2.41)),
   s-channel radial coordinate of the block: rho_s = -((1-sqrt(rho))/(1+sqrt(rho)))^2,
   heavy density (2.31):  Pihat[Delta,J] = (-1)^J/q_{Delta,J} * (the above).
Compared with (3.20) and its 1/m^2, 1/m^4 corrections (ancillary file of the paper):
   Pihat_{k,nu}[Delta,J] = gamma_{2Dphi+k-1}(nu)^2/a_{Dphi} * (2 P_{-1/2+i nu}(eta_AdS) + c_1/m^2 + c_2/m^4 + ...)/m^{2k}.

Usage: python3 s0_regge_moment.py [quick]
"""
import sys, time
import numpy as np, mpmath as mp
from blocks3d import Block3d
from chmrs_formulas import m2_eta, log_q, gamma_a, a_Dphi, Pnu
from chmrs_series import series

mp.mp.dps = 25

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

def grid(Delta, nt=(20, 28, 24), nu_nodes=(18, 30), T=(0.0, 6.0, 28.0, 85.0)):
    us, uw = [], []
    for (a, b, n) in ((0.0, 0.25, nu_nodes[0]), (0.25, 1.0, nu_nodes[1])):
        x, w = gl(a, b, n); us.append(x); uw.append(w)
    us, uw = np.concatenate(us), np.concatenate(uw)
    ts, tw = [], []
    for i in range(3):
        x, w = gl(T[i], T[i+1], nt[i]); ts.append(x); tw.append(w)
    ts, tw = np.concatenate(ts), np.concatenate(tw)
    return us, uw, ts, tw

class Moments:
    """All Regge moments Pihat_{k,nu} of one block, for arbitrary k, nu, Dphi, from one evaluation of the block
    on the (t, u) grid:  sqrt(rho) = s = t/(2 Delta (1+u)),  sqrt(rhobar) = s u."""
    def __init__(self, Delta, J, N, **gridkw):
        self.Delta, self.J = Delta, J
        us, uw, ts, tw = grid(Delta, **gridkw)
        U, Tt = np.meshgrid(us, ts, indexing='ij')
        S = Tt/(2*Delta*(1+U))
        ok = S < 0.999
        self.U, self.S, self.ok = U, S, ok
        self.W = np.outer(uw, tw)/(2*Delta*(1+U))        # du ds
        x, xb = S[ok], (S*U)[ok]
        rs = ((1-x)/(1+x))**2        # |rho_s| (smaller one, since x > xb)
        rb = ((1-xb)/(1+xb))**2
        blk = Block3d(Delta, J, N)
        h = blk.h_abs(rb, rs, tol=1e-13)
        # log of r_s^Delta * h  (the factor 4^Delta is added analytically)
        self.logg = np.full(S.shape, -np.inf)
        self.logg[ok] = Delta*0.5*(np.log(rb) + np.log(rs)) + np.log(h)
        self.us = us

    def log_moment(self, k, nu, Dphi, d=3):
        """(phase, log modulus) of  Pihat_{k,nu}[G]/(2 sin^2) / (-1)^J.  nu may be complex (needed for (C.16))."""
        U, S = self.U, self.S
        Pu = np.array([complex(Pnu(nu, (1+u*u)/(2*u))) for u in self.us])
        if abs(complex(nu).imag) < 1e-30:
            Pu = Pu.real
        rho, rhob = S**2, (S*U)**2
        # (1/2) * [d rho d rhobar = 4 s^3 u ds du] * (rho-rhobar) * (rho rhobar)^{(k-4)/2} * P * (16 rho rhobar)^Dphi/((1-rho)(1-rhob))^{2Dphi} * g
        with np.errstate(divide='ignore', invalid='ignore'):
            logw = (np.log(2.0) + 3*np.log(S) + np.log(U) + np.log(rho - rhob) + 0.5*(k-4)*np.log(rho*rhob)
                    + Dphi*np.log(16*rho*rhob) - 2*Dphi*np.log((1-rho)*(1-rhob)) + self.logg)
        logw[~self.ok] = -np.inf
        shift = np.max(logw)
        val = np.sum(np.exp(logw - shift)*Pu[:, None]*self.W)
        return val/abs(val), np.log(abs(val)) + shift + self.Delta*np.log(4.0)

def heavy_density(mom, k, nu, Dphi):
    sgn, lg = mom.log_moment(k, nu, Dphi)
    return mp.mpmathify(complex(sgn))*mp.exp(mp.mpf(lg) - log_q(mp.mpf(mom.Delta), mom.J, mp.mpf(Dphi)))

def C_heavy_density(mom, k, nu, Dphi, d=3):
    """Heavy density of the physical functional C_{k,nu} through the Regge-moment expansion (C.16):
       gamma^2/a C_{k,nu} = Pihat_{k,nu} + A Pihat_{k+2,nu} - B Pihat_{k+2,nu-2i} - conj(B) Pihat_{k+2,nu+2i} + O(Pi_{k+4})."""
    nu = mp.mpf(nu)
    A = (d-2)*(2*(k-3)*(k-1) + (k+1)*(d-4))/(4*(nu**2+1))
    B = (d+2j*nu)*(2*(k-3)*(k-1) + (k+1)*(d-2+2j*nu))/(8*nu*(nu-1j))
    tot = heavy_density(mom, k, nu, Dphi) + A*heavy_density(mom, k+2, nu, Dphi) - 2*(B*heavy_density(mom, k+2, nu-2j, Dphi)).real
    return (a_Dphi(Dphi)/gamma_a(2*Dphi+k-1, nu)**2*tot).real

def prediction(Delta, J, k, nu, Dphi, order, which='Pi'):
    m2, eta = m2_eta(mp.mpf(Delta), J)
    c = series(which, k, mp.mpf(nu), eta, mp.mpf(Dphi), 3, max_order=max(order, 0))
    tot = sum(c[n]/m2**n for n in range(order+1))
    if which == 'C':
        return (tot/m2**k).real
    return (gamma_a(2*Dphi+k-1, nu)**2/a_Dphi(Dphi)*tot/m2**k).real

def prediction_printed(Delta, J, k, nu, Dphi):
    """leading + first correction with (3.27) exactly as printed (Pochhammer (Dphi-d-1)_3)."""
    from chmrs_formulas import C_series_first
    m2, eta = m2_eta(mp.mpf(Delta), J)
    P0, P01 = C_series_first(k, mp.mpf(nu), mp.mpf(Delta), J, mp.mpf(Dphi), 3, printed=True)
    return (gamma_a(2*Dphi+k-1, nu)**2/a_Dphi(Dphi)*2*P01/m2**k).real

if __name__ == "__main__":
    quick = len(sys.argv) > 1 and sys.argv[1] == 'quick'
    out = []
    def pr(*a):
        line = " ".join(str(x) for x in a); print(line, flush=True); out.append(line)
    for Dphi in ((3.0,) if quick else (3.0, 2.0)):
        cases = [(20.3, 6), (40.3, 12)] if quick else [(20.3, 6), (30.3, 9), (40.3, 12), (60.3, 18), (80.3, 24), (120.3, 36), (40.3, 2), (40.3, 24), (60.3, 40), (80.3, 60)]
        pr(f"# d = 3, Dphi = {Dphi}.  R_n := (exact heavy density of Pihat_{{k,nu}}) / ((3.20) with corrections through 1/m^{{2n}} from the ancillary file)")
        pr("#   R_1' : the same with the 1/m^2 term of (3.27) exactly as printed.   S_2 := (exact heavy density of C_{k,nu} via (C.16)) / (seriesC through 1/m^4)")
        pr("#  Delta    J      m     beta   k   nu       exact Pihat[Delta,J]       R_0          R_1          R_2      m^6(R_2-1)      R_1'          S_2")
        for (Delta, J) in cases:
            t0 = time.time()
            N = int(min(8000, max(800, 45*Delta)))//2*2
            mom = Moments(Delta, J, N)
            m2, eta = m2_eta(mp.mpf(Delta), J)
            for k in (2, 4):
                for nu in (0.0, 0.7, 2.0, 4.0):
                    ex = heavy_density(mom, k, nu, Dphi).real
                    R = [ex/prediction(Delta, J, k, nu, Dphi, n) for n in (0, 1, 2)]
                    Rp = ex/prediction_printed(Delta, J, k, nu, Dphi)
                    S2 = C_heavy_density(mom, k, nu, Dphi)/prediction(Delta, J, k, nu, Dphi, 2, 'C') if nu > 0 else mp.nan
                    pr(f"  {Delta:6.1f} {J:4d} {float(mp.sqrt(m2)):7.3f} {float(mp.acosh(eta)):6.3f}  {k}  {nu:4.1f}  {mp.nstr(ex, 12):>22s}  "
                       + "  ".join(f"{float(r):11.8f}" for r in R) + f"  {float((R[2]-1)*m2**3):11.1f}  {float(Rp):11.8f}  {float(S2):11.8f}")
            pr(f"#   (N = {N}, {time.time()-t0:.0f} s)")
    if not quick:
        open('out_s0_regge_moment.txt', 'w').write("\n".join(out) + "\n")
