"""
lsl/lsl_smear.py (LARGESPIN_L.md §2.3, §5): (a) the exact positive plane-wave (Bernstein) smearing of a heavy state and its size;
(b) the scalar subleading Regge term (CHMRS 2106.10274 eq. (3.20)) along fixed twist up to J = 10^6 and its J -> infinity limit.

(a) Scalar heavy factor in the gauge x5 = 0, x1 = -e, x2 = inf (CHMRS eq. "neat B integral", d = 3):
        f(w) = (-n.w)^{-(J+1)} (-w^2)^{-lam},  lam = tau/2 + Dphi - 2,   w = e + x + y,
    = int dmu(p) e^{p.w}, dmu = [s^J e^{-s}/J!] ds x [(2lam-1)/(2 pi Gamma(2lam))] (-q^2)^{lam-3/2} d^3q,  p = s n + q   (exact, positive).
    The leading Laplace transform of the Regge kernel at momentum p is gamma^2 P_{-1/2+i nu}(cosh beta(p)) |p|^{-2(k+2Dphi-1)}
    (CHMRS eq. Fourier identity), beta(p) = 2 artanh(|vec p|/p^0).  "Smearing-only" image:
        S(nu) = E_mu[|p|^{-2a'} e^{-p0}... P(cosh beta(p))] / E_mu[|p|^{-2a'}]   (e^{p.e} = e^{-s-q0} is in dmu by sampling)
    sampled by Monte Carlo (s ~ Gamma(J+1); q: cosh(chi) ~ Pareto(2lam-1), |q| | chi ~ Gamma(2lam, cosh chi); azimuth uniform).
    Compared with the full CHMRS (3.20) relative correction at the same state.
(b) CHMRS (3.20), d = 3, a = 2Dphi + k - 3/2:
    (m^{2k}/2) C = P + (1/6m^2){ P'(eta)[(eta^2-1)(nu^2 + c1) + 1] + eta P [(nu^2 - a(a+1)) c2 + c3] },
    c1 = -3a^2 + 6a(2Dphi-3) - 2, c2 = 3a - 6Dphi + 19/2, c3 = -2(Dphi-4)_3 - 2(a+1/2)(1 - a^2 - a/2).
    Large-J limit at fixed tau (eta/m^2 -> 1/(2(tau-2)^2), per travelling wave d/dbeta -> (+- i nu - 1/2)):
        F_+-(nu, tau) - 1 = [ (+- i nu - 1/2)(nu^2 + c1) + (nu^2 - a(a+1)) c2 + c3 ] / (12 (tau-2)^2).
Output: cache/smear.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
from scipy.special import gammaln
import lsl_core as C

OUT = os.path.join(C.CACHE, 'smear.txt')
lines = []
rng = np.random.default_rng(7)


def pr(s=''):
    print(s, flush=True); lines.append(s)


DPHI, K = 3.0, 2
A = 2*DPHI + K - 1.5
c1 = -3*A*A + 6*A*(2*DPHI - 3) - 2
c2 = 3*A - 6*DPHI + 9.5
c3 = -2*(DPHI - 4)*(DPHI - 3)*(DPHI - 2) - 2*(A + 0.5)*(1 - A*A - A/2)


def chmrs_rel(J, tau, nu):
    """relative O(1/m^2) correction of CHMRS (3.20) divided by the local LR amplitude, and the bare ratio at nu where P != 0."""
    Delta = tau + J
    m2 = (Delta - J - 2)*(Delta + J - 1)
    eta = 1 + (2*J + 1)**2/(2*m2)
    beta = np.arccosh(eta)
    P, dP = C.P_real(beta, nu, nphi=4000)
    Pp = dP/np.sinh(beta)
    corr = (Pp*((eta*eta - 1)*(nu*nu + c1) + 1) + eta*P*((nu*nu - A*(A + 1))*c2 + c3))/(6*m2)
    amp = np.sqrt(P**2 + (dP + 0.5*P)**2/(nu**2 + 0.25))
    return corr, P, amp, beta, m2


def smear(J, tau, nu, N=400000):
    lam = tau/2 + DPHI - 2
    s = rng.gamma(J + 1, 1.0, N)
    u = rng.random(N)**(-1/(2*lam - 1))                  # cosh chi
    qm = rng.gamma(2*lam, 1.0, N)/u
    chi = np.arccosh(u); ph = rng.random(N)*2*np.pi
    q0 = qm*np.cosh(chi); qx = qm*np.sinh(chi)*np.cos(ph); qy = qm*np.sinh(chi)*np.sin(ph)
    p0 = s + q0; px = s + qx; py = qy                    # n = (1, 1, 0)
    pv = np.sqrt(px*px + py*py); M = np.sqrt(p0*p0 - pv*pv)
    bp = 2*np.arctanh(pv/p0)
    lw = -2*(K + 2*DPHI - 1)*np.log(M); w = np.exp(lw - lw.max())
    Delta = tau + J; m2 = (Delta - J - 2)*(Delta + J - 1); bcl = np.arccosh(1 + (2*J + 1)**2/(2*m2))
    wb = w/w.sum()
    mb = (wb*bp).sum(); sb = np.sqrt((wb*(bp - mb)**2).sum())
    mM = (wb*M).sum()/np.sqrt(m2)
    # smeared image on a nu grid, via sorting beta(p) into fine bins (P evaluated on bin centres)
    hb, edges = np.histogram(bp, bins=400, range=(mb - 8*sb, mb + 8*sb), weights=wb)
    bc = 0.5*(edges[1:] + edges[:-1])
    Pb = np.array([C.P_real(b, nu, nphi=2000)[0] for b in bc])
    S = (hb[:, None]*Pb).sum(0)/hb.sum()
    return dict(mean_beta=mb, sd_beta=sb, beta_cl=bcl, M_over_m=mM, S=S)


def main():
    pr(__doc__.split('Output')[0])
    pr(f"(constants) a = {A}, c1 = {c1}, c2 = {c2}, c3 = {c3}")
    nus = np.array([0.0, 2.0, 5.0, 10.0, 20.0, 30.0])
    pr("\n(a) Bernstein smearing of the heavy factor (exact, positive), tau = 1000:")
    for J in (50, 100, 344, 3000):
        r = smear(J, 1000.0, nus)
        corr, P, amp, beta, m2 = chmrs_rel(J, 1000.0, nus)
        Pcl = C.P_real(r['beta_cl'], nus, nphi=4000)[0]
        pr(f"  J = {J}: beta_cl = {r['beta_cl']:.4f}; <beta(p)> - beta_cl = {r['mean_beta'] - r['beta_cl']:+.4f}, sd(beta(p)) = {r['sd_beta']:.4f} "
           f"(sqrt(2/tau)*... = {np.sqrt(2/1000):.4f}); <|p|>/m = {r['M_over_m']:.4f}")
        pr("     nu:                         " + "  ".join(f"{v:9.1f}" for v in nus))
        pr("     smearing only (S - P)/amp:  " + "  ".join(f"{v:+9.2e}" for v in (r['S'] - Pcl)/amp))
        pr("     CHMRS (3.20) full corr/amp: " + "  ".join(f"{v:+9.2e}" for v in corr/amp))
    pr("\n(b) CHMRS (3.20) correction / local LR amplitude along fixed twist tau = 1000 (J -> inf saturates; parameter eta/m^2 <= 1.5/(tau-2)^2)")
    nus2 = np.array([0.0, 0.5, 1.0, 2.0, 5.0, 10.0])
    pr("     J        beta    eta/m^2*(tau-2)^2   " + "  ".join(f"nu={v:<6}" for v in nus2))
    for J in (26, 50, 100, 344, 1000, 3000, 10000, 100000, 1000000):
        corr, P, amp, beta, m2 = chmrs_rel(J, 1000.0, nus2)
        eta = 1 + (2*J + 1)**2/(2*m2)
        pr(f"  {J:8d}  {beta:7.4f}   {eta/m2*998**2:8.4f}          " + "  ".join(f"{v:+.2e}" for v in corr/amp))
    tp = 998.0
    Fp = lambda nu: ((1j*nu - 0.5)*(nu*nu + c1) + (nu*nu - A*(A + 1))*c2 + c3)/(12*tp*tp)
    pr("  J = inf: |F_+(nu) - 1| = " + "  ".join(f"{abs(Fp(v)):.2e}" for v in nus2) + "   (lightcone amplitude, eq. above)")
    with open(OUT, 'w') as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
