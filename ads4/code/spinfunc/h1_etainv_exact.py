"""
h1 (M3b): EXACT eta-inversion of the route-A' weight psi_2 (Mehler-Dirichlet / Abel representation), all nu including small nu.

Problem.  The D-safe weight psi_2(p) on the pure-harmonic C_2 defines the Regge-side eta-kernel  K(eta) = int_0^inf dnu phi(nu) P_{-1/2+i nu}(eta),
phi(nu) = psi_2(nu/Lam) on [0, Lam].  The CLS-D object C^CLS_2 has leading kernel eta P_nu (j = 1).  We need phit with
   int_0^inf dnu phit(nu) eta P_nu(eta) = K(eta)   for all eta >= 1.
Solution [D].  Mehler-Dirichlet: P_{-1/2+i nu}(cosh xi) = (sqrt2/pi) int_0^xi cos(nu t) dt / sqrt(cosh xi - cosh t), so K = A[Phi] with
Phi(t) = int_0^inf phi(nu) cos(nu t) dnu and A the Abel-type operator.  For any Phit with Ht' = Phit, Ht(0) = 0, integration by parts gives
   cosh(xi) A[Phit](xi) = A[cosh(t) Phit + (1/2) sinh(t) Ht](xi)   (exact).
So eta K~ = K  <=>  Phi = cosh Phit + (1/2) sinh Ht, a first-order linear ODE for Ht with solution
   Ht(t) = sech^{1/2}(t) int_0^t Phi(s) sech^{1/2}(s) ds,     Phit = Phi/cosh - (1/2) tanh(t) Ht,
and phit(nu) = (2/pi) int_0^inf Phit(t) cos(nu t) dt.  No analyticity of phi is needed (psi_2 ~ p at p -> 0 is not even in nu).
Numerics: Phi(t) = Lam Phi1(Lam t), Phi1(tau) = int_0^1 psi_2(p) cos(p tau) dp in closed form (terminating integration by parts, psi_2 a
polynomial) for tau >= 400, Gauss-Legendre below; int_0^t Phi w = Psi w - int Psi w' with Psi(t) = int_0^1 (psi_2/p) sin(p Lam t) dp (closed form),
cumulative Simpson; t-grid dt = pi/(16 Lam), t <= 60 (Phit ~ e^{-t/2}).
Checks: (V1) large-nu symbol  phit - phi = phi''/2 - phi'/(2nu) + phi/(2nu^2) + O(nu^-4) in the interior; (V2) Abel identity cosh(xi) K~ = K at
sample xi (direct quadrature); (V3) graviton-pole identity  int phit g~ = int phi g, g = 1/(nu^2 + 9/4), g~(nu) = Re[(1 + 1/(2 i nu)) g(nu - i)]
(Legendre recursion; the statement that the inversion leaves the graviton-pole light term, hence A and B_AdS, unchanged).
Outputs: dk0 = -2 int (psit_2 - psi_2) p^6 dp (flat |g3|^2 coefficient, NF normalised to k(0) = 1), tail mass beyond p = 1, phit on grids
(pickled for h2).  Usage: python3 h1_etainv_exact.py [Dg ...] -> cache/h1_etainv_exact.txt, cache/h1_<Dg>.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from numpy.polynomial import polynomial as Pl
from scipy.integrate import cumulative_simpson
from spf_common import NearGap
from common import gl01

out = []
def pr(s=""):
    print(s, flush=True); out.append(s)

def gl(a, b, n):
    x, w = np.polynomial.legendre.leggauss(n); return 0.5*(b - a)*x + 0.5*(b + a), 0.5*(b - a)*w

def Cint(h, T):
    """int_0^1 h(p) cos(T p) dp for polynomial h (monomial coeffs), T >= ~400 (terminating closed form)."""
    if len(h) == 0 or not np.any(h):
        return np.zeros_like(T)
    h1 = Pl.polyder(h, 1) if len(h) > 1 else np.array([0.0]); h2 = Pl.polyder(h, 2) if len(h) > 2 else np.array([])
    val = Pl.polyval(1.0, h)*np.sin(T)/T + (Pl.polyval(1.0, h1)*np.cos(T) - Pl.polyval(0.0, h1))/T**2
    if len(h2):
        val = val - Cint(h2, T)/T**2
    return val

def Sint(g, T):
    """int_0^1 g(p) sin(T p) dp, T >= ~400."""
    g1 = Pl.polyder(g, 1) if len(g) > 1 else np.array([])
    val = (Pl.polyval(0.0, g) - Pl.polyval(1.0, g)*np.cos(T))/T
    if len(g1):
        val = val + Cint(g1, T)/T
    return val

class Inverter:
    def __init__(self, psi_coeffs, Lam, tmax=60.0, ref=16):
        self.c = np.asarray(psi_coeffs, float); self.Lam = Lam
        assert abs(self.c[0]) < 1e-14                 # psi_2(0) = 0
        self.g = self.c[1:]                           # psi_2/p
        dt = np.pi/(ref*Lam); t = np.arange(0.0, tmax + dt, dt); self.t, self.dt = t, dt
        T = Lam*t; tau0 = 400.0
        pq, wq = gl(0.0, 1.0, 600); f = Pl.polyval(pq, self.c); gq = Pl.polyval(pq, self.g)
        Phi1 = np.empty_like(T); Psi = np.empty_like(T)
        lo = T < tau0
        Phi1[lo] = np.cos(np.outer(T[lo], pq)) @ (wq*f)
        Psi[lo] = np.sin(np.outer(T[lo], pq)) @ (wq*gq)
        Phi1[~lo] = Cint(self.c, T[~lo]); Psi[~lo] = Sint(self.g, T[~lo])
        Phi = Lam*Phi1                                  # int_0^Lam phi(nu) cos(nu t) dnu,  phi(nu) = psi_2(nu/Lam)
        w = 1/np.sqrt(np.cosh(t)); wp = -0.5*w*np.tanh(t)
        I = Psi*w - np.concatenate([[0.0], cumulative_simpson(Psi*wp, x=t)])
        Ht = w*I
        self.Phi, self.Phit = Phi, Phi/np.cosh(t) - 0.5*np.tanh(t)*Ht
        self.Psi = Psi
    def phit(self, nu, chunk=None):
        """phit(nu) = (2/pi) int_0^inf Phit cos(nu t) dt (trapezoid, even integrand)."""
        nu = np.atleast_1d(np.asarray(nu, float)); res = np.empty_like(nu)
        chunk = chunk or max(1, int(4e6/len(self.t)))
        wt = np.full_like(self.t, self.dt); wt[0] *= 0.5
        for i in range(0, len(nu), chunk):
            res[i:i + chunk] = (2/np.pi)*(np.cos(np.outer(nu[i:i + chunk], self.t)) @ (wt*self.Phit))
        return res
    def phi(self, nu):
        p = np.asarray(nu, float)/self.Lam
        return np.where((p >= 0) & (p <= 1), Pl.polyval(p, self.c), 0.0)
    def abel(self, Phi, xi, n=4000):
        """(sqrt2/pi) int_0^xi Phi(t) dt/sqrt(cosh xi - cosh t), t = xi (1 - s^2)."""
        from scipy.interpolate import CubicSpline
        s, ws = gl(0.0, 1.0, n); tt = xi*(1 - s*s)
        m = self.t <= xi + 0.01
        val = CubicSpline(self.t[m], Phi[m])(tt)
        den = 2*np.sinh(0.5*(xi + tt))*np.sinh(0.5*(xi - tt))      # cosh xi - cosh t without cancellation
        return (np.sqrt(2)/np.pi)*np.sum(ws*2*xi*s*val/np.sqrt(den))

def main(Dgs):
    pr("# h1: exact eta-inversion of psi_2 (route A', C_2, j = 1) via the Mehler-Dirichlet/Abel ODE solution")
    for Dg in Dgs:
        NG = NearGap(Dg); L = NG.Lam; F = NG.F
        c = np.array([float(v) for v in F.s2])
        inv = Inverter(c, L)
        # V1: symbol in the interior
        pv = np.array([0.2, 0.4, 0.6, 0.8]); nu = L*pv
        ex = inv.phit(nu) - inv.phi(nu)
        d1 = Pl.polyval(pv, Pl.polyder(c, 1))/L; d2 = Pl.polyval(pv, Pl.polyder(c, 2))/L**2; f0 = Pl.polyval(pv, c)
        sym = 0.5*d2 - d1/(2*nu) + f0/(2*nu**2)
        pr(f"  Delta_g = {Dg} (Lam = {L:.2f}; {len(inv.t)} t-points)")
        pr("    V1 interior, p = 0.2/0.4/0.6/0.8: exact dphi/phi = " + " ".join(f"{e/f:+.3e}" for e, f in zip(ex, f0)) +
           " ; symbol = " + " ".join(f"{s/f:+.3e}" for s, f in zip(sym, f0)))
        # V2: Abel identity cosh(xi) K~ = K
        for xi in (0.05, 0.3, 1.0, 2.5, 5.0):
            K = inv.abel(inv.Phi, xi); Kt = inv.abel(inv.Phit, xi)
            pr(f"    V2 xi = {xi}: K = {K:+.10e}, cosh(xi) K~ = {np.cosh(xi)*Kt:+.10e}, rel diff {abs(np.cosh(xi)*Kt - K)/abs(K):.1e}")
        # fine nu grid for integrals: [0, 30] dense, [30, Lam + 120]
        n1, w1 = gl(0.0, 30.0, 600); n2, w2 = gl(30.0, L, max(800, int(4*L))); n3, w3 = gl(L, L + 120.0, 600)
        nus = np.concatenate([n1, n2, n3]); wns = np.concatenate([w1, w2, w3])
        ph = inv.phi(nus); pt = inv.phit(nus)
        # V3: graviton-pole identity
        g = 1/(nus**2 + 2.25); gm = 1/((nus - 1j)**2 + 2.25); gt = np.real((1 + 1/(2j*nus))*gm)
        lhs = np.sum(wns*pt*gt); rhs = np.sum(wns*ph*g)
        pr(f"    V3 graviton pole: int phit g~ = {lhs:.12e}, int phi g = {rhs:.12e}, rel diff {abs(lhs - rhs)/abs(rhs):.1e}  "
           f"(without inversion: int phi g~ = {np.sum(wns*ph*gt):.6e}, rel {abs(np.sum(wns*ph*gt) - rhs)/rhs:.2e})")
        # dk0 (p = nu/Lam): flat |g3|^2 coefficient -2 int psi_2 p^6 (+ psi_3 part unchanged)
        p6 = (nus/L)**6/L
        dk_all = -2*np.sum(wns*(pt - ph)*p6); dk_in = -2*np.sum((wns*(pt - ph)*p6)[nus <= L])
        tail = np.sum((wns*np.abs(pt))[nus > L])/np.sum(wns*np.abs(ph))
        # symbol dk0 for comparison (as in f6)
        pf, wf = gl01(400); psf = Pl.polyval(pf, c)
        dps = (Pl.polyval(pf, Pl.polyder(c, 2)) - Pl.polyval(pf, Pl.polyder(c, 1))/pf + psf/pf**2)/(2*L**2)
        dk_sym = -2*np.sum(wf*dps*pf**6)
        pr(f"    dk0 exact = {dk_all:+.5e} (p <= 1 only: {dk_in:+.5e}); symbol estimate {dk_sym:+.5e}; tail mass int_(p>1)|psit|/int|psi| = {tail:.2e}")
        pr(f"    => A, B_AdS x {1/(1 + dk_all):.6f}; bounds x {(1 + dk_all)**-0.5:.6f}")
        mx = np.max(np.abs(pt - ph)[nus <= L])/np.max(np.abs(ph))
        lowmask = nus <= 5
        pr(f"    max|phit - phi|/max|phi| on [0, Lam] = {mx:.2e}; on nu <= 5: {np.max(np.abs(pt - ph)[lowmask])/np.max(np.abs(ph)):.2e}")
        pickle.dump(dict(Lam=L, t=inv.t, Phit=inv.Phit, nus=nus, wns=wns, phit=pt, phi=ph, dk0=dk_all, c=c),
                    open(os.path.join(HERE, 'cache', 'h1_%d.pkl' % Dg), 'wb'))
    open(os.path.join(HERE, 'cache', 'h1_etainv_exact_' + '_'.join(map(str, Dgs)) + '.txt'), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    main([int(a) for a in sys.argv[1:]] or [100, 300, 1000])
