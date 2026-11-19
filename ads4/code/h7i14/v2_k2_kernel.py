"""v2 (H7_I14_LIGHTSIDE.md §14): large-eta decay of the total de-shifted C2 kernel (K2 at eta -> oo), with a rigorous bound.

Object.  The leading-Regge eta-kernel of int_0^Lam dnu psi2(nu/Lam) C_{2,nu} (CHMRS-type, de-shifted):
     K_w(beta) = int_0^Lam dnu w(nu) P_{-1/2+i nu}(cosh beta),   eta = cosh beta,
with w(nu) = psi2(nu/Lam) (heavy / Regge-moment normalisation) or w(nu) = psi2(nu/Lam) a/gamma_7(nu)^2 (dDisc-space normalisation, CHMRS
eq. C_{k,nu}); in the CLS-D building-block form the same function is eta x [eta-inverted superposition of P_{l'}] (check (E) of v2_d3:
the magic-block kernel g P_{l-1} + g~ P_{l+1} is the heavy image of eta P_l), so the TOTAL kernel is K_w itself.

Proposition K (proved in §14; checked here).  Let w be C^4 on [0, Lam] with w(0) = 0, w'(Lam) = w''(Lam) = 0 (psi2 has (1-p)^5).  Put
  M0 = ||w||_1,  M2 = |w'(0)| + ||w''||_1,  c4 = |w'''(0)| + ||w''''||_1.  Then for beta >= 2
  | (pi/2) e^{beta/2} K_w(beta) - w'(0)/beta |  <=  eps(beta) := c4/(3 beta^3) + 2 M0 (e^{-beta/2} + beta e^{-beta}/2)
                                                       + (8 M2/beta^2)(e^{-1} + beta e^{-beta}/2) + 2.6 M2/(beta-1)^2.
Hence |K_w(eta)| <= (2/pi) eta^{-1/2} [ |w'(0)|/log eta + O(1/log^2 eta) ]  (e^{-beta/2} <= eta^{-1/2}, beta >= log eta).
The leading coefficient w'(0) = X(0)/(k0 Lam) != 0 for NF02: the decay is eta^{-1/2}/log eta, NOT eta^{-1/2}/log^2 eta (§7.2 corrected).

Numerics: (a) K_w by the Mehler-Dirichlet form (Phi(t) = int w cos(nu t) exact by terminating integration by parts for t >= 1), checked
against direct nu-integration with mpmath at beta = 1, 3, 6; (b) (pi/2)e^{beta/2}K_w beta vs w'(0) up to beta = 400; (c) the bound eps;
(d) per-nu (unsmeared) kernel P_l and the individual shifted kernels eta P_{l +- 1} for contrast; (e) the r-integrated large-eta tail.
Output: cache/v2_k2_kernel.txt
"""
import os, sys, io
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
import numpy as np
import mpmath as mp
from fractions import Fraction
from numpy.polynomial import polynomial as Pn
from numpy.polynomial.legendre import leggauss
import json

buf = io.StringIO()
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")

# ---------------- NF02 psi2 (exact rationals; read-only json)
J = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_noforward.json')))['NF02']
from nf_common import NF                       # read-only import: k0
k0 = float(NF('NF02').k0)
X = [float(Fraction(c)) for c in J['X']]
psi2 = Pn.polymul(Pn.polymul(Pn.polypow([1, -1], 5), [0, 1]), X)/k0      # ascending coefficients in p
pr(f"NF02: k0 = {k0:.10f}; psi2'(0) = X(0)/k0 = {psi2[1]:.8f}; check psi2_exact[1]/k0 = {float(Fraction(J['psi2_exact'][1]))/k0:.8f}")

def wpoly(Lam):
    """w(nu) = psi2(nu/Lam) as ascending polynomial in nu."""
    return psi2/Lam**np.arange(len(psi2))

def Phi(t, c, Lam):
    """int_0^Lam w(nu) cos(nu t) dt, w = poly c.  t >= 1: terminating integration by parts; t < 1: Gauss-Legendre."""
    t = np.atleast_1d(np.asarray(t, float)); out = np.empty_like(t)
    big = t >= 1.0
    if big.any():
        tb = t[big]; acc = np.zeros_like(tb, dtype=complex)
        ders = [c]
        while len(ders[-1]) > 0 and np.any(ders[-1] != 0):
            ders.append(Pn.polyder(ders[-1]))
            if len(ders) > len(c) + 1:
                break
        for k, dk in enumerate(ders):
            if len(dk) == 0:
                break
            val = lambda nu: Pn.polyval(nu, dk)
            term = (val(Lam)*np.exp(1j*Lam*tb) - val(0.0))*(-1)**k/(1j*tb)**(k + 1)
            acc += term
        out[big] = acc.real
    if (~big).any():
        x, wq = leggauss(400); nu = 0.5*Lam*(x + 1); wn = 0.5*Lam*wq
        out[~big] = (Pn.polyval(nu, c)[None, :]*np.cos(np.outer(t[~big], nu))*wn).sum(1)
    return out

_g20, _w20 = leggauss(20)
def _sgrid(npan):
    e = np.linspace(0, 1, npan + 1); a, b = e[:-1, None], e[1:, None]
    return (0.5*(b - a)*(_g20 + 1) + a).ravel(), (0.5*(b - a)*_w20).ravel()

def Kscaled(beta, c, Lam):
    """(pi/2) e^{beta/2} K_w(beta), via t = beta(1-s^2):  (pi/2)e^{beta/2}(sqrt2/pi) int Phi(t) (cosh b - cosh t)^{-1/2} dt
       = int Phi(t) kappa(t, beta) dt,  kappa = [(1-e^{t-b})(1-e^{-t-b})]^{-1/2}.  Composite 20-point Gauss-Legendre in s with
       ~ Lam*beta/8 panels (Phi carries e^{i Lam t} terms from the nu = Lam endpoint)."""
    s, ws = _sgrid(int(Lam*beta/8) + 80)
    t = beta*(1 - s*s)
    one_m_x = -np.expm1(-beta*s*s); one_m_y = -np.expm1(-t - beta)
    kap = 1/np.sqrt(one_m_x*one_m_y)
    return (Phi(t, c, Lam)*kap*2*beta*s*ws).sum()

def Kdirect_mp(beta, c, Lam):
    mp.mp.dps = 20
    f = lambda nu: Pn.polyval(float(nu), c)*mp.re(mp.legenp(-0.5 + 1j*nu, 0, mp.cosh(beta), type=3))
    pts = np.linspace(0, Lam, 60)
    return float(mp.quad(f, list(pts)))

def norms(c, Lam):
    x, wq = leggauss(2000); nu = 0.5*Lam*(x + 1); wn = 0.5*Lam*wq
    d = lambda k: Pn.polyder(c, k) if k else c
    L1 = lambda k: (np.abs(Pn.polyval(nu, d(k)))*wn).sum()
    M0 = L1(0); M2 = abs(Pn.polyval(0, d(1))) + L1(2); c4 = abs(Pn.polyval(0, d(3))) + L1(4)
    return M0, M2, c4, Pn.polyval(0, d(1)), Pn.polyval(Lam, d(1)), Pn.polyval(Lam, d(2))

def eps_bound(beta, M0, M2, c4):
    return (c4/(3*beta**3) + 2*M0*(np.exp(-beta/2) + beta*np.exp(-beta)/2) + (8*M2/beta**2)*(np.exp(-1) + beta*np.exp(-beta)/2)
            + 2.6*M2/(beta - 1)**2)

for Dg in (100, 300):
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); c = wpoly(Lam)
    M0, M2, c4, w1, w1L, w2L = norms(c, Lam)
    pr(f"\nDelta_g = {Dg}, Lam = {Lam:.4f}: w'(0) = {w1:.6e}; w'(Lam) = {w1L:.1e}, w''(Lam) = {w2L:.1e};  M0 = {M0:.4e}, M2 = {M2:.4e}, c4 = {c4:.4e}")
    for b in (1.0, 3.0, 6.0):
        kd = Kdirect_mp(b, c, Lam); ks = Kscaled(b, c, Lam)*2/np.pi*np.exp(-b/2)
        pr(f"   (a) beta = {b}: Mehler {ks:+.12e}, direct nu-integral (mpmath) {kd:+.12e}, rel diff {abs(ks/kd - 1):.1e}")
    pr("   (b,c) beta   (pi/2)e^{b/2}K*beta/w'(0)   |(pi/2)e^{b/2}K - w'(0)/b|*b^2/w'(0)   eps(beta)*b^2/|w'(0)|   bound holds")
    for b in (2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 400.0):
        ks = Kscaled(b, c, Lam)
        rem = abs(ks - w1/b)
        e = eps_bound(b, M0, M2, c4)
        pr(f"        {b:6.1f}   {ks*b/w1:+.8f}            {rem*b*b/abs(w1):10.4f}                     {e*b*b/abs(w1):12.2f}         {rem <= e}")

# ---------------- (d) per-nu and shifted kernels
pr("\n(d) envelopes at large beta: per-nu P_l (no smearing), and the individual CLS-D building blocks eta P_{l+-1} (shifted lines)")
mp.mp.dps = 20
for b in (5.0, 10.0, 20.0, 40.0):
    x = mp.cosh(b)
    p0 = mp.re(mp.legenp(-0.5, 0, x, type=3)); p1 = abs(mp.legenp(-0.5 + 1j*1.0, 0, x, type=3))
    pm = abs(x*mp.legenp(-0.5 + 1j*(1.0 - 1j), 0, x, type=3))
    pr(f"   beta = {b}: e^{{b/2}} P_(-1/2)(cosh b)/b = {mp.nstr(p0*mp.e**(b/2)/b, 6)} (-> 2/pi = {mp.nstr(2/mp.pi, 6)}, as (2/pi)(1 + 2 log2/b));"
       f"  e^{{b/2}}|P_(-1/2+i)| = {mp.nstr(p1*mp.e**(b/2), 6)};  e^{{-3b/2}}|eta P_(l+1)|(nu=1) = {mp.nstr(pm*mp.e**(-1.5*b), 6)}")

# ---------------- (e) r-integrated tail
pr("\n(e) r-integrated large-eta tail (leading-Regge kernel, k = 2): T(eta1) = int_{eta1}^oo 2 deta |K_w(eta)| int_0^{c e^{-beta}} dr r^{1-J0},"
   " J0 = 1, c = 1/2, Delta_g = 100 (heavy normalisation)")
Lam = np.sqrt(98.0*99.0); c = wpoly(Lam)
bs = np.linspace(2.0, 60.0, 600); Ks = np.array([abs(Kscaled(b, c, Lam))*2/np.pi*np.exp(-b/2) for b in bs])
integrand = 2*Ks*np.sinh(bs)*(0.5*np.exp(-bs))**(2 - 1)/(2 - 1)
cum = np.cumsum(integrand[::-1]*(bs[1] - bs[0]))[::-1]
for b1 in (2.0, 5.0, 10.0, 20.0, 40.0):
    i = np.searchsorted(bs, b1)
    pr(f"   eta1 = cosh({b1}): T(eta1) ~ {cum[i]:.4e}   (asymptotic (2/pi) |w'(0)| int beta^-1 e^{{-beta/2}} dbeta: integrable, margin e^{{-beta/2}} = eta^{{-1/2}})")
open(os.path.join(HERE, 'cache', 'v2_k2_kernel.txt'), 'w').write(buf.getvalue())
