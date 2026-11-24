"""
gc_lemma (GAP_C3.md §2): explicit Laplace remainder bound (Lemma L3) from certified polydisc data.

Data (all certified by gc_cert, Delta-independent at fixed xi): rho, D (Taylor-model degree), Bphi[k] = sum_{|a|=k} |Phi_a| rho^k (k = 3..D,
centre Taylor coefficients, upper bounds), remPhi (sup_P |Phi - P_D|), Bh[k], remh (same for h = g/g_env), E (Morse defect).
For eps = Delta^{-1/2} and tau in (0, 1]:
  kappa(tau) = (sum_{k=3}^D Bphi_k tau^{k-2} + remPhi tau^{D-1}) / rho^2,     Mh(tau) = sum_k Bh_k tau^k + remh tau^{D+1},
  b = (1 - 7E)/2 - kappa(tau) > 0,   U = tau rho/(2 eps),
  rel(Delta) <= [ 2 Mh (eps/(tau rho))^6 * 7*15/(8 b^3) (pi/b)^{7/2}  +  6 Mh e^{-b U^2/2} U^{-5} * 112/b^3 (2 pi/b)^3 ] / (2 pi)^{7/2}
is a bound on |I_box(Delta) - I_K2(Delta)| / (Gaussian x g_env) for the integral over the real Morse box |v|_inf <= tau rho/2 (premise
TAIL for the rest of the cycle).  All arithmetic below is done in Arb (upper bounds).
"""
import math
from flint import arb


def rel_bound(d, Delta, tau, order=6):
    """d: dict with rho, deg, Bphi (list, index k), remPhi, Bh (list), remh, E.  Returns arb upper bound (or None if b <= 0).
    order = 6: remainder after the eps^4 (K2) term; order = 2: after the leading term."""
    rho = arb(d['rho']); D = d['deg']; t = arb(tau)
    kap = sum(arb(d['Bphi'][k])*t**(k - 2) for k in range(3, D + 1)) + arb(d['remPhi'])*t**(D - 1)
    kap = kap/rho**2
    Mh = sum(arb(d['Bh'][k])*t**k for k in range(0, D + 1)) + arb(d['remh'])*t**(D + 1)
    b = (1 - 7*arb(d['E']))/2 - kap
    if not (b > 0):
        return None
    eps = 1/arb(Delta).sqrt()
    tr = t*rho
    if not (2*eps < tr):
        return None
    U = tr/(2*eps)
    pi = arb.pi()
    n = order
    # moment: int e^{-b|u|^2} |u|_inf^n <= 7 * E|N(0,1/(2b))|^n * (pi/b)^{7/2},  E|N|^n = (n-1)!!/(2b)^{n/2}
    dfact = 1
    for j in range(n - 1, 0, -2):
        dfact *= j
    mom = 7*arb(dfact)/(2*b)**(n//2)*(pi/b)**arb(3.5)
    t1 = 2*Mh*(eps/tr)**n*mom
    # tail of the truncated series outside the box: sum_{k<n} eps^k |F_k| <= n Mh e^{kappa|u|^2} (|u|_inf/U)^(n-1)
    m = n - 1
    c = b/2
    one_d = (arb(m + 1)/2).gamma()/c**(arb(m + 1)/2)          # int_R |t|^m e^{-c t^2} dt
    t2 = n*Mh*(-b*U**2/2).exp()/U**m*7*one_d*(pi/c)**3
    return (t1 + t2)/(2*pi)**arb(3.5)


def best(d, Delta, order=6, taus=None):
    taus = taus or [0.05*i for i in range(1, 21)]
    out = None
    for tau in taus:
        r = rel_bound(d, Delta, tau, order)
        if r is not None and r.is_finite():
            v = float(r.upper())
            if out is None or v < out[0]:
                out = (v, tau)
    return out


def delta1(d, target, order=6, lo=1e3, hi=1e14):
    """smallest Delta (on a log grid, 40 pts/decade) with rel(Delta) <= target; rel is checked to be decreasing beyond it."""
    import numpy as np
    grid = np.logspace(math.log10(lo), math.log10(hi), int(40*math.log10(hi/lo)) + 1)
    for D in grid:
        r = best(d, D, order)
        if r is not None and r[0] <= target:
            return float(D), r
    return None, None
