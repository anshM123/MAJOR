"""
gm_chart (GAP_MORSE.md §2): certified holomorphic Morse chart of the Taylor POLYNOMIAL of the phase, and the Laplace remainder bound
(Lemma M), in Arb, from the Taylor-model data of gm_cert (Hermitian ball B_rho, Bombieri norms).

Phi - Phi* = P(v) + R(v):  P = sum_{k=2}^D Phi_k (exact Taylor polynomial, ball coefficients), |R| <= remPhi on B_rho, R = O(v^{D+1}).
Chart of P:  W(v) = S(v) v,  S = Q^{1/2} = sum_n binom(1/2,n)(Q - I)^n,  Q(v) = -2 sum_{k=2}^D Hess Phi_k(v)/(k(k-1))   (P = -W.W/2).
Bounds for |v| <= t (operator norms; [.] = Bombieri norm = Frobenius norm of the symmetric coefficient tensor):
  eta(t) = ||E'||_F + sum_{k=3}^D 2 [Phi_k] t^{k-2}             >= ||Q(v) - I||            (||Q_k(v)|| <= 2[Phi_k]|v|^{k-2})
  dQ(t)  = sum_{k=3}^D 2 (k-2) [Phi_k] t^{k-3}                  >= ||D_x Q(v)||, |x| = 1
  c(t)   = 1 - sqrt(1 - eta)                                    >= ||S(v) - I||            (N = W - id, |N(v)| <= c|v|)
  L(t)   = c + dQ t / (2 sqrt(1 - eta))                         >= ||DN(v)||
  (the same eta, c, L also bound the FROBENIUS norms of Q - I, S - I, DN: ||Q_k(v)||_F <= 2[Phi_k]|v|^{k-2} by Cauchy-Schwarz)
  T(t)   = |tr E'|/2 + sup|p_tr| + (c - eta/2) + dQ t (sqrt7 g1 + g2 - g1)   >= |tr DN(v)|,  p_tr = -sum_k Lap(Phi_k)/k,
           g1 = (c - eta/2)/eta, g2 = 1/(2 sqrt(1-eta)) - 1/2
  lamJ(t) = T + (-log(1-L) - L)                                 >= |log det W'(v)|  (|tr A^n| <= ||A||_F^2 ||A||^{n-2})
Inverse (contraction): eta(rho2), L(rho2) < 1, r = sqrt(1 - eta(rho2)) rho2  =>  V = W^{-1}: B_r -> B_rho2 holomorphic, |V(w)| <= |w|/(1-c2).
In w: integrand e^{-Delta w.w/2} e^{Delta R(V(w))} H(w), H = h(V) det V';  sup_{B_rc}|H| <= M_H(rc) = Mh(rv) e^{lamJ(rv)}, rv = rc/(1-c2).
Laplace (u = w/eps, Cauchy in eps on |eps'| = rc/|u|, rc = tau r):  F(eps,u) = H(eps u) e^{R(V(eps u))/eps^2},
  |F| <= M_H e^{kappa|u|^2},  kappa = remPhi (rv/rho)^{D+1}/rc^2;  R enters at order eps^{D-1} >= eps^5, so K2 = Gaussian moments of H_0, H_2, H_4.
  rel(Delta) <= 693 (1-2kappa)^{-13/2} M_H (eps/rc)^6/(1 - (s/rc)^2)                    [even orders >= 6 on the real ball |u| <= s/eps]
             +  M_H sum_{k=0,2,4} (eps/rc)^k E[|u|^k ; |u| > s/eps]                      [K2 polynomial outside the ball]
             +  M_H(r) |S^6| c(s) s^7 (1+L(s))^6 e^{-Delta a} (Delta/2pi)^{7/2},  a = s^2 beta/2 - R_s   [Stokes boundary]
  beta = 1 - 2c(s) - c(s)^2, R_s = remPhi (rv_s/rho)^{D+1}, rv_s = min(rho2, s(1+c(s))/(1-c2)); real contour = real v-ball of radius s,
  s(1 + c(s)) <= r.  It bounds |I_ball/((2pi/Delta)^{7/2} e^{Delta Phi*} g_env) - (c0 + c2 eps^2 + c4 eps^4)|, c_j the exact Laplace
  coefficients.  Every term is non-increasing in Delta once Delta a >= 7/2 (checked).
  + annulus s < |v| <= s_out (real v): t4 = Mh(s_out) (1-theta)^{-7/2} P(chi^2_7 > Delta (1-theta) s^2), theta = ||E'||_F +
  sum_k 2[Phi_k] s_out^{k-2} + 2 R(s_out)/s_out^2 (Re(Phi - Phi*) <= -(1-theta)|v|^2/2 on the real ball); so the certified statement is
  for the real Morse ball of radius s_out = the full certified radius (premise TAIL only beyond it).
"""
import math
from flint import arb

PI = arb.pi()
S6 = 2*PI**arb(3.5)/arb(3.5).gamma()          # area of the unit sphere S^6 in R^7


def chart(d, t):
    """chart bounds at radius t (dict of arb), or None if eta >= 1 or L >= 1."""
    t = arb(t); D = d['deg']
    eta = arb(d['EF']) + sum(2*arb(d['bPhi'][k])*t**(k - 2) for k in range(3, D + 1))
    if not (eta < 1):
        return None
    dQ = sum(2*(k - 2)*arb(d['bPhi'][k])*t**(k - 3) for k in range(3, D + 1))
    sq = (1 - eta).sqrt()
    c = 1 - sq
    L = c + dQ*t/(2*sq)
    if not (L < 1):
        return None
    ptr = sum(arb(d['bTr'][k])*t**(k - 2) for k in range(3, D + 1))
    g1 = (c - eta/2)/eta                      # sum_{n>=2} |b_n| eta^{n-1}
    g2 = 1/(2*sq) - arb(0.5)                  # sum_{n>=2} |b_n| n eta^{n-1}
    T = arb(d['trE'])/2 + ptr + (c - eta/2) + dQ*t*(arb(7).sqrt()*g1 + (g2 - g1))
    lamJ = T + (-(1 - L).log() - L)           # |sum_{n>=2} tr(A^n)/n| <= ||A||_F^2 sum ||A||^{n-2}/n, ||A||_F <= L
    return dict(t=t, eta=eta, dQ=dQ, c=c, L=L, T=T, lamJ=lamJ)


def Mh(d, t):
    """sup_{B_t} |h|, t <= rho_h."""
    t = arb(t); D = d.get('deg_h', d['deg'])
    assert t <= arb(d.get('rho_h', d['rho']))
    rh = arb(d.get('rho_h', d['rho']))
    return sum(arb(d['bh'][k])*t**k for k in range(D + 1)) + arb(d['remh'])*(t/rh)**(D + 1)


def Rb(d, t):
    """sup_{B_t} |R| (Schwarz), t <= rho."""
    t = arb(t)
    return arb(d['remPhi'])*(t/arb(d['rho']))**(d['deg'] + 1)


def MH(d, rv):
    C = chart(d, rv)
    return Mh(d, rv)*C['lamJ'].exp()


def setup(d, rho2, tau, sfrac):
    rmax = min(d['rho'], d.get('rho_h', d['rho']))
    if rho2 > rmax:
        return None
    C2 = chart(d, rho2)
    if C2 is None:
        return None
    rho2 = arb(rho2)
    r = (1 - C2['eta']).sqrt()*rho2
    rc = arb(tau)*r
    c2 = C2['c']
    rv = rc/(1 - c2)
    if not (rv < rho2):
        rv = rho2
    MHc = MH(d, rv)
    MHr = Mh(d, rho2)*C2['lamJ'].exp()
    kap = Rb(d, rv)/rc**2
    if not (kap < arb(0.25)):
        return None
    s = arb(sfrac)*rc
    Cs = chart(d, s)
    if Cs is None or not (s*(1 + Cs['c']) <= r):
        return None
    beta = 1 - 2*Cs['c'] - Cs['c']**2
    rvs = s*(1 + Cs['c'])/(1 - c2)
    if not (rvs < rho2):
        rvs = rho2
    a = s*s*beta/2 - Rb(d, rvs)
    if not (a > 0):
        return None
    # annulus s < |v| <= s_out of the real v-ball (s_out = full certified radius): Re(Phi - Phi*) <= -(1 - theta)|v|^2/2 there
    sout = arb(rmax)
    D = d['deg']
    def th(x):
        return arb(d['EF']) + sum(2*arb(d['bPhi'][k])*x**(k - 2) for k in range(3, D + 1)) + 2*Rb(d, x)/x**2
    theta = th(sout)
    if not (theta < 1) or not (s < sout):
        return None
    NS = 12
    edges = [s + (sout - s)*j/NS for j in range(NS)] + [sout]
    shells = [(edges[j], th(edges[j + 1]), Mh(d, edges[j + 1])) for j in range(NS)]
    return dict(rho2=rho2, tau=arb(tau), sfrac=arb(sfrac), r=r, rc=rc, rv=rv, MH=MHc, MHr=MHr, kap=kap, s=s, cs=Cs['c'], Ls=Cs['L'],
                beta=beta, a=a, eta2=C2['eta'], c2=c2, L2=C2['L'], lamJ2=C2['lamJ'], sout=sout, theta=theta, Mhout=Mh(d, sout), shells=shells)


_MOM = {k: arb(2)**(arb(k)/2)*(arb(7 + k)/2).gamma()/arb(3.5).gamma() for k in range(0, 44, 2)}   # E|u|^k, u ~ N(0, I_7)


def _chi_tail(k, U):
    """E[|u|^k ; |u| > U], u ~ N(0, I_7)."""
    a = arb(7 + k)/2
    return arb(2)**(arb(k)/2)*(U*U/2).gamma_upper(a)/arb(3.5).gamma()      # Gamma(a, U^2/2) (z = self, s = arg)


def rel(st, Delta):
    Dl = arb(Delta)
    eps = 1/Dl.sqrt()
    rc, M, s = st['rc'], st['MH'], st['s']
    x = eps/rc
    # even orders k = 6..K with full Gaussian moments E_kappa|u|^k = (1-2kappa)^{-(7+k)/2} E|u|^k, tail k > K by the geometric bound;
    # the minimum over K of these (each non-increasing in Delta) bounds the interior remainder
    om = 1 - 2*st['kap']
    t1 = None
    acc = arb(0)
    for K in range(6, 41, 2):
        acc = acc + M*x**K*om**(-arb(7 + K)/2)*_MOM[K]
        cand = acc + M*x**(K + 2)*om**(-arb(9 + K)/2)*_MOM[K + 2]/(1 - (s/rc)**2)
        if t1 is None or cand < t1:
            t1 = cand
    U = s/eps
    t2 = M*sum(x**k*_chi_tail(k, U) for k in (0, 2, 4))
    t3 = st['MHr']*S6*st['cs']*s**7*(1 + st['Ls'])**6*(-Dl*st['a']).exp()*(Dl/(2*PI))**arb(3.5)
    t4 = arb(0)
    for (sj, thj, Mj) in st['shells']:           # shell s_j < |v| <= s_{j+1}: weight <= Mh(s_{j+1}) e^{-Delta (1-theta_{j+1})|v|^2/2}
        X = Dl*(1 - thj)*sj*sj/2
        t4 = t4 + Mj*(1 - thj)**arb(-3.5)*X.gamma_upper(arb(3.5))/arb(3.5).gamma()
    mono = bool(Dl*st['a'] >= arb(3.5))
    return t1 + t2 + t3 + t4, dict(t1=t1, t2=t2, t3=t3, t4=t4, mono=mono)


def scan(d, f2=(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), taus=(0.7, 0.8, 0.9, 1.0), sfracs=(0.4, 0.5, 0.6, 0.7, 0.8)):
    rmax = min(d['rho'], d.get('rho_h', d['rho']))
    out = []
    for f in f2:
        for tau in taus:
            for sf in sfracs:
                st = setup(d, f*rmax, tau, sf)
                if st is not None:
                    out.append(st)
    return out


def delta1(d, target, lo=1e2, hi=1e10, sts=None):
    """smallest Delta on a log grid (40/decade) with certified rel <= target, every term non-increasing from there on
    (the predicate is monotone in Delta for each set-up, so a bisection per set-up is exact on the grid); best set-up."""
    import numpy as np
    sts = sts if sts is not None else scan(d)
    if not sts:
        return None, None
    grid = np.logspace(math.log10(lo), math.log10(hi), int(40*math.log10(hi/lo)) + 1)

    def ok(st, i):
        v, parts = rel(st, grid[i])
        return parts['mono'] and v.is_finite() and float(v.upper()) <= target

    best = None
    for st in sts:
        if not ok(st, len(grid) - 1):
            continue
        a, b = -1, len(grid) - 1                 # ok(b) true; find the first true
        while b - a > 1:
            m = (a + b)//2
            if ok(st, m):
                b = m
            else:
                a = m
        if best is None or b < best[0]:
            best = (b, st)
    if best is None:
        return None, None
    b, st = best
    v, parts = rel(st, grid[b])
    return float(grid[b]), dict(rel=float(v.upper()), parts={k: (float(p.upper()) if k != 'mono' else p) for k, p in parts.items()},
                                st={k: float(x.upper()) for k, x in st.items() if k != 'shells'})
