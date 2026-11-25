"""
epswcert/ew_cert.py (EPSW_CERT.md §2): ball-arithmetic certificate for the published k = 2 continuum weight of Omega^EE at one gap.

Certifies, for r(mu) = E(mu) - psi(mu) (ew_common), the numbers that control W[chi] - Phi = int_0^oo r(mu) cos(mu t) dmu:
  r(0)                     (ball)                          -> int_0^oo (W - Phi) dt = (pi/2) r(0)
  r'(0+)                   (ball)
  eps_W  >= int_0^oo |r|                                   -> |W - Phi|(t) <= eps_W                for all t >= 0
  D2     >= int_0^oo |r''|  (r'' in L^1: log singularity at 0)
  B := |r'(0+)| + D2                                       -> |W - Phi|(t) <= B/t^2              for all t > 0
     (two integrations by parts: W - Phi = -[r'(0+) + int_0^oo r'' cos(mu t) dmu]/t^2; r, r' -> 0 at infinity)
  sup |r|.
Method (nothing sampled):
  * [0, A0], A0 = 1/5: f = F1 + F2 log mu with F1, F2 analytic in |mu| < 1 (ew_common); their Taylor coefficients at 0 to order N0
    are computed exactly in ball arithmetic (power-series arithmetic), the truncation is bounded by Cauchy with a maximum-modulus
    enclosure on |mu| = RHO0 (ring of boxes).  2 Re F2 has no mu^0, mu^1 terms (F2(0) = 0; F2'(0) = -2i Y(-i) R_-(-i), with
    Y(-i), R_-(-i) = 18/25 real) [P].  Cells of width 1/500, polynomial Taylor shifts, explicit log factors.
  * [A0, mu_end], mu_end = ceil(Lam) + 60: panels (width 1/40 .. 1/10; split at Lam where psi stops) with Taylor series of
    f_near to order N at the panel centre (exact series arithmetic), truncation by Cauchy (radius rho, disc sup M by explicit RBF
    bounds x box enclosures), far RBF terms bounded on the disc (F), subcells (5 per panel) by polynomial Taylor shift.
  * (mu_end, oo): explicit Gaussian tail (erfc) and Cauchy for r''.
Usage: python3 ew_cert.py DG [re|im] [pub|<wtag>]    -> cache/ew_cert_<DG>_<part>[_<wtag>].{txt,json}
"""
import os, sys, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, acb, arb_poly
import ew_common as C
from ew_common import up

N = 40            # panel Taylor order
Q = 5             # subcells per panel
N0 = 64           # order at mu = 0
RHO0 = 0.8
A0 = Fr(1, 5)
H0 = Fr(1, 500)   # cell width on [0, A0]


def fa(x):
    x = Fr(x)
    return arb(x.numerator)/x.denominator


def abs_poly_bounds(cf, hs):
    """for p(s) = sum cf[n] s^n on |s| <= hs: (sup |p| upper, int_{-hs}^{hs} |p| upper) as arbs."""
    sup = arb(0); itg = arb(0); hp = arb(1); hs = arb(hs)
    for n, c in enumerate(cf):
        a = abs(c).upper()
        sup += a*hp
        itg += a*2*hp*hs/(n + 1)
        hp = hp*hs
    return sup, itg


def coeffs(p, n):
    return [p[i] for i in range(n)]


def shift(p, c):
    """Taylor shift of arb_poly p to centre c."""
    return p(arb_poly([c, 1]))


def deriv2(cf):
    return [cf[n]*n*(n - 1) for n in range(2, len(cf))] or [arb(0)]


class Acc:
    def __init__(self):
        self.I1 = arb(0); self.I2 = arb(0); self.sup = arb(0); self.by = {}

    def add(self, i1, i2, sp, a):
        self.I1 += i1; self.I2 += i2
        self.sup = C.amax(self.sup, sp)
        key = 'mu<1' if a < 1 else ('1-40' if a < 40 else ('40-Lam' if a < self.Lam else '>Lam'))
        b = self.by.setdefault(key, [arb(0), arb(0)])
        b[0] += i1; b[1] += i2


def near_zero(W, pr, part):
    """[0, A0]: returns r(0), r'(0+), int|r|, int|r''|, sup|r| (arbs)."""
    F1, F2 = C.F12_series(W, N0)
    M1, M2 = C.F12_disc_sup(W, RHO0)
    P1 = [2*c.real for c in F1]
    P2 = [2*c.real for c in F2]
    assert F2[0] == 0
    assert abs(F2[1].real) < arb('1e-40') and F2[1].real.contains(0)
    P2[0] = arb(0); P2[1] = arb(0)          # exact [P]: F2(0) = 0, Re F2'(0) = 0
    psi = W.psi_mu + [arb(0)]*(N0 - len(W.psi_mu)) if part == 're' else [arb(0)]*N0
    A = [P1[n] - psi[n] for n in range(N0)]
    r0 = A[0]; r1 = A[1]
    pr(f'  near 0: F1, F2 to order {N0}; disc sups on |mu| = {RHO0}: M1 <= {M1.str(5)}, M2 <= {M2.str(5)}')
    pr(f'  r(0) = {r0.str(15)};  r\'(0+) = {r1.str(15)}  (kappa = psi\'(0) = {W.kappa.str(15)})')
    # A2 = A'' + sum_{n>=2} (2n-1) P2_n mu^(n-2)  (the analytic part of r''), P2'' (multiplies log mu)
    A2 = deriv2(A)
    for n in range(2, N0):
        A2[n - 2] += (2*n - 1)*P2[n]
    P2d2 = deriv2(P2)
    pA, pP2, pA2, pP2d2 = (arb_poly(x) for x in (A, P2, A2, P2d2))
    rho = arb(RHO0)
    I1 = arb(0); I2 = arb(0); sup = arb(0)
    # cells: [0, U0], geometric [U0 2^j, U0 2^(j+1)] up to H0, then uniform width H0 to A0; on every cell [u, v] with u > 0 the
    # centre c and half-width hc satisfy hc/c <= 1/3, and log mu = log c + log(1 + s/c) is replaced by its degree-KL Taylor
    # polynomial (|remainder| <= x^(KL+1)/((KL+1)(1-x)), x = hc/c) so that the cancellation between A and P2 log mu is kept.
    U0 = Fr(H0, 2**31)
    edges = [Fr(0), U0]
    while edges[-1] < H0:
        edges.append(edges[-1]*2)
    assert edges[-1] == H0
    k = 1
    while Fr(k + 1)*H0 <= A0:
        edges.append(Fr(k + 1)*H0); k += 1
    assert edges[-1] == A0
    KL = 48
    for i in range(len(edges) - 1):
        u = fa(edges[i]); v = fa(edges[i + 1]); c = (u + v)/2; hc = (v - u)/2
        x = v/rho
        tr1 = 2*M1*x**N0/(1 - x)                                         # |T1|, T1 = (2 Re) F1 tail
        tr2 = 2*M2*x**N0/(1 - x)                                         # |T2|
        g2 = N0*(N0 - 1)*x**(N0 - 2)*(1 + x)/(1 - x)**3/rho**2           # sum_{n>=N0} n(n-1) mu^(n-2)/rho^n
        logv = -(v.log())
        if i == 0:
            # [0, U0]: |P2(mu) log mu| <= v^2 |log v| sum_{n>=2} |P2_n| v^(n-2)   (mu^2 |log mu| increasing on (0, e^-1/2))
            sA, iA = abs_poly_bounds(coeffs(shift(pA, c), N0), hc)
            s2 = sum((abs(P2[n]).upper()*v**(n - 2) for n in range(2, N0)), arb(0))
            sL = v*v*logv*s2
            eT = tr1 + tr2*logv                                          # mu^N0 |log mu| increasing on (0, e^(-1/N0))
            supc = sA + sL + eT
            I1 += iA + 2*hc*(sL + eT)
            sA2, iA2 = abs_poly_bounds(coeffs(shift(pA2, c), N0), hc)
            sPd, _ = abs_poly_bounds(coeffs(shift(pP2d2, c), N0), hc)
            intlog = v - v*v.log()                                       # int_0^v |log mu|
            e2 = 2*M1*g2 + 2*M2*2*(logv + 2)*g2
            I2 += iA2 + sPd*intlog + 2*hc*e2
        else:
            xl = hc/c
            assert xl <= arb(1)/3 + arb('1e-30')
            L = arb_poly([c.log()] + [(-1)**(j + 1)/(arb(j)*c**j) for j in range(1, KL + 1)])
            el = xl**(KL + 1)/((KL + 1)*(1 - xl))
            lu = -(u.log())
            sP2, _ = abs_poly_bounds(coeffs(shift(pP2, c), N0), hc)
            sPd, _ = abs_poly_bounds(coeffs(shift(pP2d2, c), N0), hc)
            rp = shift(pA, c) + shift(pP2, c)*L
            r2p = shift(pA2, c) + shift(pP2d2, c)*L
            s1, i1 = abs_poly_bounds(coeffs(rp, rp.degree() + 1), hc)
            s2_, i2 = abs_poly_bounds(coeffs(r2p, r2p.degree() + 1), hc)
            eT = tr1 + tr2*lu + sP2*el
            e2 = 2*M1*g2 + 2*M2*2*(lu + 2)*g2 + sPd*el
            supc = s1 + eT
            I1 += i1 + 2*hc*eT
            I2 += i2 + 2*hc*e2
        sup = C.amax(sup, supc.upper())
    ncell = len(edges) - 1
    pr(f'  [0, {float(A0)}]: int|r| <= {up(I1):.4e}, int|r\'\'| <= {up(I2):.4e}, sup|r| <= {up(sup):.4e}  ({ncell} cells)')
    return r0, r1, I1, I2, sup


def panels(W):
    Lam = W.Lam
    E = [Fr(1, 5) + Fr(k, 40) for k in range(13)]                   # 0.2 .. 0.5 by 1/40
    E += [Fr(1, 2) + Fr(k, 20) for k in range(1, 11)]               # .. 1 by 1/20
    Lf = Fr(int(float(Lam.lower())*10), 10)                         # < Lam
    k = 1
    while Fr(1) + Fr(k, 10) < Lf:
        E.append(Fr(1) + Fr(k, 10)); k += 1
    E1 = [fa(e) for e in E] + [Lam]
    mu_end = Fr(int(float(Lam.upper())) + 1 + 60)
    Lc = Fr(int(float(Lam.upper())*10) + 1, 10)                     # > Lam
    E2 = [Lam]
    x = Lc
    while x < mu_end:
        E2.append(fa(x)); x += Fr(1, 10)
    E2.append(fa(mu_end))
    P = [(E1[i], E1[i + 1], True) for i in range(len(E1) - 1)] + [(E2[i], E2[i + 1], False) for i in range(len(E2) - 1)]
    return P, mu_end


def tail(W, mu_end, part):
    """int_{mu_end}^oo |r| and |r''| (psi = 0 there).  On |z - mu| <= 1/2 (mu >= mu_end):
    |f| <= qbar * sum_k |x_k| exp((9/4 - (mu - 1/2 - c_k)^2)/(2 s_k^2)) + Ybound,  |Q R_-| <= qbar(R0) (product of decreasing
    ratios, R0 = mu_end - 1/2 <= |nu|);  |E| <= 2|f|, |E''| <= 2 * 2 sup|f|/(1/2)^2 = 16 sup|f| (Cauchy)."""
    R0 = fa(mu_end) - arb('0.5')
    qbar = 32*R0*(2*R0 + 3)*(2*R0 + 5)*(2*R0 + 11)/((2*R0 - 1)*(2*R0 - 7)*(2*R0 - 3)*(2*R0 - 7))
    tot = arb(0)
    for (c, s, x), (ca, sa, xa, inv2s2) in zip(W.sc, W.sc_arb):
        z = (fa(mu_end) - arb('0.5') - ca)/(sa*arb(2).sqrt())
        assert z > 0
        tot += abs(xa)*(arb(9)/4*inv2s2).exp()*sa*(arb.pi()/2).sqrt()*z.erfc()
    rbf = qbar*tot
    # Y-terms: |e^{-nu^2/2}| <= exp(-((mu-1/2)^2 - 9/4)/2); |nu| <= mu + 2; |nu^2 + 1| <= (mu+2)^2 + 1; |Log| <= log((mu+2)^2+1) + pi;
    # |R_-| <= 2 (2|nu|+11)^6 / (|nu| (2|nu|-7)^3) with |nu| >= R0.  For mu >= mu_end >= 358 the product P(mu) e^{-(mu-1/2)^2/2 + 9/8}
    # is <= e^{-(mu - 1/2)^2/4} (P(mu) <= (3 mu)^20 << e^{(mu-1/2)^2/4}), and int_{mu_end}^oo e^{-(mu-1/2)^2/4} <= (2/R0) e^{-R0^2/4}.
    ysum = sum((abs(y) for y in W.y), arb(0))
    m2 = fa(mu_end) + 2
    Pb = ysum*m2**4*(m2*m2 + 1)*((m2*m2 + 1).log() + arb.pi())*2*(2*m2 + 11)**6/(R0*(2*R0 - 7)**3)
    assert Pb*(arb(9)/8).exp() < (R0*R0/4).exp()        # and log P grows only logarithmically: holds for all mu >= mu_end
    ytail = 2/R0*(-(R0*R0)/4).exp()
    I1 = 2*(rbf + ytail); I2 = 16*(rbf + ytail)
    return I1, I2


def main():
    Dg = int(sys.argv[1]); part = sys.argv[2] if len(sys.argv) > 2 else 're'
    wtag = sys.argv[3] if len(sys.argv) > 3 else 'pub'          # 'pub' = prl/data/k2 (published); else cache/k2_continuum_D_<Dg>_<wtag>.json
    path = None if wtag == 'pub' else os.path.join(C.CACHE, f'k2_continuum_D_{Dg}_{wtag}.json')
    t0 = time.time()
    lines = []; pr = C.pr_to(lines)
    pr(__doc__.split('Usage')[0])
    assert C.check_kernels()
    W = C.Weight(Dg, path)
    pr(f'weight file: {os.path.relpath(W.path, C.ADS4)}')
    if part == 'im':
        d = W.json
        xs = [float(s) for s in d['x_im']]
        sc = []
        for (c, s), x in zip(d['centres_sigma'], xs):
            sc.append((c, s, x))
            if c > 0: sc.append((-c, s, x))
        sc.sort(); W.sc = sc; W.cs_sorted = [c for c, s, x in sc]
        W.sc_arb = [(arb(c), arb(s), arb(x), 1/(2*arb(s)**2)) for c, s, x in sc]
        W.sumabs = {}
        for c, s, x in sc:
            W.sumabs[s] = W.sumabs.get(s, arb(0)) + abs(arb(x))
        W.y = [arb(float(s)) for s in d['y_im']]
    pr(f'Dg = {Dg}, part = {part}: Lam = {W.Lam.str(20)}, {len(W.sc)} signed RBF terms, y = {[float(y.mid()) for y in W.y]}')
    r0, r1, I1, I2, sup = near_zero(W, pr, part)
    acc = Acc(); acc.Lam = float(W.Lam.mid())
    P, mu_end = panels(W)
    pr(f'  panels on [{float(A0)}, {float(mu_end)}]: {len(P)} (N = {N}, {Q} subcells each)')
    maxratio = 0.0; maxM = 0.0; maxerr = arb(0)
    for ip, (a, b, inside) in enumerate(P):
        m = (a + b)/2; h = arb(((b - a)/2).upper())
        mf = float(m.mid()); af = float(a.lower())
        rho = min(0.5, 0.6*mf)
        rad = float(m.rad()) + 1e-30
        Mf, Ff = C.f_disc_sup(W, mf, rho + rad)
        x = h/rho
        maxratio = max(maxratio, up(x)); maxM = max(maxM, up(Mf))
        fc = C.f_series(W, m, N)
        E = [2*c.real for c in fc] if True else None
        if inside and part == 're':
            ps = W.psi_taylor(m, N)
            rc = [E[n] - ps[n] for n in range(N)]
        else:
            rc = E
        Mn = Mf + Ff
        er = 2*(Mn*x**N/(1 - x) + Ff)
        er2 = 2*(Mn*N*(N - 1)*x**(N - 2)*(1 + x)/(1 - x)**3/arb(rho)**2 + 2*Ff/(arb(rho) - h)**2)
        maxerr = C.amax(maxerr, er.upper())
        pr_ = arb_poly(rc); pr2 = arb_poly(deriv2(rc))
        hs = h/Q
        i1 = arb(0); i2 = arb(0); sp = arb(0)
        for j in range(Q):
            cj = -h + (2*j + 1)*hs
            s1, j1 = abs_poly_bounds(coeffs(shift(pr_, cj), N), hs)
            s2, j2 = abs_poly_bounds(coeffs(shift(pr2, cj), N - 2), hs)
            i1 += j1 + 2*hs*er; i2 += j2 + 2*hs*er2
            sp = C.amax(sp, (s1 + er).upper())
        acc.add(i1, i2, sp, af)
        if ip % 500 == 0:
            pr(f'    panel {ip}/{len(P)} at mu = {mf:.3f}: int|r| so far {up(I1 + acc.I1):.4e}, int|r\'\'| so far {up(I2 + acc.I2):.4e} '
               f'[{time.time() - t0:.0f}s]')
    T1, T2 = tail(W, mu_end, part)
    eps = I1 + acc.I1 + T1
    D2 = I2 + acc.I2 + T2
    supr = C.amax(sup, acc.sup)
    Bc = abs(r1).upper() + D2
    pr(f'  panels: max h/rho = {maxratio:.3f}, max disc sup M = {maxM:.3e}, max truncation+far error per point <= {up(maxerr):.2e}')
    for k_, (a_, b_) in sorted(acc.by.items()):
        pr(f'    range {k_:>7}: int|r| <= {up(a_):.4e}, int|r\'\'| <= {up(b_):.4e}')
    pr(f'  tail mu > {float(mu_end)}: int|r| <= {up(T1):.3e}, int|r\'\'| <= {up(T2):.3e}')
    pr('')
    pr(f'RESULT Dg = {Dg} ({part}):  r(0) = {r0.str(12)},  r\'(0+) = {r1.str(12)}')
    pr(f'RESULT Dg = {Dg} ({part}):  eps_W = int|r| <= {up(eps):.6e};  int|r\'\'| <= {up(D2):.6e};  B = |r\'(0+)| + int|r\'\'| <= {up(Bc):.6e};'
       f'  sup|r| <= {up(supr):.4e}   [{time.time() - t0:.0f}s]')
    out = dict(Dg=Dg, part=part, r0_mid=float(r0.mid()), r0_rad=float(r0.rad()), r0_lo=float(r0.lower()), r0_hi=float(r0.upper()),
               r1_lo=float(r1.lower()), r1_hi=float(r1.upper()), epsW=up(eps), D2=up(D2), B=up(Bc), supr=up(supr),
               kappa=float(W.kappa.mid()), Lam=float(W.Lam.mid()))
    tag = f'{Dg}_{part}' + ('' if wtag == 'pub' else '_' + wtag)
    json.dump(out, open(os.path.join(C.CACHE, f'ew_cert_{tag}.json'), 'w'), indent=1)
    open(os.path.join(C.CACHE, f'ew_cert_{tag}.txt'), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
