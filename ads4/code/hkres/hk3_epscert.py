r"""
hk3_epscert (HK_RESIDUAL.md §3.3): INTERVAL-RIGOROUS upper bound on eps_W = int_0^oo |E_even(mu)| dmu for the published continuum C4
de-shift weights (prl/data/q_deshift_continuum_{c4v,c4d}.json: chi = Q B, B = Gaussian RBF sum with the 60-digit decimal coefficients,
h, sigma = the IEEE doubles nearest 0.15, 0.45, exactly as
ab_solve.py used them; with the decimal values the residue part moves by ~1e-16), plus rigorous enclosures of the residue part, of a2, and of the tail beyond mu = L.

Ball arithmetic: python-flint (arb/acb), 384 bits.  Nothing is a midpoint estimate: every quantity is an arb ball or an explicit bound.
(A) Laurent data at mu = 0 of f_{+-3} (needed for S_j and for the residue part): Taylor coefficients of g = mu^p f_j by the K0-point
    trapezoid rule on |mu| = r0, with the aliasing error bounded by M(rho0) rho0^-n (r0/rho0)^K0 / (1 - (r0/rho0)^K0), where M(rho0) is an
    enclosure of max |g| on |mu| = rho0 (circle covered by boxes, interval evaluation).  g is analytic in |mu| < 1/2 (the u-poles at
    nu = +-5i/2, +-i(3+2 sqrt2)/2 sit at zeros of the crossed R_j; the second cancellation is done EXACTLY below).
(B) eps_W on [0, L]: panels of width w; on each, E_even(m + s) = sum_n c_n s^n with c_n from a K-point trapezoid on |z - m| = r
    (aliasing bound as in (A)) and |c_n| <= M_rho / rho^n (Cauchy) for the truncation; M_rho = termwise interval enclosure of max |E_even|
    on the boundary of a chunk rectangle containing all circles |z - m| = rho (maximum modulus).  M is large (1e9-1e13: the termwise
    enclosure does not see the de-shift cancellation), so K is chosen per chunk with M (r/rho)^K <= 1e-28.  E_even is analytic in |Im mu| < 1/2 (target pole at 0 removed by S_j with the ball-valued a_k):
        int_panel |E_even| <= sum_{n<=N} (|c^_n| + e_n) 2 (w/2)^{n+1}/(n+1) + w M_rho (w/(2 rho))^{N+1}/(1 - w/(2 rho)).
(C) tail mu > L: |E| <= sum_j |Q R_j| |B| + |S_j|, |B(mu + i j)| <= sum_k |x_k| 2 exp((j^2 - (mu - K h)^2)/(2 sigma^2)), |Q R_j| <= A mu^d
    (coefficient bounds) -> explicit log10 bound.
Exact cancellation used: 4z^2 + (3+2 sqrt2)^2 in Q against 4z^2 -+ 12 i z - 1 = (2z -+ i a)(2z -+ i b) in R_{-+3} (a, b = 3 -+ 2 sqrt2),
so f_{-+3} = Q0 Rc_{-+3} (2z -+ i b)/(2z +- i a) B, with Rc = R/(4z^2 -+ 12iz - 1) reduced exactly (sympy) -- no 0/0 near mu = +-0.086i.
Output: cache/hk3_epscert_<piece>.txt.   ~ 5-10 min per piece, one core.
Usage: python3 hk3_epscert.py c4v | c4d
"""
import os, sys, time, json, math
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'deshift'))
import sympy as sp
from flint import acb, arb, ctx
import kern                                     # read-only: exact R_j as Gaussian-rational coefficient lists

ctx.prec = 384                                  # the RBF sum cancels by ~1e-30 at |Im z| = 3: 192 bits is not enough
I = acb(0, 1)
PI = arb.pi()
nu_s = sp.Symbol('nu')
OUT = []
T0 = time.time()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


TINY = 5e-324                                   # smallest positive double: an upper bound for any underflowing positive arb


def up(x):
    """upper bound (float) of a nonnegative arb; never rounds a positive enclosure down to 0 (float underflow beyond ~1e-308)."""
    if not hasattr(x, 'upper'):
        return float(x)
    v = float(x.upper())
    return v if (v > 0 or x.upper() <= 0) else TINY


class Weight:
    def __init__(self, piece):
        d = json.load(open(os.path.join(ADS4, 'prl', 'data', f'q_deshift_continuum_{piece}.json')))
        assert d['piece'] == piece and d['h'] == 0.15 and d['sigma'] == 0.45
        self.piece = piece
        self.p = 1 if piece == 'c4v' else 3
        self.Phi = (arb(1), arb(0)) if piece == 'c4v' else (arb(0), arb('-0.5'))
        self.h = arb(d['h']); self.sig = arb(d['sigma']); self.K = d['K']     # the IEEE doubles, exactly as ab_solve built them
        self.s2 = 2*self.sig**2
        lab = d['labels']; assert lab[0] == ['g', 0] and all(l[0] == 're' for l in lab[1:])
        self.x = [arb(s) for s in d['coefficients']]
        self.R = kern.load_C4()
        a = 3 + 2*arb(2).sqrt(); b = 3 - 2*arb(2).sqrt()
        self.a, self.b = a, b
        # exact reductions R_{-+3} / (4 nu^2 -+ 12 i nu - 1)
        self.Rc = {}
        for j, fac in ((-3, 4*nu_s**2 + 12*sp.I*nu_s - 1), (3, 4*nu_s**2 - 12*sp.I*nu_s - 1)):
            red = sp.cancel(sp.together(self.R[j].expr/fac))
            assert sp.simplify(red*fac - self.R[j].expr) == 0
            num, den = sp.fraction(red)
            self.Rc[j] = kern.Rat(red)
        self.absx = sum(abs(xi) for xi in self.x)

    def Q0(self, z):
        z2 = z*z
        return z2*(z2 + 1)*(z2 + 4)/((z2 + 9)**self.p*(4*z2 + 1)*(4*z2 + 25))

    def B(self, z):
        g0 = (-(z*z)/self.s2).exp()
        e = (z*self.h/(self.sig**2)).exp(); ei = 1/e
        tot = self.x[0]*g0
        ep, em = acb(1), acb(1)
        for k in range(1, self.K + 1):
            ep *= e; em *= ei
            dk = (-(self.h*k)**2/self.s2).exp()
            tot += self.x[k]*g0*dk*(ep + em)
        return tot

    def B_box(self, z):
        """B on a box: termwise exp((-(z - c)^2)/s2) (no shared factors: avoids the dependency blow-up of g0 * e^k)."""
        tot = acb(0)
        for k in range(0, self.K + 1):
            c = self.h*k
            t = (-((z - c)**2)/self.s2).exp()
            if k > 0:
                t += (-((z + c)**2)/self.s2).exp()
            tot += self.x[k]*t
        return tot

    def f(self, j, mu, box=False):
        z = mu + I*j
        Bv = self.B_box(z) if box else self.B(z)
        if abs(j) == 3:
            s = 1 if j < 0 else -1          # j = -3: (2z + i b)/(2z - i a);  j = +3: (2z - i b)/(2z + i a)
            return self.Q0(z)*self.Rc[j](z)*(2*z + s*I*self.b)/(2*z - s*I*self.a)*Bv
        z2 = z*z
        Q = self.Q0(z)/(4*z2 + self.a**2)
        return Q*self.R[j](z)*Bv

    def S(self, j, mu):
        a = self.LA[j]; e = (-(mu*mu)/2).exp()
        if self.p == 1:
            return a[0]*e/mu
        return e*(a[2]/mu**3 + a[1]/mu**2 + (a[0] + a[2]/2)/mu)

    def E(self, mu, box=False):
        tot = acb(0)
        for j in (-3, -1, 1, 3):
            v = self.f(j, mu, box)
            if abs(j) == 3:
                v -= self.S(j, mu)
            tot += v
        return tot

    def Eeven(self, z, box=False):
        return (self.E(z, box) + self.E(-z, box))/2


def circle_max(fun, m, rho, nb=48):
    """enclosure (upper bound, float) of max |fun| on the circle |z - m| = rho: nb boxes covering the arcs."""
    hs = arb(rho)*math.pi/nb*1.02
    mx = 0.0
    for k in range(nb):
        th = 2*math.pi*(k + 0.5)/nb
        cx = arb(m) + arb(rho)*arb(math.cos(th)); cy = arb(rho)*arb(math.sin(th))
        Z = acb(arb(cx.mid(), hs), arb(cy.mid(), hs))
        v = abs(fun(Z))
        u = float(v.upper())
        if not math.isfinite(u):
            return float('inf')
        mx = max(mx, u)
    return mx


def taylor(fun, m, r, K):
    """c^_n, n = 0..K-1 (acb balls) from the K-point trapezoid on |z - m| = r."""
    vals = []
    for k in range(K):
        w = (2*PI*I*k/K).exp()
        vals.append((fun(m + r*w), w))
    out = []
    for n in range(K):
        acc = acb(0)
        for v, w in vals:
            acc += v*w**(-n)
        out.append(acc/(K*arb(r)**n))
    return out


def laurent(W):
    """ball-valued Laurent coefficients a_1..a_p of f_{+-3} at mu = 0 (a_k = Taylor coeff of order p - k of g = mu^p f_j)."""
    r0, rho0, K0 = arb('0.2'), 0.4, 160
    W.LA = {}
    for j in (-3, 3):
        g = lambda mu, j=j: mu**W.p*W.f(j, mu)
        gb = lambda Z, j=j: Z**W.p*W.f(j, Z, box=True)
        M = circle_max(gb, 0.0, rho0, nb=96)
        cs = taylor(g, acb(0), r0, K0)
        q = (0.2/rho0)**K0
        a = []
        for k in range(1, W.p + 1):
            n = W.p - k
            err = M*rho0**(-n)*q/(1 - q)
            c = cs[n]
            a.append(acb(arb(c.real.mid(), c.real.rad() + err), arb(c.imag.mid(), c.imag.rad() + err)))
        W.LA[j] = a
        pr(f'    j = {j:+d}: max|mu^p f_j| on |mu| = {rho0} <= {M:.3e};  aliasing bound {M*q/(1-q):.1e};  a_k = '
           + ', '.join(f'{v.real.str(12)} + i {v.imag.str(3)}' for v in a))
    s = {-3: -1, 3: 1}
    const = sum(I*PI*s[j]/2*W.LA[j][0] for j in (-3, 3))
    t2 = sum(I*PI*s[j]/2*(-W.LA[j][2]/2) for j in (-3, 3)) if W.p == 3 else acb(0)
    d0 = abs(const - W.Phi[0]); d2 = abs(t2 - W.Phi[1])
    a2 = abs(W.LA[-3][1] + W.LA[3][1]) if W.p == 3 else arb(0)
    pr(f'    residue part: const = {const.real.str(15)} (+ i {const.imag.str(3)}), target {W.Phi[0].str(3)} -> |const - Phi(0)| <= {up(d0):.2e}')
    if W.p == 3:
        pr(f'    t^2 coefficient = {t2.real.str(15)} (+ i {t2.imag.str(3)}), target -1/2 -> |dev| <= {up(d2):.2e};  |a2(-3) + a2(+3)| <= {up(a2):.2e}')
    return up(d0), up(d2), up(a2)


def rect_max(fun, x0, x1, y, bs=0.03):
    """enclosure of max |fun| on the boundary of [x0, x1] x [-y, y] (boxes of side bs); by the maximum-modulus principle this bounds
    |fun| on the whole rectangle when fun is analytic there (E_even is, for y < 1/2)."""
    pts = []
    nx = int(math.ceil((x1 - x0)/bs)); ny = int(math.ceil(2*y/bs))
    for k in range(nx):
        cx = x0 + (k + 0.5)*(x1 - x0)/nx
        pts += [(cx, -y), (cx, y)]
    for k in range(ny):
        cy = -y + (k + 0.5)*2*y/ny
        pts += [(x0, cy), (x1, cy)]
    hx = arb((x1 - x0)/nx/2*1.01 + 1e-12); hy = arb(2*y/ny/2*1.01 + 1e-12)
    mx = 0.0
    for (cx, cy) in pts:
        Z = acb(arb(cx, hx), arb(cy, hy))
        u = up(abs(fun(Z)))
        if not math.isfinite(u):
            return float('inf')
        mx = max(mx, u)
    return mx


def eps_cert(W, L=110.0, w=0.15, r=0.09, rho=0.3, chunk=2.0, target=1e-28):
    """panels of width w; Taylor data from a K-point trapezoid on |z - m| = r; Cauchy/maximum-modulus constant M from the boundary of
    the chunk rectangle [A - rho, B + rho] x [-rho, rho] (contains every circle |z - m| = rho, m in the chunk)."""
    tot = 0.0; peak = 0.0; Mmax = 0.0; mid_sum = 0.0; Kmax = 0
    hw = w/2
    assert hw <= r + 1e-15 and r < rho and hw < rho
    nch = int(round(L/chunk)); npc = int(round(chunk/w))
    for c in range(nch):
        A = c*chunk
        M = rect_max(lambda Z: W.Eeven(Z, box=True), A - rho - 0.02, A + chunk + rho + 0.02, rho + 0.02)
        if not math.isfinite(M):
            raise RuntimeError(f'rectangle enclosure failed at chunk {c}')
        Mmax = max(Mmax, M)
        q1 = r/rho; q2 = hw/rho
        K = max(24, int(math.ceil(math.log(target/max(M, 1.0))/math.log(max(q1, q2)))))
        Kmax = max(Kmax, K)
        N = K - 1
        q = q1**K
        for i in range(npc):
            m = A + (i + 0.5)*w
            cs = taylor(lambda z: W.Eeven(z), acb(m), arb(r), K)
            b = 0.0
            for n in range(N + 1):
                e_n = M*rho**(-n)*q/(1 - q)
                b += (up(abs(cs[n])) + e_n)*2*hw**(n + 1)/(n + 1)
            b += w*M*q2**(N + 1)/(1 - q2)
            b *= (1 + 1e-12)               # outward safety for the float accumulation
            tot += b
            mid_sum += abs(complex(cs[0].mid()))*w
            peak = max(peak, b/w)
        if c % 5 == 0:
            pr(f'      chunk {c}/{nch} [{A:.0f}, {A + chunk:.0f}]: M <= {M:.2e}, K = {K}, cumulative bound {tot:.4e}  [{time.time() - T0:.0f}s]')
    return tot, mid_sum, Mmax, peak, Kmax


def tail(W, L=110.0):
    """log10 of a rigorous bound on int_L^oo |E_even|."""
    Kh = float(W.h.upper())*W.K
    assert L > Kh + 10
    # |B(mu + i j)| <= 2 sum|x_k| exp((j^2 - (mu - Kh)^2)/s2) for mu >= L; |Q R_j(z)| <= A |z|^d for |z| >= L (coefficient bounds)
    best = -1e9
    for j in (-3, -1, 1, 3):
        num = [abs(c) for c in (W.Rc[j] if abs(j) == 3 else W.R[j])._acb()[0]]
        den = [abs(c) for c in (W.Rc[j] if abs(j) == 3 else W.R[j])._acb()[1]]
        dn, dd = len(num) - 1, len(den) - 1
        z = arb(L)
        numb = sum(c*z**(dn - k) for k, c in enumerate(num))/z**dn                     # |num(z)| <= numb |z|^dn
        denb = den[0] - sum(c*z**(-(k)) for k, c in enumerate(den) if k > 0)            # |den(z)| >= denb |z|^dd
        assert float(denb.lower()) > 0
        # Q (or Q0 x cancelled factor) <= 1 x |z|^{-2..} for |z| >= L: generous: |Q| <= 2 |z|^{6 - 2p - 6 - 2}
        A = numb/denb*4
        d = dn - dd + 6 - 2*W.p - 6 - 2 + 2
        # int_L^oo A mu^d 2 sum|x| exp((9 - (mu - Kh)^2)/s2) dmu <= A L^d 2 sum|x| e^{9/s2} e^{-(L-Kh)^2/s2} x (s2/(L-Kh)) x 2 (d small)
        s2 = float(W.s2.upper())
        lg = (math.log10(float(A.upper())) + d*math.log10(L) + math.log10(2*float(W.absx.upper())) + (9 - (L - Kh)**2)/s2/math.log(10)
              + math.log10(s2/(L - Kh)*2 + 1))
        best = max(best, lg)
    # S_j: |a| e^{-mu^2/2}/mu^k: e^{-L^2/2}
    lgS = math.log10(sum(up(abs(a)) for j in (-3, 3) for a in W.LA[j]) + 1e-300) - L*L/2/math.log(10)
    return max(best, lgS)


if __name__ == '__main__':
    piece = sys.argv[1]
    W = Weight(piece)
    pr(f'# hk3_epscert {piece}: rigorous (ball-arithmetic, {ctx.prec} bits) bound on eps_W for prl/data/q_deshift_continuum_{piece}.json')
    pr('(A) Laurent data of f_{-+3} at mu = 0 (ball-valued, aliasing included):')
    d0, d2, a2 = laurent(W)
    pr('(B) eps_W on [0, 110]:')
    epsB, midsum, Mmax, peak, Kmax = eps_cert(W)
    lt = tail(W)
    pr(f'    certified int_0^110 |E_even| <= {epsB:.4e}   (midpoint-rule estimate {midsum:.4e}; max panel average {peak:.2e}; '
       f'max rectangle enclosure M {Mmax:.2e}; trapezoid points up to K = {Kmax})')
    pr(f'(C) tail: int_110^oo |E_even| <= 10^{lt:.0f}')
    eps = epsB + 10**max(lt, -300)
    pr(f'RESULT {piece}: eps_W <= {eps:.4e} [certified];  |W(t) - Phi(t)| <= eps_W + {d0:.1e} + {d2:.1e} t^2 + {a2:.1e} (sqrt(2 pi) + pi t)/2 '
       f'for all t >= 0  [{time.time() - T0:.0f}s]')
    with open(os.path.join(HERE, 'cache', f'hk3_epscert_{piece}.txt'), 'w') as fh:
        fh.write('\n'.join(OUT) + '\n')
