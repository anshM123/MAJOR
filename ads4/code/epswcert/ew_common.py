"""
epswcert/ew_common.py (EPSW_CERT.md): exact data and Arb (python-flint) series machinery for the published k = 2 continuum weights
of Omega^EE (prl/data/k2/k2_continuum_D_<Dg>.json), read-only.

Objects (CERT_CLOSE §3.2, k2/README §2):
  chi(nu) = Q(nu) B(nu) + sum_{j<=2} y_j nu^{2j} l(nu) e^{-nu^2/2},  Q = nu^2/(nu^2+1/4)^2,  l = (nu^2+1) Log(nu^2+1),
  B(nu) = sum_k x_k [g_k(nu - c_k) + g_k(nu + c_k)]  (one term for c_k = 0),  g_k(z) = exp(-z^2/(2 sigma_k^2)).
  Coefficients x_k, y_j, centres c_k, widths sigma_k: the IEEE doubles of the JSON, taken EXACTLY (arb(float)).
  f(mu) := chi(mu - i) R_-(mu - i),   E(mu) := 2 Re f(mu),   r(mu) := E(mu) - psi(mu)  (mu >= 0),
  psi(nu) = psi2^NF02(nu/Lam)/k0 (nu <= Lam), 0 beyond;  Lam = sqrt((Dg-2)(Dg-1)) exactly;  k0 = 27.48509817 exactly.
  W[chi](t) - Phi(t) = int_0^oo r(mu) cos(mu t) dmu   (shifted form, valid for every real t).

Exact algebra used:  Q R_- = -32 nu (2nu+3i)(2nu+5i)(2nu+11i) / ((2nu-i)(2nu-7i)(2nu-3i)(2nu+7i))   (asserted with sympy),
so f = [Q R_-](nu) B(nu) + Y(nu) l(nu) R_-(nu),  Y(nu) = sum_j y_j nu^{2j} e^{-nu^2/2}.
Near mu = 0 (nu = mu - i):  nu^2 + 1 = mu (mu - 2i)  and, for mu > 0, Log(mu(mu - 2i)) = log mu + Log(mu - 2i), so
  f = F1(mu) + F2(mu) log mu,   F1 = [QR_-]B + mu(mu-2i)Log(mu-2i) Y R_-,   F2 = mu(mu-2i) Y R_-,
F1, F2 analytic in |mu| < 1 (nearest singularity: the pole of R_- at nu = 0, i.e. mu = i).
"""
import os, sys, json, math
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import flint
from flint import arb, acb, arb_poly, acb_series, ctx

ctx.prec = 192
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
CACHE = os.path.join(HERE, 'cache')
K0 = Fr('27.48509817')
I = acb(0, 1)


def up(x):
    """float upper bound of a real arb (never rounds a positive enclosure to 0)."""
    v = float(x.upper())
    if v <= 0 and x.upper() > 0:
        return 5e-324
    return v


def absup(x):
    """upper bound of |x| as an arb with zero radius (x arb or acb)."""
    if isinstance(x, acb):
        a = abs(x)
    else:
        a = abs(x)
    return arb(a.upper())


def amax(a, b):
    """exact max of two exact (zero-radius) arbs, as an upper bound."""
    a = arb(a.upper()); b = arb(b.upper())
    return a if a > b else b


def pr_to(lines):
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    return pr


# ---------------------------------------------------------------- exact kernel R_- (prl/data/kernels.json, D_k2_EE)
def _gauss_poly(factors, const):
    """expand const * prod (a nu + b i) with integer a, b into a list of (re, im) Fraction coefficients, lowest first."""
    poly = [(Fr(const), Fr(0))]
    for a, b in factors:                       # multiply by (a nu + i b)
        new = [(Fr(0), Fr(0))]*(len(poly) + 1)
        for k, (pr_, pi_) in enumerate(poly):
            # (pr + i pi) * (i b) -> degree k ; * a nu -> degree k+1
            r0, i0 = new[k]; new[k] = (r0 - pi_*b, i0 + pr_*b)
            r1, i1 = new[k + 1]; new[k + 1] = (r1 + pr_*a, i1 + pi_*a)
        poly = new
    return poly


# R_- = -2 (2nu - i)(2nu + i)^2 (2nu + 3i)(2nu + 5i)(2nu + 11i) / ( nu (2nu - 7i)(2nu - 3i)(2nu + 7i) )
RM_NUM = _gauss_poly([(2, -1), (2, 1), (2, 1), (2, 3), (2, 5), (2, 11)], -2)
RM_DEN = _gauss_poly([(1, 0), (2, -7), (2, -3), (2, 7)], 1)
# Q R_- = -32 nu (2nu+3i)(2nu+5i)(2nu+11i) / ((2nu - i)(2nu - 7i)(2nu - 3i)(2nu + 7i))
QR_NUM = _gauss_poly([(1, 0), (2, 3), (2, 5), (2, 11)], -32)
QR_DEN = _gauss_poly([(2, -1), (2, -7), (2, -3), (2, 7)], 1)


def check_kernels():
    """sympy: R_- equals prl/data/kernels.json's string, and Q R_- equals the cancelled form above (exact)."""
    import sympy as sp
    nu = sp.Symbol('nu')
    d = json.load(open(os.path.join(DATA, 'kernels.json')))['D_k2_EE']
    Rm = sp.sympify(d['R_minus'], locals={'nu': nu, 'I': sp.I})
    tolist = lambda P: sum((sp.Rational(a.numerator, a.denominator) + sp.I*sp.Rational(b.numerator, b.denominator))*nu**k
                           for k, (a, b) in enumerate(P))
    assert sp.simplify(tolist(RM_NUM)/tolist(RM_DEN) - Rm) == 0
    Q = nu**2/(nu**2 + sp.Rational(1, 4))**2
    assert sp.simplify(tolist(QR_NUM)/tolist(QR_DEN) - Q*Rm) == 0
    assert sp.simplify(Rm.subs(nu, -sp.I) - sp.Rational(18, 25)) == 0
    return True


def _acbc(P):
    return [acb(arb(a.numerator)/a.denominator, arb(b.numerator)/b.denominator) for (a, b) in P]


RMN, RMD, QRN, QRD = (_acbc(P) for P in (RM_NUM, RM_DEN, QR_NUM, QR_DEN))


def horner(coeffs, z):
    out = coeffs[-1]*1
    for c in reversed(coeffs[:-1]):
        out = out*z + c
    return out


# ---------------------------------------------------------------- the published weight
class Weight:
    def __init__(self, Dg):
        self.Dg = Dg
        d = json.load(open(os.path.join(DATA, 'k2', f'k2_continuum_D_{Dg}.json')))
        assert d['Delta_g'] == Dg
        self.json = d
        self.LAM2 = Fr((Dg - 2)*(Dg - 1))
        self.Lam = arb(self.LAM2.numerator).sqrt()
        cs = d['centres_sigma']; xs = [float(s) for s in d['x_re']]
        assert len(cs) == len(xs)
        self.y = [arb(float(s)) for s in d['y_re']]
        # signed centres: (c, sigma, x) with c = +c_k and (c_k > 0) -c_k
        sc = []
        for (c, s), x in zip(cs, xs):
            assert c >= 0 and s in (0.45, 0.9)
            sc.append((c, s, x))
            if c > 0:
                sc.append((-c, s, x))
        sc.sort()
        self.sc = sc
        self.sc_arb = [(arb(c), arb(s), arb(x), 1/(2*arb(s)**2)) for c, s, x in sc]
        self.cs_sorted = [c for c, s, x in sc]
        self.sumabs = {}
        for c, s, x in sc:
            self.sumabs[s] = self.sumabs.get(s, arb(0)) + abs(arb(x))
        self.cmax = max(c for c, s, x in sc)
        nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
        assert Fr(nf['k0']) == K0 if isinstance(nf['k0'], str) else abs(nf['k0'] - 27.48509817) < 1e-15
        self.psi2 = [Fr(c) for c in nf['psi2_exact']]           # psi2^NF02(p), unnormalised; psi = psi2(nu/Lam)/k0
        # psi in mu: coefficient psi2_n / (k0 Lam^n)
        self.psi_mu = [arb(c.numerator)/c.denominator/(arb(K0.numerator)/K0.denominator)/self.Lam**n
                       for n, c in enumerate(self.psi2)]
        self.kappa = self.psi_mu[1]
        assert self.psi2[0] == 0

    # ---- psi(m + s) Taylor coefficients (exact polynomial; valid only for m + s <= Lam)
    def psi_taylor(self, m, N):
        P = arb_poly(self.psi_mu)
        q = P(arb_poly([m, 1]))
        c = [q[n] if n < q.degree() + 1 else arb(0) for n in range(N)]
        return c

    # ---- near / far split of the RBF sum around the real point m
    CUT = 15.0                                  # RBF terms with |m - c| > CUT*sigma are "far" (bounded, not expanded)

    def _near_idx(self, m):
        import bisect
        lo_ = bisect.bisect_left(self.cs_sorted, m - self.CUT*0.9 - 1e-6)
        hi_ = bisect.bisect_right(self.cs_sorted, m + self.CUT*0.9 + 1e-6)
        return [i for i in range(lo_, hi_) if abs(m - self.sc[i][0]) <= self.CUT*self.sc[i][1]]

    def far_bound(self, rho):
        """for |z - m| <= rho: sum over far terms of |x_k g_k(z - i - c_k)| <= sum_sigma exp(((1+rho)^2 - (CUT sigma - rho - 1e-6)^2)/(2 sigma^2))
        * sum_{k: sigma_k = sigma} |x_k|   (|g(w)| = exp(-Re(w^2)/(2 s^2)), Re(w^2) >= (|Re w| )^2 - (Im w)^2, |Im w| <= 1 + rho)."""
        tot = arb(0)
        for s, sx in self.sumabs.items():
            sa = arb(s)
            tot += sx*((arb(1 + rho)**2 - (self.CUT*sa - rho - arb('1e-6'))**2)/(2*sa*sa)).exp()
        return tot

    def near_far(self, m, rho):
        idx = self._near_idx(m)
        near = [(self.sc_arb[i][0], self.sc_arb[i][2], self.sc_arb[i][3]) for i in idx]
        return near, self.far_bound(rho)

    def B_disc_bound(self, m, rho):
        """sup over |z - m| <= rho of |B(z - i)|: near terms one by one (arb), far terms by far_bound."""
        tot = self.far_bound(rho)
        ma = arb(m)
        for i in self._near_idx(m):
            ca, sa, xa, inv2s2 = self.sc_arb[i]
            d = abs(ma - ca) - rho
            d2 = arb(0) if not d > 0 else d*d
            tot += abs(xa)*((arb(1 + rho)**2 - d2)*inv2s2).exp()
        return tot


def series_var(z0, N):
    return acb_series([z0, 1], prec=N)


def rat_series(num, den, S):
    return horner(num, S)/horner(den, S)


def gauss_sum_series(near, S):
    """sum_k x_k exp(-(S - c_k)^2/(2 s_k^2)) as acb_series."""
    out = None
    for ca, xa, inv2s2 in near:
        D = S - ca
        t = (-(D*D)*inv2s2).exp()*xa
        out = t if out is None else out + t
    return out


def Y_series(W, S):
    S2 = S*S
    e = (-S2/2).exp()
    return (W.y[0] + W.y[1]*S2 + W.y[2]*S2*S2)*e


def f_series(W, m, N):
    """Taylor coefficients (acb, length N) at s = 0 of f_near(m + s) for real m (arb ball allowed) with m >= 0.2: the near RBF
    terms (|m - c_k| <= CUT sigma_k) and the y-terms, principal Log (analytic for Re >= 0.2 > 0, see module doc)."""
    nu0 = acb(m, -1)
    S = series_var(nu0, N)
    near, _ = W.near_far(float(m.mid()) if isinstance(m, arb) else m, 0.5)
    QR = rat_series(QRN, QRD, S)
    B = gauss_sum_series(near, S)
    L1 = S*S + 1
    l = L1*L1.log()
    Rm = rat_series(RMN, RMD, S)
    f = QR*B + Y_series(W, S)*l*Rm
    return [f[n] for n in range(N)]


def disc_box(m, rho):
    """acb box containing {z - i : |z - m| <= rho} (m float, rho float)."""
    return acb(arb(m, rho), arb(-1, rho))


def qr_sup(m, rho):
    box = disc_box(m, rho)
    return abs(horner(QRN, box)/horner(QRD, box)).upper()


def f_disc_sup(W, m, rho):
    """(M, F): M >= sup_{|z-m|<=rho} |f(z)| (explicit RBF bound times box enclosure of |Q R_-|, plus box enclosure of |Y l R_-|);
    F >= sup of the far-RBF part |Q R_- B_far| on the same disc.  Requires m - rho > 0 (Log analytic)."""
    assert m - rho > 0.05
    box = disc_box(m, rho)
    qr = abs(horner(QRN, box)/horner(QRD, box)).upper()
    nu2 = box*box
    ybox = (W.y[0] + W.y[1]*nu2 + W.y[2]*nu2*nu2)*((-nu2/2).exp())
    L1 = nu2 + 1
    yl = abs(ybox*L1*L1.log()*horner(RMN, box)/horner(RMD, box))
    return (qr*W.B_disc_bound(m, rho)).upper() + yl.upper(), (qr*W.far_bound(rho)).upper()


# ---------------------------------------------------------------- near mu = 0: f = F1 + F2 log mu
def F12_series(W, N):
    S = acb_series([0, 1], prec=N)                      # mu
    nu = S - I
    QR = rat_series(QRN, QRD, nu)
    near = [(ca, xa, inv2s2) for (c, s, x), (ca, sa, xa, inv2s2) in zip(W.sc, W.sc_arb)]
    B = gauss_sum_series(near, nu)                     # all terms (exact; no far part)
    Rm = rat_series(RMN, RMD, nu)
    YR = Y_series(W, nu)*Rm
    mm = S*(S - 2*I)
    F2 = mm*YR
    F1 = QR*B + mm*(S - 2*I).log()*YR
    return [F1[n] for n in range(N)], [F2[n] for n in range(N)]


def F12_disc_sup(W, rho):
    """sup over |mu| <= rho (rho < 1) of |F1|, |F2| (box enclosures on a ring of boxes covering the circle |mu| = rho;
    maximum modulus principle)."""
    nb = 256
    s1 = arb(0); s2 = arb(0)
    two_pi = 2*arb.pi()
    for k in range(nb):
        th = two_pi*(arb(k) + arb(0.5, 0.5))/nb
        z = acb(rho*th.cos(), rho*th.sin())
        nu = z - I
        qr = horner(QRN, nu)/horner(QRD, nu)
        nu2 = nu*nu
        Y = (W.y[0] + W.y[1]*nu2 + W.y[2]*nu2*nu2)*((-nu2/2).exp())
        YR = Y*horner(RMN, nu)/horner(RMD, nu)
        mm = z*(z - 2*I)
        Bb = arb(0)
        for (c, s, x), (ca, sa, xa, inv2s2) in zip(W.sc, W.sc_arb):
            Bb += abs(xa)*abs((-((nu - ca)**2)*inv2s2).exp())
        a1 = abs(qr)*Bb + abs(mm*(z - 2*I).log()*YR)
        a2 = abs(mm*YR)
        s1 = amax(s1, a1.upper()); s2 = amax(s2, a2.upper())
    return s1, s2
