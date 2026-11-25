"""
fullpos3/fp3_htdata.py (FULLPOS3.md §4): certified static inputs of the high-twist (h-expansion) bound.  Arb (python-flint) balls;
every function returns an UPPER bound (or a ball containing the value) unless stated.

Read-only inputs: prl/data/nf02_psi.json (psi2, psi3), prl/data/k2/k2_continuum_D_1000_rep.json (realised weight, only for the
Gaussian tail bound beyond mu = 1059), code/epswcert/cache/ew_cert_1000_re_rep.txt (certified cumulative int|r| on mu-panels),
code/lrtight/cache/lrt3_rcert.json (certified F_cert(beta) cells).

Weights psi_k (k = 2, 3): NF02 exact polynomial / k0 on [0,1]; psi_k(0) = 0 and (1-p)^5 | psi_k (asserted).
  A_s^k(x) = sum_{j even >= 2} |psi_k^(j)(0)| / (Lam x)^{j+1}          (p = 0 endpoint, sine transform)
  A_e^k(x) = sum_{j >= 5}      |psi_k^(j)(1)| / (Lam x)^{j+1}          (p = 1 endpoint, carries e^{Lam |Im t|})
  C_cor^k(eta) >= int_0^1 |psi_k(p)| 2 cosh(Lam p eta)/(e^{2 pi Lam p} - 1) dp   (coth correction)
               <= sum_n |a_n| Lam^{-n-1} 2 n! sum_{m>=1} (2 pi m - eta)^{-n-1}  (|psi| <= sum |a_n| p^n, extend to infinity)
  J_X(beta) >= int_0^oo A_X(beta+u) (2 sinh(beta+u/2) sinh(u/2))^{-1/2} du   (X = s, e) and J_C(beta) (A = 1): decreasing in beta.
Lemma (FULLPOS3 §4.2): for Re xi > 0, |Im xi| < pi/2,
  |int_0^1 psi_k(p) P_{-1/2+i Lam p}(cosh xi) dp| <= (sqrt2/pi)[J_s(beta) + C_cor(eta) J_C(beta) + e^{Lam eta} J_e(beta)],
  beta = Re xi, eta = |Im xi|.
"""
import os, sys, json, math, re
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
from flint import arb, acb, ctx, fmpq

ctx.prec = 256
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
DATA = os.path.join(ADS4, 'prl', 'data')
K0 = Fr('27.48509817')
LAM2 = Fr(998*999)


def A(q):
    q = Fr(q)
    return arb(fmpq(q.numerator, q.denominator))


LAM = A(LAM2).sqrt()
PI = arb.pi()


def up(x):
    return arb(x.upper())


def lo(x):
    return arb(x.lower())


def psi_coeffs(k):
    d = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
    k0 = Fr(d['k0'])
    assert k0 == K0
    return [Fr(c)/k0 for c in d['psi%d_exact' % k]]


def derivs(c, p0):
    D = len(c) - 1; out = []
    for j in range(D + 1):
        s = Fr(0)
        for n in range(j, D + 1):
            ff = 1
            for i in range(j):
                ff *= (n - i)
            s += c[n]*ff*Fr(p0)**(n - j)
        out.append(s)
    return out


class Psi:
    def __init__(self, k):
        self.k = k
        self.c = psi_coeffs(k)
        assert self.c[0] == 0
        self.d0 = derivs(self.c, 0); self.d1 = derivs(self.c, 1)
        assert all(self.d1[j] == 0 for j in range(5))
        self.sterms = [(j + 1, A(abs(self.d0[j]))) for j in range(2, len(self.c), 2) if self.d0[j] != 0]
        self.eterms = [(j + 1, A(abs(self.d1[j]))) for j in range(5, len(self.c)) if self.d1[j] != 0]
        self.absc = [A(abs(a)) for a in self.c]

    def As(self, x):            # x arb > 0 (use lower bound of x for an upper bound: decreasing)
        X = LAM*x
        return sum((a/X**p for p, a in self.sterms), arb(0))

    def Ae(self, x):
        X = LAM*x
        return sum((a/X**p for p, a in self.eterms), arb(0))

    def Ccor(self, eta):
        """upper bound of int_0^1 |psi| 2 cosh(Lam p eta)/(e^{2 pi Lam p}-1) dp, 0 <= eta < pi."""
        eta = arb(eta)
        assert eta < PI
        tot = arb(0)
        for n, a in enumerate(self.absc):
            if n == 0 or a == 0:
                continue
            # sum_{m>=1} (2 pi m - eta)^{-n-1} <= (2pi-eta)^{-n-1} + int_1^oo (2 pi x - eta)^{-n-1} dx
            S = (2*PI - eta)**(-n - 1) + (2*PI - eta)**(-n)/(2*PI*n)
            tot += a*LAM**(-n - 1)*2*arb.fac_ui(n)*S
        return up(tot)

    def absexp(self, kappa, ncell=2000):
        """upper bound of int_0^1 |psi(p)| e^{kappa Lam p} dp (kappa >= 0): |psi| <= max on cells (ball evaluation)."""
        tot = arb(0); kap = arb(kappa)
        for i in range(ncell):
            pa, pb = Fr(i, ncell), Fr(i + 1, ncell)
            pball = A(pa).union(A(pb))
            v = sum((A(a)*pball**n for n, a in enumerate(self.c)), arb(0))
            m = abs(v).upper()
            # int_pa^pb e^{kappa Lam p} dp
            if kappa == 0:
                I = A(pb - pa)
            else:
                I = ((kap*LAM*A(pb)).exp() - (kap*LAM*A(pa)).exp())/(kap*LAM)
            tot += arb(m)*I
        return up(tot)


def jtable(psi, betas, u0=Fr(1, 10**6), U=60, ng=600):
    """upper bounds of J_s, J_C, J_e at each beta in the grid (list of Fractions)."""
    out = []
    # geometric u-grid on [u0, U]
    us = [A(u0)*(A(U)/A(u0))**(arb(i)/ng) for i in range(ng + 1)]
    for b in betas:
        b = A(b)
        def kern(u):        # (2 sinh(b+u/2) sinh(u/2))^{-1/2}, decreasing in u
            return 1/(2*(b + u/2).sinh()*(u/2).sinh()).sqrt()
        Js = Jc = Je = arb(0)
        # [0, u0]: sinh(b+u/2) >= sinh b, sinh(u/2) >= u/2  ->  int <= A(b) (sinh b)^{-1/2} * 2 sqrt(u0)
        head = 2*us[0].sqrt()/b.sinh().sqrt()
        Js += psi.As(b)*head; Jc += head; Je += psi.Ae(b)*head
        for i in range(ng):
            ua, ub = us[i], us[i + 1]
            kk = kern(ua)*(ub - ua)
            Js += psi.As(b + ua)*kk; Jc += kk; Je += psi.Ae(b + ua)*kk
        # tail u >= U: kernel <= sqrt2 e^{-(b+u)/2}/(1-e^{-U})  ->  int_U^oo <= 2 sqrt2 e^{-(b+U)/2}/(1-e^{-U})
        UU = arb(U)
        tail = 2*arb(2).sqrt()*(-(b + UU)/2).exp()/(1 - (-UU).exp())
        Js += psi.As(b + UU)*tail; Jc += tail; Je += psi.Ae(b + UU)*tail
        out.append((up(Js), up(Jc), up(Je)))
    return out


# ------------------------------------------------------------------ realised-weight residual r = f_theta - psi(nu/Lam)
def eps_cells():
    """certified cells (mu_lo, mu_hi, upper bound of int_cell |r|) from the EPSW log (cumulative sums of panel bounds; each printed
    cumulative value rounded to 5 significant digits -> +2e-10 per cell; cell edges widened by one panel (0.1) for misalignment)."""
    txt = open(os.path.join(ADS4, 'code', 'epswcert', 'cache', 'ew_cert_1000_re_rep.txt')).read()
    m0 = re.search(r'\[0, 0\.2\]: int\|r\| <= ([0-9.e+-]+)', txt)
    pts = [(0.2, float(m0.group(1)))]
    for mm in re.finditer(r'panel (\d+)/\d+ at mu = ([0-9.]+): int\|r\| so far ([0-9.e+-]+)', txt):
        pts.append((float(mm.group(2)), float(mm.group(3))))
    epsW = float(re.search(r'eps_W = int\|r\| <= ([0-9.e+-]+)', txt).group(1))
    cells = [(0.0, 0.2 + 0.1, pts[0][1] + 2e-10)]
    for (ma, ca), (mb, cb) in zip(pts[:-1], pts[1:]):
        cells.append((ma - 0.1, mb + 0.1, max(cb - ca, 0) + 2e-10))
    cells.append((pts[-1][0] - 0.1, 1059.0, max(epsW - pts[-1][1], 0) + 2e-10))
    return cells, epsW


def gauss_tail_params():
    d = json.load(open(os.path.join(ADS4, 'prl', 'data', 'k2', 'k2_continuum_D_1000_rep.json')))
    cs = d['centres_sigma']; xs = d['x_re']
    S = sum(abs(Fr(float(x))) for x in xs)
    cmax = max(Fr(float(c)) for c, s in cs); smin = min(Fr(float(s)) for c, s in cs); smax = max(Fr(float(s)) for c, s in cs)
    return S, cmax, smin, smax


def Rint(eta):
    """upper bound of int_0^oo |r_theta(nu)| e^{eta nu} dnu (eta >= 0), r_theta = f_theta - psi(nu/Lam).
    Cells from the EPSW certificate (r without theta) + theta correction + Gaussian tail beyond 1059."""
    eta = arb(eta)
    cells, epsW = eps_cells()
    tot = arb(0)
    for (a, b, w) in cells:
        tot += arb(w)*(eta*arb(b)).exp()
    # tail nu >= 1059: |f| <= |Q R_-|(nu - i) |B(nu - i)| <= 16 * 2 S e^{1/(2 smin^2)} e^{-(nu - cmax)^2/(2 smax^2)}
    S, cmax, smin, smax = gauss_tail_params()
    pref = 16*2*A(S)*(1/(2*A(smin)**2)).exp()
    # int_{1059}^oo e^{eta nu - (nu-cmax)^2/(2 s^2)} dnu <= sqrt(2 pi) s e^{eta cmax + eta^2 s^2/2} * [complete Gaussian], and
    # since 1059 - cmax - eta s^2 >= 30 - 3 > 0 we bound by the Gaussian tail factor e^{-(1059-cmax-eta s^2)^2/(2 s^2)}
    s = A(smax); x0 = arb(1059) - A(cmax) - eta*s*s
    assert x0 > 20
    tail = pref*(2*PI).sqrt()*s*(eta*A(cmax) + eta*eta*s*s/2).exp()*(-(x0*x0)/(2*s*s)).exp()
    # theta correction: |f_theta - f| <= |1 - theta(nu - i)| |f|; for nu <= Lam+30, |1-theta| <= erfc-type <= e^{-40};
    # int|f| <= Lam int|psi2| + epsW <= 1.1e4;  for nu > Lam + 30 the realised |f| <= int_{>Lam}|r| (EPSW range > Lam) <= 3.2e-10,
    # with |1 - theta| <= 2.1.  Weighted: e^{eta (Lam+130)} for the Gaussian-supported part (centres <= 1028.7, theta cut at 1098.5).
    # theta correction (r_theta - r = f_theta - f = 2 Re[chi R_- (theta - 1)](nu - i)):
    #  nu <= Lam + 30: |1 - theta(nu - i)| <= e^{0.01} erfc(7)/... <= e^{-49}  (theta: centre Lam+100, width 10), and
    #     int 2|chi R_-|(nu - i) dnu <= 1e10 (crude: 32 S e^{1/(2 smin^2)} sqrt(2 pi) smax ~ 1e9 for the RBF part, plus log terms);
    #  nu > Lam: r = f (psi = 0) and int_{>Lam} |r| <= 3.2e-10 (EPSW range > Lam), |theta| <= 1.1 -> |f_theta - f| <= 2.1|f|.
    thc = arb(1e10)*arb(-49).exp()*(eta*(LAM + 30)).exp() + arb(2.1)*arb(3.2e-10)*(eta*arb(1059)).exp()
    return up(2.1*tail + tot + thc), epsW


def P_mhalf(beta):
    """P_{-1/2}(cosh beta) = 2F1(1/2, 1/2; 1; -sinh^2(beta/2)) (ball)."""
    b = arb(beta)
    z = -(b/2).sinh()**2
    return z.hypgeom_2f1(arb(1)/2, arb(1)/2, arb(1))


def fcert_cells():
    d = json.load(open(os.path.join(ADS4, 'code', 'lrtight', 'cache', 'lrt3_rcert.json')))
    return [(Fr(r['ba']), Fr(r['bb']), r['Fcert']) for r in d['rows']]
