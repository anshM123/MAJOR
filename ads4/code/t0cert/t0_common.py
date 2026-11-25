"""
t0cert/t0_common.py (T0_CERT.md): shared exact / Arb helpers for the certificate of Theorem T0 (leading-Regge hybrid action of
Omega^EE >= c N_s on every heavy state, every Delta_g >= DGMIN).

Read-only imports (nothing of theirs is written; no bytecode):
  lrtight/lrt_common.py   exact psi_k (prl/data/nf02_psi.json, /k0), Arb helpers;
  gaplarge/gl_common.py   HatWeight W[k] (exact psi^(j)(0), psi^(j)(1), a1, certified int|psi|), m1_cells() (certified M1 lower
                          bounds on beta-cells of [0.1, 64]), build_I / E_up / tail_majorant (Lemma LR-inf pieces, Lam-parametrised).
Everything here is exact rational or Arb ball arithmetic.

Objects (k = 2, 3; psi = psi_k; p in [0, 1]):
  chi   = 2 psi - p psi'           (exact polynomial; chi(1) = chi'(1) = 0, chi'(0) = a1, chi''(0) = 0)
  F_k(b) = int_0^1 chi''(p) J0(b p) dp          (gap-unit curvature: Lam^2 curv_k = a1 M1(beta) + F_k(bhat)/24 + eps h, T0_CERT §2)
  hbar_k(x): the Taylor-remainder majorant of T0_CERT Lemma 2 (x <= 2 Lam), H_k = sup_x hbar_k(x).
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
CODE = os.path.abspath(os.path.join(HERE, '..'))
sys.path.insert(0, os.path.join(CODE, 'gaplarge')); sys.path.insert(0, os.path.join(CODE, 'lrtight'))
from fractions import Fraction as Fr
from math import factorial
import flint
from flint import arb, ctx, fmpq, fmpq_poly
import lrt_common as LC                    # read-only
import gl_common as G                      # read-only

ctx.prec = 128
CACHE = os.path.join(HERE, 'cache')
A, up, lo = LC.A, LC.up, LC.lo
PI = arb.pi()

DGMIN = Fr(os.environ.get('T0_DGMIN', '1000'))
LAM2MIN = (DGMIN - 2)*(DGMIN - 1)            # Lam^2 = (Dg-2)(Dg-1), exact; Lam increasing in Dg
EPS0 = 1/A(LAM2MIN)                           # eps = 1/Lam^2 <= EPS0 for every Dg >= DGMIN
LAMMIN = A(LAM2MIN).sqrt()
C = Fr(os.environ.get('T0_C', '199/1000'))    # the constant c of T0
B_POS = Fr(os.environ.get('T0_BPOS', '16'))   # split: bhat <= B_POS flat-dominated (Part F), bhat >= B_POS curvature-positive (Part K)


def pr_to(lines):
    def pr(s=''):
        print(s, flush=True); lines.append(s)
    return pr


# ------------------------------------------------------------------ exact polynomials
def pder(c):
    return [i*c[i] for i in range(1, len(c))] or [Fr(0)]


def pval(c, x):
    r = Fr(0)
    for a in reversed(c):
        r = r*x + a
    return r


PSI = {k: LC.psi_coeffs(k) for k in (2, 3)}


def chi_of(psi):
    d = pder(psi); out = [2*x for x in psi]
    for i, x in enumerate(d):
        out[i + 1] -= x                        # - p psi'
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


CHI = {k: chi_of(PSI[k]) for k in (2, 3)}
CHI2 = {k: pder(pder(CHI[k])) for k in (2, 3)}
A1 = {k: PSI[k][1] for k in (2, 3)}
for k in (2, 3):
    ch, d1 = CHI[k], pder(CHI[k])
    assert pval(ch, 1) == 0 and pval(d1, 1) == 0 and d1[0] == A1[k] and CHI2[k][0] == 0
    assert pval(CHI2[k], 1) == 0 and pval(pder(CHI2[k]), 1) == 0      # chi'' has (1-p)^2


def _fq(x):
    return fmpq(x.numerator, x.denominator)


def taylor_shift(c, x0):
    """exact coefficients of P(x0 + t) in t."""
    P = fmpq_poly([_fq(a) for a in c])
    S = P(fmpq_poly([_fq(Fr(x0)), 1]))
    n = S.degree()
    return [Fr(int(S.coeffs()[i].p), int(S.coeffs()[i].q)) for i in range(n + 1)] if n >= 0 else [Fr(0)]


def abs_sup_cell(c, a, b):
    """rigorous upper bound of sup_{[a,b]} |P| (exact Taylor shift to the centre; |P(m+t)| <= sum |c_j| h^j)."""
    m = (Fr(a) + Fr(b))/2; h = (Fr(b) - Fr(a))/2
    S = taylor_shift(c, m)
    return sum(abs(x)*h**j for j, x in enumerate(S))


def abs_int_upper(c, ncell=2000, wpow=Fr(0)):
    """rigorous upper bound of int_0^1 |P(p)| p^wpow dp (wpow >= 0 rational; p^wpow <= b^wpow on the cell) or with wpow < 0
    (p^wpow integrable, decreasing): sup|P| * int_cell p^wpow exactly."""
    tot = arb(0)
    for i in range(ncell):
        a, b = Fr(i, ncell), Fr(i + 1, ncell)
        s = A(abs_sup_cell(c, a, b))
        if wpow >= 0:
            tot += s*A(b)**A(wpow)*A(b - a)
        else:
            e = A(wpow + 1)
            tot += s*(A(b)**e - (A(a)**e if a > 0 else arb(0)))/e
    return up(tot)


# ------------------------------------------------------------------ F_k(b) = int chi'' J0(b p) dp (power series, exact moments)
MOM = {}


def moments(k, M):
    if (k, M) not in MOM:
        c = CHI2[k]
        MOM[(k, M)] = [sum(ci/(i + 2*m + 1) for i, ci in enumerate(c)) for m in range(M + 1)]
    return MOM[(k, M)]


L0 = {}; L1 = {}; L2 = {}


def init_L():
    for k in (2, 3):
        if k not in L0:
            L0[k] = abs_int_upper(CHI2[k], 2000, Fr(0))
            L1[k] = abs_int_upper(CHI2[k], 2000, Fr(1))
            L2[k] = abs_int_upper(CHI2[k], 2000, Fr(2))


def F_point(k, b):
    """enclosures of F_k(b) and F_k'(b) at a rational point b >= 0 (power series of J0, tail bounded by |mu_{2m}| <= L0)."""
    init_L()
    b = Fr(b)
    old = ctx.prec
    ctx.prec = 128 + int(2*float(b)) + 32
    try:
        M = int(2*float(b)) + 60
        mu = moments(k, M + 2)
        X = A(b)**2/4
        F = arb(0); dF = arb(0); term = arb(1)          # term = X^m/(m!)^2
        for m in range(M + 1):
            if m > 0:
                term = term*X/(m*m)
            mum = A(mu[m])
            F += (term*mum if m % 2 == 0 else -term*mum)
            if m > 0:
                dt = term*2*m/A(b)                       # d/db X^m/(m!)^2 = 2m X^m/(b (m!)^2)
                dF += (dt*mum if m % 2 == 0 else -dt*mum)
        # tail m > M: ratio X/(m+1)^2 <= 1/4 since M+1 >= 2b >= 2 sqrt(X)*... (X/(M+1)^2 = b^2/(4(M+1)^2) <= 1/16)
        t1 = term*X/((M + 1)**2)
        tailF = up(t1*L0[k]*4/3)
        tailD = up(t1*L0[k]*4/3*2*(M + 2)/A(max(b, Fr(1, 10**6))))
        F = F + arb(0, tailF.upper()); dF = dF + arb(0, tailD.upper())
    finally:
        ctx.prec = old
    return F, dF


def F_cell(k, ba, bb):
    """enclosure of F_k on [ba, bb]: F(bc) + F'(bc)[-h, h] + [-1, 1] L2 h^2/2   (|F''| <= int |chi''| p^2 since |J0''| <= 1)."""
    init_L()
    bc = (Fr(ba) + Fr(bb))/2; h = (Fr(bb) - Fr(ba))/2
    F, dF = F_point(k, bc)
    hh = A(h)
    return F + dF*arb(0, hh.upper()) + arb(0, (L2[k]*hh*hh/2).upper())


# tail: chi'' = p f, f(1) = f'(1) = 0 ;  F(b) = int f J0(bp) p dp = -(1/b^2) int (p f'' + f') J0(bp) dp  =>  |F| <= CT/b^2
def tail_const(k):
    f = CHI2[k][1:]
    assert pval(f, 1) == 0 and pval(pder(f), 1) == 0
    df = pder(f); d2 = pder(df)
    q = [Fr(0)]*max(len(df), len(d2) + 1)
    for i, x in enumerate(df):
        q[i] += x
    for i, x in enumerate(d2):
        q[i + 1] += x
    return abs_int_upper(q, 2000, Fr(0))


# ------------------------------------------------------------------ hbar majorant (Lemma 2) and its sup
class RB:
    """R_n(x) >= |r^(n)(x)| (r = Phi + a1/x^2), n = 0, 1, 2: min of the exact finite integration-by-parts sum (|.| termwise; each
    term decreasing in x) and the crude form int|psi| + c_n a1/x^(2+n) (c = 1, 2, 6).  Copied from gaplarge/gl7_routea_uniform.RB."""

    def __init__(self, w):
        n = len(w.c)
        self.nonosc = [(j, abs(w.d0[j])) for j in range(3, n, 2) if w.d0[j] != 0]
        self.osc = [(j, abs(w.d1[j])) for j in range(5, n) if w.d1[j] != 0]
        self.a1 = A(w.a1); self.I = A(w.absint_up)

    def asym(self, x, n):
        s = arb(0)
        for j, c in self.nonosc:
            f = 1
            for i in range(1, n + 1):
                f *= (j + i)
            s += A(c)*f/x**(j + 1 + n)
        for j, c in self.osc:
            a = j + 1
            fac = arb(1) if n == 0 else (1 + a/x if n == 1 else 1 + 2*a/x + a*(a + 1)/x**2)
            s += A(c)*fac/x**a
        return s

    def crude(self, x, n):
        return self.I + [1, 2, 6][n]*self.a1/x**(2 + n)

    def R(self, x, n):
        a, b = self.asym(x, n), self.crude(x, n)
        return a if a.upper() <= b.upper() else b


def hbar_cell(rb, xa, xb, asym_only=False):
    """sup over x in [xa, xb] of hbar(x) = (3x^5/640 + x^5/192)|r'(x)| + (3x^4/128)|r(x)| + (x^6/1152) sup_{[5x/6, x]}|r''|."""
    R = (lambda x, n: rb.asym(x, n)) if asym_only else rb.R
    q = A(Fr(5, 6))*xa
    return up((3*xb**5/640 + xb**5/192)*R(xa, 1) + 3*xb**4/128*R(xa, 0) + xb**6/1152*R(q, 2))


PSIMOM = {}


def r_derivs_point(k, x, nmax=3):
    """enclosures of r^(n)(x), n = 0..nmax, at a rational point x > 0: Phi(x) = sum_m (-1)^m mu_2m x^2m/(2m)! (exact moments of psi;
    tail |mu_j| <= int|psi|), r = Phi + a1/x^2."""
    x = Fr(x)
    old = ctx.prec; ctx.prec = 128 + int(3*float(x)) + 64
    try:
        M = int(2*float(x)) + 40
        if (k, M) not in PSIMOM:
            c = PSI[k]
            PSIMOM[(k, M)] = [sum(ci/(i + j + 1) for i, ci in enumerate(c)) for j in range(2*M + nmax + 2)]
        mu = PSIMOM[(k, M)]
        X = A(x); I = A(G.W[k].absint_up)
        out = []
        for n in range(nmax + 1):
            # Phi^(n)(x) = sum_{j even} (-1)^{j/2} mu_j d^n/dx^n x^j/j!
            s = arb(0)
            for j in range(n, 2*M + 1):
                if j % 2:
                    continue
                t = A(mu[j])*X**(j - n)/factorial(j - n)
                s += t if (j//2) % 2 == 0 else -t
            # tail j > 2M: |.| <= I x^(j-n)/(j-n)! summed, ratio x^2/((j-n+1)(j-n+2)) <= 1/4
            j0 = 2*M + 2
            tl = I*X**(j0 - n)/factorial(j0 - n)*4/3
            s += arb(0, up(tl).upper())
            ff = 1
            for i in range(n):
                ff *= -(2 + i)
            s += A(G.W[k].a1)*ff/X**(2 + n)
            out.append(s)
    finally:
        ctx.prec = old
    return out


def D_bound(k, n, xa):
    """sup_{x >= xa} |r^(n)(x)| <= int p^n |psi| + (n+1)! a1 / xa^(n+2) (crude; int p^n|psi| <= int|psi|)."""
    return A(G.W[k].absint_up) + factorial(n + 1)*A(G.W[k].a1)/A(xa)**(n + 2)


def R_taylor_cells(k, x0, x1, h):
    """cells [a, a+h] of [x0, x1]: sup |r^(n)| <= |r^(n)(c)| + (h/2)|r^(n+1)(c)| + (h/2)^2/2 sup|r^(n+2)|, n = 0, 1, 2."""
    cells = []
    a = Fr(x0)
    while a < x1:
        b = min(Fr(x1), a + h); c = (a + b)/2; hh = A((b - a)/2)
        d = r_derivs_point(k, c, 4)
        R = [up(abs(d[n]) + hh*abs(d[n + 1]) + hh*hh/2*D_bound(k, n + 2, a)) for n in range(3)]
        cells.append((a, b, R))
        a = b
    return cells


def H_sup(k, xT0=Fr(1, 2), xT1=Fr(48), hT=Fr(1, 100), xbig=Fr(4000)):
    """certified H_k >= sup_{x > 0} hbar_k(x).  (0, 1e-3]: crude closed form (every term increasing in x -> value at 1e-3);
    [1e-3, xT0]: crude/asymptotic majorants on geometric cells; [xT0, xT1]: Taylor cells of r (R_taylor_cells) for r, r' and r''
    (the r'' sup over [5x/6, x] is the max over the covering Taylor cells, or the R-majorant below xT0); [xT1, xbig]: R-majorants;
    x >= xbig: every term of the asymptotic majorant is non-increasing -> its value at xbig."""
    w = G.W[k]; rb = RB(w)
    a1, I = rb.a1, rb.I
    x0 = A(Fr(1, 1000)); c65 = A(Fr(6, 5))
    best = up((3*x0**5/640 + x0**5/192)*(I + 2*a1/x0**3) + 3*x0**4/128*(I + a1/x0**2) + x0**6/1152*(I + 6*a1*c65**4/x0**4))
    where = ('closed', 0)
    def upd(v, tag):
        nonlocal best, where
        if v.upper() > best.upper():
            best, where = v, tag
    edges = [Fr(1, 1000)]
    while edges[-1] < xT0:
        edges.append(min(xT0, edges[-1]*Fr(101, 100)))
    for a, b in zip(edges[:-1], edges[1:]):
        upd(hbar_cell(rb, A(a), A(b)), ('crude', float(a)))
    TC = R_taylor_cells(k, xT0, xT1, hT)
    for i, (a, b, R) in enumerate(TC):
        qa = Fr(5, 6)*a
        R2 = R[2]
        if qa < xT0:
            R2 = R2.max(rb.R(A(qa), 2)) if hasattr(R2, 'max') else R2
            if rb.R(A(qa), 2).upper() > R2.upper():
                R2 = up(rb.R(A(qa), 2))
        for (a2, b2, Rb) in TC:
            if b2 > qa and a2 < b and Rb[2].upper() > R2.upper():
                R2 = Rb[2]
        xb = A(b)
        upd(up((3*xb**5/640 + xb**5/192)*R[1] + 3*xb**4/128*R[0] + xb**6/1152*R2), ('taylor', float(a)))
    e = [xT1]
    while e[-1] < xbig:
        e.append(min(xbig, e[-1] + (Fr(1, 10) if e[-1] < 200 else 1)))
    for a, b in zip(e[:-1], e[1:]):
        if Fr(5, 6)*a < xT1:
            # r'' sup over [5a/6, b]: Taylor cells below xT1, majorant above
            R2 = up(rb.R(A(xT1), 2))
            for (a2, b2, Rb) in TC:
                if b2 > Fr(5, 6)*a and Rb[2].upper() > R2.upper():
                    R2 = Rb[2]
            xb = A(b)
            upd(up((3*xb**5/640 + xb**5/192)*rb.R(A(a), 1) + 3*xb**4/128*rb.R(A(a), 0) + xb**6/1152*R2), ('maj', float(a)))
        else:
            upd(hbar_cell(rb, A(a), A(b)), ('maj', float(a)))
    tail = hbar_cell(rb, A(xbig), A(xbig), asym_only=True)
    upd(tail, ('tail', float(xbig)))
    return best, tail, where


# ------------------------------------------------------------------ M1 lower bounds
CELLS = G.m1_cells()


def M1_small(beta):
    """[P] (GAP_LARGESPIN §3.2) lower bound of M1 on (0, beta] for beta <= 0.1."""
    nu1 = arb(2).sqrt()/beta
    return lo((arb(1)/24 - 7*beta**2/1920)/(beta/2).cosh() - 2*(-2*PI*nu1).exp()*(nu1/(2*PI) + 1/(4*PI**2)))


def M1min_upto(bmax):
    """certified lower bound of min_{0 < beta <= bmax} M1(beta)."""
    v = M1_small(A(Fr(1, 10)))
    for a, b, m in CELLS:
        if A(a) < bmax and m.lower() < v.lower():
            v = m
    return v


# ------------------------------------------------------------------ flat-Regge remainder constants (Part L)
def flat_E_consts(k):
    """psi = p f; Omega_flat(b) = -f1/b^3 + E_fl, |E_fl| <= e(b)/b^3,
    e(b) = |f1|(1 + 2/b) + |g(1)| + (|g'(1)| + int|p g'' + g'|)/b,  g = sum_{j>=2} j^2 f_j p^{j-2}  (T0_CERT Lemma 5)."""
    f = PSI[k][1:]
    assert pval(f, 1) == 0 and pval(pder(f), 1) == 0
    f1 = f[1]
    g = [j*j*f[j] for j in range(2, len(f))]
    dg = pder(g); d2g = pder(dg)
    q = [Fr(0)]*max(len(dg), len(d2g) + 1)
    for i, x in enumerate(dg):
        q[i] += x
    for i, x in enumerate(d2g):
        q[i + 1] += x
    return dict(f1=abs(f1), g1=abs(pval(g, 1)), dg1=abs(pval(dg, 1)), Iq=abs_int_upper(q, 2000, Fr(0)), f1_signed=f1)
