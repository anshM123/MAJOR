"""
fullpos3/fp3_ht_large.py (FULLPOS5.md): the FULLPOS3 §4 high-twist certificate (fp3_ht.py) at large impact parameter, beta > 8.4.

fp3_ht.py is imported unchanged; this module replaces, IN THIS PROCESS ONLY, the pieces that break or lose at large beta
(c = tanh^2(beta/2), d = 1 - c ~ 4 e^{-beta}).  Every replacement is a valid upper bound by itself (FULLPOS5 §1):
  * precision: ctx.prec = 256 + 2*bb/ln2 bits (1 - z ~ d loses log2(1/d) bits);
  * circle_sups: 1 - z = (1 - c_l) - c_l expm1(t) (no cancellation); the arc-subdivision floor is RELATIVE to the piece length
    (the absolute 1e-12 floor of fp3_ht is larger than the whole stadium once d < 1e-11: the reported 'stadium geometry failed');
  * jval: fine beta-grid near the cell (step 0.02) + the scaling lemma J_X(beta) <= e^{-(beta-g)/2} J_X(g) for g <= beta
    (sinh(x+delta) >= e^delta sinh x and A_X decreasing); fp3_ht's grid stops at 9.9 and would freeze J at its value there;
  * residual (EPSW-cell) components and the m = 0 residual G^r: min(Lemma B, Lemma B') per frequency cell, where
    Lemma B': |P_{-1/2+i nu}(cosh xi)| <= (sqrt2/pi) coth(pi nu) gamma_beta e^{-beta/2} e^{nu eta} Q(nu, delta; beta)   (FULLPOS5 §1.4)
    (Mehler-sin + one integration by parts; beta = Re xi, eta = |Im xi| < pi/2; any delta > 0);
  * norm_major_float: written in d (float c rounds to 1 for beta > 37);  [choice of eta only; the bound itself is the Arb one]
  * far_remainder: rho-grid relative to h_max (fp3_ht uses the absolute grid 1e-20..1, which has no admissible rho when
    h_max < 1e-20, i.e. beta > ~20), convergence cap relative to 1/d.
  * F_cert: LRINF cells (beta <= 64) as before.
Usage:
  fp3_ht_large.py b ba bb a0 [M]                 single run
  fp3_ht_large.py scan "ba,bb;ba,bb;..." out.jsonl   (out.jsonl in cache/, same line format as fp3_ht scan)
"""
import sys, os, math, json, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, acb, ctx
import numpy as np
import fp3_ht as H
import fp3_htdata as D
from fp3_htdata import A, LAM, PI, up, lo


def set_prec(bb):
    ctx.prec = max(256, 256 + int(2*float(bb)/math.log(2)))


# ---------------------------------------------------------------- J tables: fine grid + scaling lemma
_JG = {}
JSTEP = Fr(1, 50)

def _jgrid_point(k, g):
    key = (k, g)
    if key not in _JG:
        _JG[key] = D.jtable(H.PSI2 if k == 2 else H.PSI3, [g], ng=400)[0]
    return _JG[key]

def jval_large(k, beta_lo):
    """upper bounds of J_s, J_C, J_e at beta >= beta_lo: grid g = JSTEP*floor(beta_lo/JSTEP) <= beta_lo and
    J_X(beta) <= e^{-(beta-g)/2} J_X(g)."""
    bl = float(beta_lo.lower())
    if bl < 4.0:
        return H._jval_orig(k, beta_lo)
    g = JSTEP*int(math.floor(bl/float(JSTEP)) - 1)
    while A(g) > beta_lo:           # make sure g <= beta_lo rigorously
        g -= JSTEP
    fac = up((-(lo(beta_lo) - A(g))/2).exp())
    Js, Jc, Je = _jgrid_point(k, g)
    return (up(Js*fac), up(Jc*fac), up(Je*fac))


# ---------------------------------------------------------------- Lemma B' (FULLPOS5 §1.4)
_SQ2 = arb(2).sqrt()

def lemmaBp_factor(nu_lo, beta):
    """returns K with
    |P_{-1/2+i nu}(cosh xi)| <= K e^{-beta_xi/2} e^{nu eta} for every nu >= nu_lo > 0 and Re xi >= beta > 0, |Im xi| = eta < pi/2.
    K = (sqrt2/pi) coth(pi nu_lo) gamma Q,  gamma = (1 - e^{-2 beta})^{-1/2},
    Q = min_delta { 2 sqrt2 asin(sqrt(1 - e^{-delta})) + [2 s(delta) + (sqrt2/2)(coth beta - 1) asin(e^{-delta/2})]/nu_lo },
    s(delta) = sqrt2 e^{-delta/2}/sqrt(1 - e^{-delta})."""
    nu = arb(nu_lo); beta = arb(beta)
    gam = 1/(1 - (-2*beta).exp()).sqrt()
    cthb = 1/(beta.tanh())
    best = None
    for i in range(-40, 41):
        de = arb(2)**(arb(i)/4)/nu if i <= 36 else arb(2)**(arb(i - 36))
        x = (-de).exp()
        sd = _SQ2*(-de/2).exp()/(1 - x).sqrt()
        Q = 2*_SQ2*((1 - x).sqrt()).asin() + (2*sd + (_SQ2/2)*(cthb - 1)*((-de/2).exp()).asin())/nu
        if best is None or Q < best:
            best = Q
    cth = 1/((PI*nu).tanh())
    return up((_SQ2/PI)*cth*gam*best)


_EPSLO = None
def eps_nu_lo():
    """lower frequency edge of each residual component returned by rcomp_split (EPSW cells, then the three pieces of fp3_htdata's
    'rest': Gaussian tail (nu >= 1059), theta part 1 (all nu: Lemma B only), theta part 2 (nu > Lam))."""
    global _EPSLO
    if _EPSLO is None:
        cells, _ = D.eps_cells()
        _EPSLO = [a for (a, b, w) in cells] + [1059.0, -1.0, float(LAM.lower())]
    return _EPSLO


def rint_components_split(eta):
    """fp3_htdata.rint_components with its last entry 'rest' = 2.1 tail + thc split into its three pieces (same formulas, so the
    three sum to 'rest'):  2.1 tail (nu >= 1059), th1 = 1e10 e^{-49} e^{eta (Lam+30)} (nu <= Lam + 30), th2 = 2.1*3.2e-10 e^{1059 eta}
    (Lam < nu <= 1059: theta correction where psi = 0, r = f)."""
    eta = arb(eta)
    cells, epsW = D.eps_cells()
    out = [up(arb(w)*(eta*arb(b)).exp()) for (a, b, w) in cells]
    S, cmax, smin, smax = D.gauss_tail_params()
    pref = 16*2*A(S)*(1/(2*A(smin)**2)).exp()
    s = A(smax); x0 = arb(1059) - A(cmax) - eta*s*s
    assert x0 > 20
    tail = pref*(2*PI).sqrt()*s*(eta*A(cmax) + eta*eta*s*s/2).exp()*(-(x0*x0)/(2*s*s)).exp()
    th1 = arb(1e10)*arb(-49).exp()*(eta*(LAM + 30)).exp()
    th2 = arb(2.1)*arb(3.2e-10)*(eta*arb(1059)).exp()
    out += [up(2.1*tail), up(th1), up(th2)]
    return out

_RCS = {}
def rcomp_split(eta_hi):
    e = math.ceil(float(eta_hi.upper())*1e4 + 1e-9)/1e4
    if e not in _RCS:
        _RCS[e] = rint_components_split(arb(e))
    return _RCS[e]

NU_MIN = 0.05          # Lemma B' is used on EPSW cells with lower edge >= NU_MIN (it is valid for any nu_lo > 0)
_KBP = {}
def kbp(i, beta_lo_grid):
    key = (i, beta_lo_grid)
    if key not in _KBP:
        nl = eps_nu_lo()[i]
        _KBP[key] = lemmaBp_factor(lo(arb(nl)*(1 - arb(1e-9))), A(beta_lo_grid)) if nl >= NU_MIN else None
    return _KBP[key]


def residual_kernel(i, beta, axi):
    """upper bound of |P_{-1/2+i nu}(cosh xi)| e^{-nu_hi eta} for nu in EPSW cell i (or 'rest' if i is None):
    min(Lemma B: (|xi|/beta) P_{-1/2}(cosh beta), Lemma B')."""
    lb = (axi/beta)*D.P_mhalf(beta)
    if i is None or eps_nu_lo()[i] < NU_MIN:
        return up(lb)
    bg = Fr(math.floor(float(beta.lower())*50) - 1, 50)  # grid point <= beta (factor monotone decreasing in beta)
    if bg > 7:
        bg = Fr(7)
    K = kbp(i, bg)
    v = K*(-beta/2).exp()
    return up(v if v < lb else lb)


# ---------------------------------------------------------------- Lemma A' (FULLPOS5 §1.5): the coth correction, integrated by parts
def ccor1_value(psi, eta):
    """C1(eta) with |int_0^1 psi(p) (coth(pi Lam p) - 1) sin(Lam p t) dp| <= C1(eta)/(Lam |t|) for |Im t| <= eta < pi/2:
    g(p) = 2 psi(p)/(e^{2 pi Lam p} - 1), g(1) = 0, g(0) = c_1/(pi Lam); one integration by parts and |cos(Lam p t)| <= cosh(Lam p eta):
      C1 = |c_1|/(pi Lam) + sum_n 2|c_n| (2 pi Lam)^{-n} I_n(kappa),   kappa = eta/(2 pi) < 1/4,
      I_n >= int_0^oo |d/dx (x^n/(e^x - 1))| cosh(kappa x) dx:
      I_1 <= 1 + (kappa/2)[(1-kappa)^{-2} + (1-kappa)^{-1}]            (x/(e^x-1) decreasing; sum_m (m-kappa)^{-2})
      I_n <= 2 F_n cosh(kappa n) + (kappa/2) n! [(1-kappa)^{-n-1} + (1-kappa)^{-n}/n],  F_n = (2(n-1)/e)^{n-1}  (n >= 2)."""
    eta = arb(eta)
    assert eta < PI/2
    ka = eta/(2*PI)
    tot = psi.absc[1]/(PI*LAM)
    for n, a in enumerate(psi.absc):
        if n == 0 or a == 0:
            continue
        if n == 1:
            In = 1 + (ka/2)*((1 - ka)**(-2) + 1/(1 - ka))
        else:
            Fn = (arb(2*(n - 1))/arb(1).exp())**(n - 1)
            In = 2*Fn*(ka*n).cosh() + (ka/2)*arb.fac_ui(n)*((1 - ka)**(-n - 1) + (1 - ka)**(-n)/n)
        tot += 2*a*(2*PI*LAM)**(-n)*In
    return up(tot)

_CC1 = {}
def ccor1(k, eta_hi):
    e = math.ceil(float(eta_hi.upper())*1e3 + 1e-9)/1e3
    if (k, e) not in _CC1:
        _CC1[(k, e)] = ccor1_value(H.PSI2 if k == 2 else H.PSI3, arb(e))
    return _CC1[(k, e)]


def coth_part(k, eta, beta, Jc):
    """bound on the coth-correction contribution: min(C_cor(eta) J_C, C1(eta) J_C/(Lam beta))  (|t| >= Re t >= beta)."""
    a = H.ccor(k, eta)*Jc
    b = ccor1(k, eta)*Jc/(LAM*beta)
    return a if a < b else b


# ---------------------------------------------------------------- geometry
def _piece_geom_large(cl, oml, Dl, R, kind, pa, pb):
    pb_ = arb((pa + pb)/2, (pb - pa)/2)
    if kind == 1:
        t = acb(pb_, R)
    elif kind == 0:
        t = acb(R*pb_.cos(), R*pb_.sin())
    else:
        t = acb(Dl + R*pb_.cos(), R*pb_.sin())
    z = acb(cl)*t.exp()
    omz = acb(oml) - acb(cl)*t.expm1()          # 1 - z without cancellation
    X = (2 - omz)/omz
    xi = X.acosh()
    return z, omz, xi


def circle_sups_large(cl, ch, R, b, tol=2e-4):
    """as fp3_ht.circle_sups (same components and bounds), large-beta safe; residual cells use min(Lemma B, Lemma B')."""
    out = {}
    def mx(key, v):
        if key not in out or v > out[key]:
            out[key] = v
    R = arb(R); cl = arb(cl); oml = 1 - cl; Dl = up((arb(ch)/cl).log())
    stack = [(0, PI/2, PI), (1, arb(0), Dl), (2, arb(0), PI/2)]
    plen = {0: float((PI/2).upper()), 1: float(Dl.upper()), 2: float((PI/2).upper())}
    nlo = eps_nu_lo()
    while stack:
        kind, pa, pb = stack.pop()
        z, omz, xi = _piece_geom_large(cl, oml, Dl, R, kind, pa, pb)
        ok = xi.real.is_finite() and xi.imag.is_finite()
        if ok:
            beta = lo(xi.real); eta = up(abs(xi.imag))
            ok = beta > 0 and eta < PI/2 and float(xi.imag.rad()) < tol and float(xi.real.rad()) < 0.01*float(beta.lower())
        if not ok:
            if float((pb - pa).upper()) < 1e-12*plen[kind]:
                raise RuntimeError('stadium geometry failed at R=%s' % R)
            m = (pa + pb)/2
            stack.append((kind, pa, m)); stack.append((kind, m, pb))
            continue
        axi = up(abs(xi)); om = lo(abs(omz))
        Js, Jc, Je = jval_large(2, beta)
        k2 = 2*om**(-3)*(2**0.5/PI)
        mx('As', up(k2*(Js + coth_part(2, eta, beta, Jc))))
        mx('Ae', up(k2*(LAM*eta).exp()*Je))
        kr0 = (2/LAM)*om**(-3)
        rc = rcomp_split(eta)
        lb = up((axi/beta)*D.P_mhalf(beta))                 # Lemma B kernel
        eb = up((-beta/2).exp())
        bg = Fr(math.floor(float(beta.lower())*50) - 1, 50)  # grid point <= beta (Lemma B' factor decreases in beta)
        if bg > 7:
            bg = Fr(7)                                       # K(beta) <= K(7) for beta >= 7 (relative change < e^{-13})
        for i, v in enumerate(rc):
            kern = lb
            if i < len(nlo) and nlo[i] >= NU_MIN:
                kb = kbp(i, bg)*eb
                if kb < kern:
                    kern = kb
            mx('r%d' % i, up(kr0*kern*v))
        Js3, Jc3, Je3 = jval_large(3, beta)
        k3 = LAM**2*om**(-4)*(2**0.5/PI)
        mx('Ds', up(k3*(Js3 + coth_part(3, eta, beta, Jc3))))
        mx('De', up(k3*(LAM*eta).exp()*Je3))
        mx('zinv', up(1/abs(z)))
    return out


# ---------------------------------------------------------------- majorant (float, choice of eta only) written in d
_DCELL = [0.0]

def _logA_vec_d(m, rho, s, dl, logeta):
    l = np.arange(1, m + 1, dtype=float)
    yb = (2*l + 1)/rho
    Bh = 2 + yb/2 - (yb/2)/np.tan(yb/2)
    Lc = -math.log1p(-dl)                                   # -log c_lo
    x = l*Lc
    lcoef = x + np.log(-np.expm1(-x)) - np.log(l)           # log((c^{-l} - 1)/l), stable for small and large l*Lc
    base = lcoef + np.log(Bh) + (l + 1)*np.log(rho) + l*np.log(rho + s) - np.log(l + 1)
    lt = base[None, :] + np.outer(logeta, l)
    mxv = lt.max(axis=1)
    return mxv + np.log(np.exp(lt - mxv[:, None]).sum(axis=1))

def norm_major_float_large(m, R, cl, s):
    dl = _DCELL[0]
    rho = 3*m/(math.e*R)
    if (2*m + 1)/rho >= 2*math.pi*0.99:
        return math.inf, None
    hi = (1 - dl)/(rho*(rho + s))
    logeta = math.log(hi) + np.linspace(-12, 2, 141)
    LA = _logA_vec_d(m, rho, s, dl, logeta)
    v = -m*logeta + np.exp(np.minimum(LA, 700))
    i = int(np.argmin(v))
    lf = math.lgamma(3*m + 1) - 3*m*math.log(R*rho)
    return v[i] + max(0.0, lf), math.exp(logeta[i])


# ---------------------------------------------------------------- far remainder: rho grid relative to h_max
def far_remainder_large(cb, hmax, nmin, M, b, tau, d_lo):
    """identical bound to fp3_ht.far_remainder (FULLPOS3 §4.4); only the candidate rho grid (now h_max * 10^{i/4}, i = -24..56,
    i.e. the same 81-point span around h_max that the absolute grid gives at beta ~ 8) and the iteration cap (now 10^4/d) differ.
    The candidate set is arbitrary in the lemma (any admissible rho), so the bound stays valid."""
    s = 2*b + 1
    kap = A(abs(Fr(9, 2) - b*(b - 1)))
    KcA = (PI/2)*(2*H.PSI2.absexp(PI/2, 2000) + (2/LAM)*D.Rint(PI/2)[0]) if 'KcA' not in _FRC else _FRC['KcA']
    KcD = LAM**2*(PI/2)*H.PSI3.absexp(PI/2, 2000) if 'KcD' not in _FRC else _FRC['KcD']
    _FRC['KcA'] = KcA; _FRC['KcD'] = KcD
    e3 = arb(3).exp(); e4 = arb(4).exp()
    def BA(k):
        return KcA*e3*arb(k + 3)**3/27
    def Bk(k):
        return BA(k) + tau*(arb(k + 1)**2/2*BA(k + 1) + (arb(k)**2/2 + 4*k + kap)*BA(k)) + tau*KcD*e4*arb(k + 4)**4/256
    cm = arb(cb.lower()); cp = arb(cb.upper()); h = arb(hmax)
    def logfac(x):
        v = (cp + x)/(1 + x)
        v2 = abs(cm - x)/(1 - x); v3 = abs(cp - x)/(1 - x)
        m = v
        if v2 > m: m = v2
        if v3 > m: m = v3
        return up(m.log()) if m > 0 else arb(-10**6)
    hm = arb(hmax.upper()) if isinstance(hmax, arb) else arb(hmax)
    rhos = [up(hm*arb(10)**(arb(i)/4)) for i in range(-24, 57)]
    LM = [arb(0) for _ in rhos]; valid = [True]*len(rhos)
    tot = arb(0); k0 = 0
    th = arb(3)/4
    th2 = (1 + cm)/2*arb(0.9)
    if th2 < th:
        th = th2
    KT = None
    kcap = int(1e4/float(d_lo.lower())) + 10**9
    while True:
        blk = max(1, k0//400)
        k1 = k0 + blk
        ua = arb(k0*(k0 + s)); ub = arb((k1 - 1)*(k1 - 1 + s))
        best = None
        for i, r in enumerate(rhos):
            if not valid[i]:
                continue
            if not (r*ub < 1):
                valid[i] = False; continue
            lf = logfac(r*ua); lf2 = logfac(r*ub)
            if lf2 > lf: lf = lf2
            LMblk = LM[i] + (blk*lf if lf > 0 else arb(0))
            Mk = LMblk.exp()
            q = h/r
            fb = cp**k0 + Mk*((M + 1) if q <= 1 else (M + 1)*q**M)
            cand = fb
            if q < 1:
                tr = Mk*q**(M + 1)/(1 - q)
                if tr < cand:
                    cand = tr
            if best is None or cand < best:
                best = arb(cand.upper())
            LM[i] = LM[i] + blk*lf
        if best is None:
            raise RuntimeError('no admissible rho at k=%d' % k1)
        tot += up(arb(blk)*Bk(k1)*best)
        k0 = k1
        if k0 > 1000:
            K = arb(k0)
            q_th = (cp + th)/(1 + th)
            for v in (cp, abs(cm - th)/(1 - th), abs(cp - th)/(1 - th)):
                if v > q_th: q_th = v
            hu = h*K*(K + s)/th
            rr = q_th*(1 + 1/K)**(2*M + 8)
            rr2 = cp*(1 + 1/K)**8
            if rr < 1 and rr2 < 1:
                fK = Bk(k0)*(cp**k0 + (M + 1)*q_th**k0*(hu**M if hu > 1 else arb(1)))
                if fK < arb(10)**(-400):
                    rmax = rr if rr > rr2 else rr2
                    tot += up(2*fK/(1 - rmax))
                    KT = k0
                    break
        if k0 > kcap:
            raise RuntimeError('far remainder did not converge')
    return up(tot), KT

_FRC = {}


# ---------------------------------------------------------------- m = 0 residual term with Lemma B'
def Gr_large(cb, beta_lo):
    """upper bound of |G^r(c)| on the cell: (2/Lam)(1 - c_hi)^{-3} sum_i w_i min(Lemma B, Lemma B')(beta_lo) at eta = 0
    (fp3_ht uses (2/Lam)(1 - c_hi)^{-3} P_{-1/2}(cosh beta_lo) Rint(0) = the same sum with Lemma B only)."""
    rc = rcomp_split(arb(0))
    nlo = eps_nu_lo()
    assert len(rc) == len(nlo)
    tot = arb(0)
    b = arb(beta_lo)
    for i, v in enumerate(rc):
        tot += v*residual_kernel(i, b, b)
    return up((2/LAM)*(1 - cb.upper())**(-3)*tot)


# ---------------------------------------------------------------- LRINF floor beyond the last cell (beta >= 64)
def fcert_tail(beta):
    """float lower bound of LRINF's floor at beta >= 64 (LRINF_TIGHT §1.5, [P + certified] there):
    F >= (1 - r_tail) A1min M1lo/Lam^2 (+ the flat term A2'/(2 bhat^3) > 0, dropped), r_tail <= 2.03e-7 (taken as 2.1e-7),
    A1min = 2 psi2'(0) - psi3'(0) (exact), M1lo = (sqrt2/pi) K0 ell/beta, K0 = (cosh beta - 1)^{-1/2},
    ell = 64 L(64), L(b) = 1/b - 1/(2 sinh(b/2)) - sqrt2 (e^{-b/2}/4.7 + 4/b^2) - 2 sqrt(2 ln2)/(b - ln2)^2 (b L(b) increasing).
    Decreasing in beta, so its value at the cell's right end bounds the cell."""
    b = A(Fr(beta)); assert b >= 64
    b64 = arb(64)
    L64 = 1/b64 - 1/(2*(b64/2).sinh()) - arb(2).sqrt()*((-b64/2).exp()/arb(4.7) + 4/b64**2) \
        - 2*(2*arb(2).log()).sqrt()/(b64 - arb(2).log())**2
    ell = 64*L64
    A1 = 2*A(H.PSI2.c[1]) - A(H.PSI3.c[1])
    K0 = 1/(b.cosh() - 1).sqrt()
    F = (1 - arb(2.1e-7))*A1*(arb(2).sqrt()/PI)*K0*ell/(b*LAM**2)
    return float(F.lower())*(1 - 1e-12)


# ---------------------------------------------------------------- install
def install():
    if getattr(H, '_large_installed', False):
        return
    H._jval_orig = H.jval
    H.jval = jval_large
    H.circle_sups = circle_sups_large
    H.norm_major_float = norm_major_float_large
    H.far_remainder = far_remainder_large
    H._large_installed = True


def run_large(b, ba, bb, a0, M=None, m1=70, verbose=True, geo=None):
    """fp3_ht.run, line for line, except: precision set from bb; _DCELL for the float eta choice; G^r by Gr_large
    (min of Lemma B and Lemma B' per EPSW cell); far_remainder_large; norm_major_float_large."""
    t0 = time.time()
    set_prec(Fr(bb))
    _DCELL[0] = float(lo(1 - up((A(Fr(bb))/2).tanh()**2)).lower())
    s = 2*b + 1; sig = 1 if b == 0 else -1; kap = Fr(9, 2) - b*(b - 1)
    ba = Fr(ba); bb = Fr(bb)
    c_lo = (A(ba)/2).tanh()**2; c_hi = (A(bb)/2).tanh()**2
    cb = lo(c_lo).union(up(c_hi))
    d_lo = 1 - up(c_hi); w3 = (1 - lo(c_lo))**3
    # floor
    fc = [F for (x, y, F) in D.fcert_cells() if not (y <= ba or x >= bb)]
    if bb > 64:
        fc.append(fcert_tail(bb))
    Fcert = min(fc)
    nmin, hmax = H.state_params(c_lo, a0, s)
    tau = hmax                                            # tau = 1/T <= h
    if geo is None:
        geo = H.geometry(ba, bb, verbose)
    sups, Rmax = geo
    keysA = [k for k in sups[0][1] if k == 'As' or k == 'Ae' or k[0] == 'r']
    keysD = ['Ds', 'De']
    # per component: S2(R) = S(R) + tau Dfac(R) S(R2) for G_A components; tau S(R) for D components
    def at_least(Rq):
        for (R, S) in sups:
            if R >= Rq:
                return R, S
        return sups[-1]
    comp = {k: [] for k in keysA + keysD}
    for i, (R, S) in enumerate(sups):
        for k in keysA:
            best = None
            for (R2, S2_) in sups[i + 1:]:
                if not S2_[k].is_finite():
                    break
                dl = arb(R2 - R); zinv = S2_['zinv']
                Dfac = (1 + zinv)/2*(2/dl**2 + 1/dl) + (abs(A(Fr(7, 2) - 2*b)) + zinv/2)/dl + abs(A(kap))
                v = S[k] + tau*Dfac*S2_[k]
                if best is None or v < best:
                    best = v
            comp[k].append((R, up(best) if (best is not None and best.is_finite()) else arb('inf')))
        for k in keysD:
            comp[k].append((R, up(tau*S[k]) if S[k].is_finite() else arb('inf')))
    # m = 0 extras: G^r(c) and tau D_b G_A(c)
    beta_lo = A(ba)
    Gr = Gr_large(cb, beta_lo)
    best_tD = arb(0)
    for k in keysA:
        bk = None
        for (R, S) in sups:
            Rr = arb(R)
            v = tau*((1 - lo(c_lo))/(2*lo(c_lo))*(2/Rr**2 + 1/Rr) + (abs(A(Fr(7, 2) - 2*b)) + 1/(2*lo(c_lo)))/Rr + abs(A(kap)))*S[k]
            if bk is None or v < bk:
                bk = up(v)
        best_tD += bk
    # local orders, exact P_m for m <= m1: sum over components of min over R
    Pex = H.Pm_exact(cb, s, m1)
    loc = []; loc_comp = {}
    for m in range(1, m1 + 1):
        tot = arb(0)
        nr = {R: H.normR(Pex[m], R) for (R, _S) in sups}
        for k, lst in comp.items():
            bm = None
            for (R, St) in lst:
                if not St.is_finite():
                    continue
                v = hmax**m*nr[R]*St
                if bm is None or v < bm:
                    bm = up(v)
            tot += bm
            loc_comp[k] = loc_comp.get(k, 0.0) + float((w3*bm).upper())
        loc.append(up(tot))
        if verbose:
            print(f'      exact m={m}: {float(tot.upper()):.3e}', flush=True)
    # choose M: needs the far remainder; try list
    if M:
        Ms = [M]
    else:
        Ms = []
        for MM in [100, 200, 300, 450, 650, 900, 1200, 1600]:
            try:
                fr, _k = far_remainder_large(cb, hmax, int(float(nmin.lower())), MM, b, tau, d_lo)
            except (RuntimeError, AssertionError):
                continue
            if w3*fr < arb(Fcert)*1e-3:
                Ms = [MM]; break
        if not Ms:
            out = dict(b=b, ba=str(ba), bb=str(bb), a0=a0, Fcert=Fcert, ok=False, err=float('inf'), why='far remainder',
                       nmin=float(nmin.lower()), hmax=float(hmax.upper()), secs=time.time() - t0)
            return out
    best_total = None
    for MM in Ms:
        tail_loc = arb(0); mm_detail = []
        cl = float(cb.lower())
        for m in range(m1 + 1, MM + 1):
            # float log-majorant per radius (shared by components)
            nm = {}
            for (R, _S) in sups:
                if R > 3:
                    continue
                lv, eta = norm_major_float_large(m, R, cl, s)
                if eta is not None:
                    nm[R] = (lv, eta)
            bm_tot = arb(0); percomp = {}
            for k, lst in comp.items():
                cand = [(nm[R][0] + float(St.log().upper()), R) for (R, St) in lst if R in nm and St.is_finite() and St > 0]
                if not cand:
                    continue
                lvb, R = min(cand)
                St = dict(lst)[R]
                vv = hmax**m*H.norm_major(m, R, cb, s, nm[R][1])*St
                bm_tot += vv; percomp[k] = (float(vv.upper()), R)
            bm = up(bm_tot)
            if verbose and (m <= m1 + 2 or m % 100 == 0):
                worst = max(percomp.items(), key=lambda kv: kv[1][0])
                print(f'      m={m}: {float(bm.upper()):.3e}  worst component {worst[0]} {worst[1][0]:.2e} at R={worst[1][1]:.2e}', flush=True)
            tail_loc += bm
            if m <= m1 + 3 or m % 50 == 0:
                mm_detail.append((m, float(bm.upper())))
        try:
            far, kend = far_remainder_large(cb, hmax, int(float(nmin.lower())), MM, b, tau, d_lo)
        except RuntimeError as e:
            far, kend = arb('inf'), None
        err = w3*(Gr + best_tD + sum(loc, arb(0)) + tail_loc + far)
        res = dict(M=MM, err=float(err.upper()), far=float((w3*far).upper()), tail_loc=float((w3*tail_loc).upper()), kend=kend,
                   mdet=mm_detail)
        if best_total is None or err < best_total[0]:
            best_total = (err, res)
        if verbose:
            print(f'   M={MM}: err<={float(err.upper()):.3e}  far<={res["far"]:.3e}  higher local<={res["tail_loc"]:.3e}  kend={kend}', flush=True)
        if err < arb(Fcert)/2:
            break
    err, res = best_total
    ok = err < arb(Fcert)
    out = dict(b=b, ba=str(ba), bb=str(bb), a0=a0, Fcert=Fcert, ok=bool(ok), loc_comp=loc_comp,
               Gr=float((w3*Gr).upper()), tD=float((w3*best_tD).upper()), loc=[float((w3*x).upper()) for x in loc],
               nmin=float(nmin.lower()), hmax=float(hmax.upper()), secs=time.time() - t0, **res)
    return out


def geometry_large(ba, bb, verbose=False):
    set_prec(Fr(bb))
    return H.geometry(ba, bb, verbose)


install()
H.run = run_large
# finer twist ladder for the bisection (fp3_ht.scan_cell reads the module-level A0LIST; must contain 12500 and 2200000)
H.A0LIST = [998, 1100, 1250, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2125, 2250, 2375, 2500, 2650, 2800, 2900, 3000, 3100, 3200,
            3400, 3600, 3800, 4000, 4250, 4500, 4750, 5000, 5500, 6000, 7000, 8500, 10000, 12500, 15000, 20000, 25000, 30000, 40000,
            50000, 70000, 100000, 200000, 500000, 2200000]


def scan_cell_large(ba, bb, out):
    set_prec(Fr(bb))
    _DCELL[0] = float(lo(1 - up((A(Fr(bb))/2).tanh()**2)).lower())
    return H.scan_cell(ba, bb, out)


if __name__ == '__main__':
    if sys.argv[1] == 'scan':
        cells = [l.split(',') for l in sys.argv[2].split(';')]
        out = os.path.join(HERE, 'cache', sys.argv[3])
        done = set()
        if os.path.exists(out):
            for l in open(out):
                r = json.loads(l); done.add((str(Fr(r['ba'])), str(Fr(r['bb']))))
        for ba, bb in cells:
            if (str(Fr(ba)), str(Fr(bb))) in done:
                continue
            scan_cell_large(ba, bb, out)
    else:
        b = int(sys.argv[1]); ba = sys.argv[2]; bb = sys.argv[3]; a0 = int(sys.argv[4])
        M = int(sys.argv[5]) if len(sys.argv) > 5 else None
        r = run_large(b, ba, bb, a0, M)
        print(json.dumps(r))
