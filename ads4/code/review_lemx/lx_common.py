"""lx_common.py (review of Lemma X, ANALYTIC_SLIVER §5.5): shared evaluators.

b = 0 (s = 1) dual-Hahn spin recurrence at fixed real N = n + a:
    S_{k+1} = S_{k-1} + eps_k alpha S_k,  S_{-1} = S_0 = 1,  eps_k = 2(2k+1)/((N-k)(N+k+1)),
    S^ = d^alpha S_n,  d = 1 - n(n+1)/(N(N+1)) = a(2N-a+1)/(N(N+1)).
* S_hat_arb: rigorous Arb ball (python-flint), precision chosen by the caller; returns (log|S^| lower, upper).
  NOTE: Arb's interval wrapping makes the radius grow like prod(1 + |eps_k alpha|) ~ exp(2|alpha| log(1/d)),
  i.e. ~ 10^4 bits at L ~ 3.7, tau ~ 1000.  Too low a precision returns a wide (but valid) ball.
* S_hat_vec: double-precision screening (vectorised over alphas, log-rescaled).
* bound_indep: an independent re-implementation (written from the proof, not copied) of the Lemma X bound, phase by phase.
"""
import math, sys
import numpy as np

sys.path.insert(0, '/home/user/MAJOR/ads4/code/analytic')
LAM2 = 998*999
LAMf = math.sqrt(LAM2)


def d_of(n, a):
    N = n + a
    return a*(2*N - a + 1)/(N*(N + 1))


def n_at(a, dt):
    """smallest n with d(n, a) <= dt (d decreasing in n)."""
    lo, hi = 1, 10**12
    while hi - lo > 1:
        m = (lo + hi)//2
        if d_of(m, a) > dt:
            lo = m
        else:
            hi = m
    return hi


def S_hat_vec(n, a, alphas):
    al = np.asarray(alphas, dtype=complex)
    N = n + a
    Sm = np.ones_like(al); S = np.ones_like(al); lg = np.zeros(al.shape)
    for m in range(n):
        eps = 2*(2*m + 1)/((N - m)*(N + m + 1))
        Sm, S = S, Sm + eps*al*S
        if m % 64 == 63:
            sc = np.maximum(np.abs(S), np.abs(Sm))
            S = S/sc; Sm = Sm/sc; lg += np.log(sc)
    d = d_of(n, a)
    return np.log(np.abs(S)) + lg + al.real*math.log(d)


def S_hat_arb(n, a, sig, tau, prec):
    from flint import arb, acb, ctx
    old = ctx.prec
    ctx.prec = prec
    try:
        aa = arb(a) if not isinstance(a, str) else arb(a)
        N = aa + n
        al = acb(arb(sig), arb(tau))
        Sm = acb(1); S = acb(1)
        for m in range(n):
            eps = arb(2*(2*m + 1))/((N - m)*(N + m + 1))
            Sm, S = S, Sm + eps*al*S
        d = 1 - arb(n*(n + 1))/(N*(N + 1))
        v = abs(S*(d.log()*al).exp())
        lo = v.lower(); hi = v.upper()
        rel = float((v.rad()/v.mid())) if v.mid() > 0 else float('inf')
        return float(lo.log()) if lo > 0 else -float('inf'), float(hi.log()), rel
    finally:
        ctx.prec = old


def bound_indep(sig, t_lo, t_hi, c, a0=998.0, parts=False):
    """Independent re-derivation of the Lemma X bound (log), uniform over a >= a0, sigma' in [0, sig], tau in [t_lo, t_hi].
    Phase I (eigen-coordinates, d^sigma cancellation, excess integral, adiabatic integral) up to y' = yA + dy, dy = 2 t_hi/(a0+1)^2
    (proved step bound Delta y_k <= 2 tau/(m+1)^2 <= 2 tau/(a+1)^2, from d eps/dk <= 4/m^2), last adiabatic term f(y')dy;
    window: eta = 2(delta + dy + u) (covers the first window step taken from an index with y in (yA - dy, yA]),
    W <= t_hi(1/yA - 1/yB) + 1 interior steps plus the entry step; phase III as in rt_pp.cross but with |r| >= sqrt(y^2 - 1 - u^2) everywhere and
    |lam+| <= (y + u) + sqrt(y^2 - 1 + u^2 + 2uy).  Returns log bound (and the parts)."""
    s = sig; tl = t_lo; th = t_hi; a = a0
    al = math.hypot(s, th); q = al/tl
    dl = c*tl**(-2/3); yA = 1 - dl; yB = 1 + dl
    u = s/(a + 1)
    dy = 2*th/(a + 1)**2
    yAp = yA + dy
    assert yAp < 1
    L1 = 1 - yAp*yAp
    c0 = 1 + 2*al/a**2
    fA = L1**-0.5 + q*yAp/L1
    # final conversion if the walk ends in phase I: |lam+| + |lam-| <= 2 cosh(u/sqrt(1 - y^2)) at y <= yA
    endI = math.log(2*math.cosh(u/math.sqrt(L1)))
    cI = math.log(c0) + s*math.log(2/(1 + math.sqrt(L1))) + q*(math.asin(yAp) + q/2*math.log(1/L1)) + q*fA*dy
    lI_end = cI + endI
    # window
    eta = 2*(dl + u)
    etap = 2*(dl + dy + u)          # first window step
    e1 = 2*(u + 1 - yA + dy + math.sqrt(1 - (yA - dy)**2 + u*u + 2*u))
    W = th*(1/yA - 1/yB) + 1          # steps with y in (yA, yB); the entry step from y_K <= yA is step1
    step = math.log(1 + math.sqrt(eta) + eta)
    step1 = math.log(max(1 + etap + etap/math.sqrt(eta), 1 + math.sqrt(eta)))
    cII = cI + math.log(2 + e1 + e1/math.sqrt(eta)) + step1 + W*step     # bound of the window norm N at exit / end
    lII_end = cII
    # phase III
    ymax = th/(a + 1)
    ym = max(ymax, yB*1.000001)
    rB = math.sqrt(yB*yB - 1 - u*u)
    e2 = u + yB + dy - 1 + math.sqrt((yB + dy)**2 - 1 + u*u + 2*u*(yB + dy))
    conv = math.log((math.sqrt(eta)*(1 + e2) + e2)/(2*rB))
    G = th*math.acos(min(1.0, a/th)) - a*math.acosh(max(1.0, th/a)) if th > a else 0.0
    rr = lambda y: math.sqrt(y*y - 1 - u*u)
    # adiabatic: q int_{yB}^{ym} [1/r + q y/r^2] dy  (r = sqrt(y^2 - 1 - u^2)), plus one boundary term
    yb2 = yB*yB - 1 - u*u; ym2 = ym*ym - 1 - u*u
    lad = q*(math.log((ym + math.sqrt(ym2))/(yB + math.sqrt(yb2))) + q/2*math.log(ym2/yb2)) + q*(1/rB + q*(yB + dy)/rB**2)*dy
    lsig = s*(math.log(max(th/a, 1.0)) + 1/a)/rB
    lmax = (ym + u) + math.sqrt(ym*ym - 1 + u*u + 2*u*ym)
    fin = math.log(lmax*lmax + 1)     # (|lam+| + |lam-|)|lam+| = |lam+|^2 + 1, |lam+| <= lmax
    lIII_end = cII + conv + G + lad + lsig + fin
    res = max(lI_end, lII_end, lIII_end)
    if parts:
        return res, dict(I=lI_end, II=lII_end, III=lIII_end)
    return res


def cross_float(sig, t_lo, t_hi, c, a=998.0):
    import rt_pp
    return math.log(rt_pp.cross(sig, t_lo, t_hi, c, a))


def PP0_log(sig, t_lo, t_hi):
    import rt_pp
    return math.log(rt_pp.PP0(sig, t_lo, t_hi))
