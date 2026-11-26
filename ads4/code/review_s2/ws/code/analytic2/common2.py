"""analytic2/common2.py: shared helpers for ANALYTIC_SLIVER2 (DV + certified)."""
import os, sys, math
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
AN = os.path.abspath(os.path.join(HERE, '..', 'analytic'))
sys.path.insert(0, AN)
from flint import arb, acb, ctx


def n_of_L(L, a, b=0):
    """integer n (spin - b) with -log d closest to L at twist a."""
    s = 2*b + 1
    d = math.exp(-L)
    # d = a(2n+a+s)/((n+a)(n+a+s)) ~ 2a/(n+a)
    lo, hi = 1, 10**12
    def dd(n):
        N = n + a
        return a*(2*n + a + s)/(N*(N + s))
    while hi - lo > 1:
        m = (lo + hi)//2
        if dd(m) > d: lo = m
        else: hi = m
    return lo


def S_rec(n, a, s, alpha, prec=256):
    """S_n(alpha) = L_{n,N}((1-z)^{-alpha}) by the dual-Hahn recurrence in n at fixed N = n + a (exact recurrence)."""
    old = ctx.prec; ctx.prec = prec
    N = arb(n) + arb(a)
    al = acb(alpha) if not isinstance(alpha, acb) else alpha
    if s == 1:
        Sm, S = acb(1), acb(1)
        for m in range(n):
            eps = 2*(2*m + 1)/((N - m)*(N + m + 1))
            Sm, S = S, Sm + eps*al*S
    else:
        # general s: alpha S_m = B_m S_{m+1} - (B_m + D_m) S_m + D_m S_{m-1}
        Sm = acb(0); S = acb(1)
        for m in range(n):
            Bm = (N + m + s)*(m + s)*(N - m)/((2*m + s)*(2*m + s + 1))
            Dm = -m*(N + m + s)*(N - m)/((2*m + s - 1)*(2*m + s)) if m > 0 else arb(0)
            Sm, S = S, (al*S + (Bm + Dm)*S - Dm*Sm)/Bm
    ctx.prec = old
    return S


def dstate(n, a, s):
    N = arb(n) + arb(a)
    return 1 - arb(n*(n + s))/(N*(N + s))
