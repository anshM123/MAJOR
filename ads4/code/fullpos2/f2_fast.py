"""[DV] fast exact bare kernel K_p(nu) = sum_k U_k g^{(p)}_k(nu) in python-flint arb (ball arithmetic, so the printed value carries a
rigorous enclosure of the floating evaluation; DV because nu-grid use is sampled)."""
import sys, os
sys.dont_write_bytecode = True
from flint import arb, acb, ctx

def Kp_arb(p, nu, J, b, T, prec=None):
    nu = arb(nu)
    if prec is None:
        prec = int((0.75*float(nu.mid()) + 60)*3.33) + 64
    old = ctx.prec; ctx.prec = prec
    s = nu*nu + arb(1)/4
    C = J*(J + 1); T = arb(T); n = J - b
    g0, gm = arb(1), arb(0); tot = arb(0); U = arb(1)
    for k in range(n + 1):
        tot += U*g0
        g1 = ((2*k*k + 2*p*k + p - s)*g0 - (k + p - 1)**2*gm)/((k + 1)**2)
        gm, g0 = g0, g1
        u = (b + k)*(b + k + 1); U = U*(C - u)/(T - u)
    ctx.prec = old
    return tot

def Nfrak(a, al):
    """N_a(alpha) = (2a)^alpha B(alpha, a+1) 2F1(a, alpha; a+alpha+1; -1)/Gamma(alpha)  (acb)"""
    al = acb(al); a = acb(a)
    B = (al.lgamma() + (a + 1).lgamma() - (al + a + 1).lgamma()).exp()
    return (2*a)**al*B*acb(-1).hypgeom_2f1(a, al, a + al + 1)*al.rgamma()

def Bcon(nu):
    nu = acb(0, nu)
    return nu.gamma()/(arb.pi().sqrt()*(nu + acb(0.5)).gamma())*acb(4)**(nu - acb(0.5))


def Kp_mid(p, nu, J, b, T, prec=None, kmax=None):
    """[DV] same sum with radii discarded every step (ball radii of the forward recurrence blow up like 3^k although the true
    rounding error is benign in the oscillatory regime).  Cross-checked against mpmath (f2_num) at J <= 2000."""
    nu = arb(nu)
    if prec is None:
        prec = int((0.75*float(nu.mid()) + 60)*3.33) + 64
    old = ctx.prec; ctx.prec = prec
    s = (nu*nu + arb(1)/4).mid()
    C = J*(J + 1); T = arb(T).mid(); n = J - b
    if kmax is None: kmax = n
    g0, gm = arb(1), arb(0); tot = arb(0); U = arb(1)
    for k in range(min(n, kmax) + 1):
        tot = (tot + U*g0).mid()
        g1 = (((2*k*k + 2*p*k + p - s)*g0 - (k + p - 1)**2*gm)/((k + 1)**2)).mid()
        gm, g0 = g0, g1
        u = (b + k)*(b + k + 1); U = (U*(C - u)/(T - u)).mid()
    ctx.prec = old
    return tot
