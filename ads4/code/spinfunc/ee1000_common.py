"""
ee1000_common (EE_DG1000_LARGEJ.md): shared helpers for the Delta_g = 1000 large-J near-gap completion of Omega^EE.
  * imports w2_chain (EE_NEARGAP.md; all of its loader / basis / Hampel patches) read-only, in EE mode;
  * hyb(Dg, J, mfl, h, npn): the hybrid-state pieces of q3_analysis.neargap (flat NF02 + leading-Regge curvature, scalar remainder 0),
    on a Gauss-Legendre p-grid of npn nodes (q3 uses 64; large J needs more, checked in ee1000_chain);
  * mgap(Dg, J): the lowest twist-allowed mass (tau = Delta - J >= Delta_g, alpha = Delta - J - 2 >= Delta_g - 2);
  * sens(...): weight of Gamma_{k,h}(x) in dF_hh (per unit Gamma, per state), for tolerance estimates.
Nothing in the project is modified; caches are written only by the ee1000_* drivers (cache/ee1000_*).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
os.environ.setdefault('L1_DEG', '6')
os.environ.setdefault('W2_BASIS', 'ee')
os.environ.setdefault('W2_MODES', 'ee')
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import w2_chain as W
Q = W.Q
C = os.path.join(HERE, 'cache')
W.set_mode(True)


# ------------------------------------------------------------------------------------------------ Jacobi fix (EE_DG1000_LARGEJ §1)
# l9_analysis.Dc / p4_odd_analysis.Dc evaluate the +- LO shape ((1+c)/2)^4 P^(0,8)_{J-4}(c) with scipy.special.eval_jacobi; the local
# de-shift (q3._shift_generic) calls it at COMPLEX c = 1 - 2(x +- i/m)^2.  scipy's complex eval_jacobi (hypergeometric series) loses all
# accuracy for n = J - 4 >~ 24 at x >~ 0.5 (e.g. n = 36, x = 0.8: 2.4e6 instead of 0.98): this, not the K2 Laplace series, is the
# "K2 +- breakdown" of EE_NEARGAP §5 / S1_S4 at J >= 32.  Replacement: the standard three-term recurrence (stable for complex c near
# [-1, 1]), checked against mpmath below (ee1000_checks.py jacobi).
def jacobi_rec(n, a, b, x):
    x = np.asarray(x, complex)
    p0 = np.ones_like(x)
    if n == 0:
        return p0
    p1 = 0.5*(a - b + (a + b + 2)*x)
    for k in range(2, n + 1):
        c1 = 2*k*(k + a + b)*(2*k + a + b - 2)
        c2 = (2*k + a + b - 1)*((2*k + a + b)*(2*k + a + b - 2)*x + a*a - b*b)
        c3 = 2*(k + a - 1)*(k + b - 1)*(2*k + a + b)
        p0, p1 = p1, (c2*p1 - c3*p0)/c1
    return p1


def legendre_rec(n, x):
    return jacobi_rec(n, 0, 0, x)


def Dc_fixed(h, J, x):
    c = 1 - 2*np.asarray(x)**2
    if h == '++':
        v = legendre_rec(J, c)
    else:
        v = ((1 + c)/2)**4*jacobi_rec(J - 4, 0, 8, c)
    return v if np.iscomplexobj(x) else v.real


_orig_Dc = Q.A.Dc
JACFIX = os.environ.get('EE1000_JACFIX', '1') == '1'
if JACFIX:
    Q.A.Dc = Dc_fixed
    Q.P4.Dc = Dc_fixed
    Q.A.shift = Q._shift_generic(Q.A.env, Dc_fixed, Q.A.SIG)
    Q.P4.shift = Q._shift_generic(Q.P4.env, Dc_fixed, Q.P4.SIG)


def lam_of(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))


def mphys(Delta, J):
    al = Delta - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def mgap(Dg, J):
    """lowest twist-allowed m/Lam at spin J (alpha = Dg - 2)."""
    al = Dg - 2.0
    return np.sqrt(al*(al + 2*J + 1))/lam_of(Dg)


def hyb(Dg, J, mfl, h, npn=64):
    """hybrid near-gap state pieces exactly as q3_analysis.neargap (DF tag 'hyb'), on gl01(npn)."""
    Lm = lam_of(Dg)
    pn, pw = Q.gl01(npn)
    F = Q.NearGap(Dg).F
    u = 1/mfl**2; s = Q.SIGMA[h]
    x = pn*np.sqrt(u); y = x*x
    eta = 1.0 + ((2*J + 1)/mfl)**2/(2*Lm*Lm)
    P = Q.remk.Pders_vec(Lm*pn, eta, nmax=0)[0]; jz = Q.j0((2*J + 1)*x)
    pf2, wf2 = Q.gl01(max(400, 2*J + 100))
    Js, B2, B3, _ = Q.flat_parts_nf(F, h, u, J, pf2, wf2); i = int(np.searchsorted(Js, J)); B2, B3 = B2[i], B3[i]
    p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
    cv2 = (p2*2*(P - jz)).sum(); cv3 = (p3*u*s*(P - jz)).sum()
    Ns = abs(B2 + cv2) + abs(B3 + cv3)
    return dict(base=B2 + B3 + cv2 + cv3, Ns=Ns, B2=B2, B3=B3, cv2=cv2, cv3=cv3, x=x, y=y, p2=p2, p3=p3, u=u, m=mfl*Lm, pn=pn)


def dF_state(st, G, h):
    """diagonal correction of state st (from hyb) for helicity h, given G = {('C2',h): Gamma(x nodes), ('C3H',h): ...}."""
    x, m = st['x'], st['m']
    d = (st['p2']*Q.A.env(h, x, 2)*G[('C2', h)]/m**2).sum()
    d += (st['p3']*st['u']*Q.A.env(h, x, 3)*G[('C3H', h)]/m**2).sum()
    return d


def dF_off(stp, stm, G):
    """off-diagonal correction (c_off = m Gamma_off model already divided by m in G[('C2','off')], as q3 OFFSCALE)."""
    x, m, y, u = stp['x'], stp['m'], stp['y'], stp['u']
    eo2 = np.sqrt(Q.A.env('++', x, 2)*Q.A.env('+-', x, 2))
    return (stp['p2']*eo2*G[('C2', 'off')]).sum()/m + (stp['p3']*u*eo2/(2 - y)*G[('C3H', 'off')]).sum()/m


def lam2(base, Ns, dF):
    Fm = np.diag(base) + dF; n = 1/np.sqrt(Ns)
    return np.linalg.eigvalsh(n[:, None]*Fm*n[None, :])[0]
