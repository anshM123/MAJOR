"""v1: symbolic / numerical checks for H7_I14_LIGHTSIDE.md (H7: graviton-pole residue of C2, C3 on <TTTT>).

(a) [sympy, exact] heavy/light transform factor identity: CHMRS gamma_{2Df+k-1}(nu)^2 (heavy Regge moment, sec_heavy eq.
    finalreggemomentofablock) == CHP gamma(nu)gamma(-nu) (position-space <-> impact-parameter factor alpha = gamma gamma beta,
    arXiv:1707.07689 eqs. (gamma), (AlphaBeta)) at Regge spin j = k, Delta_1 = Delta_3 = Df, h = d/2, for general d.
(b) [sympy, exact] AdS shock profile (KKSDZ arXiv:1904.05905 eq. shockprofile) z^d/(z^2+y^2)^{d-1} = z * K_{d-1}, with
    K_D = (z/(z^2+y^2))^D the H_{d-1} bulk-to-boundary propagator: (Box_{H_{d-1}} - D(D-d+2)) K_D = 0, and D(D-d+2) = d-1 for D = d-1.
(c) [mpmath] harmonic transform of the H_2 propagator of dimension D0: int_1^oo Q_{D0-1}(x) P_{-1/2+i nu}(x) dx = 1/(nu^2+(D0-1/2)^2);
    D0 = 2 (graviton, d = 3) gives 1/(nu^2 + 9/4).  A double-pole (Box-2)^{-2} profile gives (nu^2+9/4)^{-2} (control).
(d) [sympy, exact] residue counting: (1/4) oint dr r^{k-2} G_+(r), G_+ = sum_{n>=-1} c_n r^n (Regge spin 2), = (2 pi i/4) c_{-1} for
    k = 2 and 0 for k = 3 (and for every k >= 3).
(e) [sympy, exact] d = 3 one-point (shock) matrices H(n) = a 1 + b (e^{4i phi} s+ + e^{-4i phi} s-) (+ parity-odd phase):
    [H(n1), H(n3)] has no a-dependence (Einstein part and Einstein x non-Einstein cross terms drop out of the bare k = 3 commutator).
(f) [numpy/scipy] regression: scalar -C_{2,nu}|_grav/(8 pi G) = 1/(nu^2+9/4) for several Df (independent referee route,
    ads4/code/referee_adslight/r1_scalar_grav.py, re-implemented here with a smaller pole sum).
(g) [sympy/mpmath] headline constants with the NF02 weights: A = 16 x0/k0, B_AdS = 16 L_fin/k0 - A log(3/2) (route-A' flat values
    3022.3, -6552.5), and the shape sensitivity d B_AdS / d log c = -A for 1/(nu^2 + c^2).
Output: cache/v1_lightside_checks.txt
"""
import os, sys, io
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs(os.path.join(HERE, 'cache'), exist_ok=True)
buf = io.StringIO()
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")

import sympy as sp
import mpmath as mp
import numpy as np

# ---------------- (a)
nu, Df, d, k = sp.symbols('nu Delta_phi d k', positive=True)
def gam_chmrs(a):   # CHMRS (gammafactors)
    return sp.gamma((1 + a - d/2 - sp.I*nu)/2)*sp.gamma((1 + a - d/2 + sp.I*nu)/2)
def gam_chp(x, J):  # CHP (gamma) with Delta_1 = Delta_3 = Df, h = d/2
    return sp.gamma((2*Df + J + sp.I*x - d/2)/2)**2
for kk in (2, 3, 4):
    lhs = gam_chmrs(2*Df + kk - 1)**2
    rhs = gam_chp(nu, kk)*gam_chp(-nu, kk)
    pr(f"(a) k = j = {kk}: gamma_{{2Df+k-1}}(nu)^2 / [gamma_CHP(nu) gamma_CHP(-nu)] = {sp.simplify(lhs/rhs)}")
lhs = gam_chmrs(2*Df + 1)**2; rhs = gam_chp(nu, 3)*gam_chp(-nu, 3)
pr(f"(a) mismatch control k = 2, j = 3: ratio = {sp.simplify(sp.gammasimp(lhs/rhs))}  (nu-dependent -> the cancellation needs k = j)")

# ---------------- (b)
z, y, D = sp.symbols('z y D', positive=True)
for dd in (3, 4, 5):
    n = dd - 2   # transverse y-dimensions of H_{d-1}: (z, y_1..y_{d-2})
    ys = sp.symbols('y1:%d' % (n + 1), real=True)
    r2 = sum(v**2 for v in ys)
    K = (z/(z**2 + r2))**(dd - 1)
    # Laplacian on H_{n+1} (upper half space, metric (dz^2+dy^2)/z^2): z^2 (d_z^2 + sum d_y^2) - (n-1) z d_z
    box = z**2*(sp.diff(K, z, 2) + sum(sp.diff(K, v, 2) for v in ys)) - (n - 1)*z*sp.diff(K, z)
    m2 = sp.simplify(box/K)
    prof = z**dd/(z**2 + r2)**(dd - 1)
    pr(f"(b) d = {dd}: Box_H{dd-1} K_{{d-1}} / K_{{d-1}} = {m2}  (= d-1 = {dd-1}); KKSDZ profile / (z K_(d-1)) = {sp.simplify(prof/(z*K))}")

# ---------------- (c)
mp.mp.dps = 20
def harm(f, nuv):
    g = lambda x: f(x)*mp.re(mp.legenp(-0.5 + 1j*nuv, 0, x, type=3))
    return mp.quad(g, [1, 1.001, 1.1, 2, 10, 100, mp.inf])
for D0 in (1, 2, 3):
    Q = lambda x, D0=D0: mp.re(mp.legenq(D0 - 1, 0, x, type=3))
    rows = []
    for nuv in (0, 0.5, 1.3, 3.0, 7.0):
        v = harm(Q, nuv); t = 1/(nuv**2 + (D0 - 0.5)**2)
        rows.append(f"nu={nuv}: {mp.nstr(v, 12)} vs {mp.nstr(t, 12)} (rel {mp.nstr(abs(v/t - 1), 2)})")
    pr(f"(c) D0 = {D0}: int Q_(D0-1) P_(-1/2+i nu) = 1/(nu^2+(D0-1/2)^2):  " + "; ".join(rows))
# control: (Box - 2)^{-2} profile = -d/d(m^2) of the propagator -> (nu^2+9/4)^{-2}; build it as -dQ_lambda/dlambda/(2 lambda+1) at lambda = 1
def Qdbl(x):
    return -mp.diff(lambda lam: mp.re(mp.legenq(lam, 0, x, type=3)), 1)/3
for nuv in (0.0, 1.3, 3.0):
    v = harm(Qdbl, nuv); t = 1/(nuv**2 + 2.25)**2
    pr(f"(c) control double pole, nu={nuv}: {mp.nstr(v, 12)} vs (nu^2+9/4)^-2 = {mp.nstr(t, 12)} (rel {mp.nstr(abs(v/t-1), 2)})")

# ---------------- (d)
r = sp.symbols('r')
cs = sp.symbols('c_m1 c_0 c_1 c_2 c_3')
G = sum(c*r**(i - 1) for i, c in enumerate(cs))
for kk in (2, 3, 4, 5):
    res = sp.Rational(1, 4)*2*sp.pi*sp.I*sp.residue(r**(kk - 2)*G, r, 0)
    pr(f"(d) k = {kk}: (1/4) oint dr r^(k-2) G_+ = {sp.simplify(res)}")

# ---------------- (e)
a, b, ph1, ph3, th = sp.symbols('a b phi1 phi3 theta', real=True)
sp_ = sp.Matrix([[0, 1], [0, 0]]); sm_ = sp.Matrix([[0, 0], [1, 0]]); I2 = sp.eye(2)
def H(ph):   # helicity basis (+2, -2); flip amplitude b e^{i theta} (theta: parity-odd phase, alpha_T)
    return a*I2 + b*(sp.exp(sp.I*(4*ph + th))*sp_ + sp.exp(-sp.I*(4*ph + th))*sm_)
Cm = sp.simplify((H(ph1)*H(ph3) - H(ph3)*H(ph1)).applyfunc(lambda e: sp.simplify(sp.expand(e.rewrite(sp.sin)))))
pr(f"(e) [H(n1), H(n3)] = {Cm.tolist()}; d/da = {sp.simplify(sp.diff(Cm, a)).tolist()}")

# ---------------- (f)
from scipy.special import loggamma, gammaln
def minusC2(Dfv, nuv, M=40000, h=0.01, Y=35.0):
    G8 = 1/(8*np.pi)
    av = 1.5 - Dfv; c = Dfv + 0.25 - 0.5j*nuv
    m = np.arange(M)
    shape = np.exp(2*gammaln(av + m) - gammaln(m + 1) - gammaln(2.5 + m))
    S0 = np.exp(2*gammaln(av) + gammaln(2.5 - 2*av) - 2*gammaln(2.5 - av))
    yv = np.arange(-Y, Y + h/2, h); zz = Dfv - 0.25 + 1j*yv; U = 2*Dfv - 2*zz
    At = np.zeros_like(zz)
    for i0 in range(0, M, 2000):
        mm = m[i0:i0 + 2000]; sh = shape[i0:i0 + 2000]
        At += (sh[None, :]/(U[:, None] - 1 - 2*mm[None, :])).sum(1)
    pexp = 2*av - 3.5; Cm_ = shape[-1]/(M - 1)**pexp
    At += -0.5*Cm_*M**pexp/(-pexp)
    kern = np.exp(2*loggamma(zz) + loggamma(c - zz) + loggamma(np.conj(c) - zz))
    I = (kern*At).sum()*h/(2*np.pi)
    pref = np.exp(gammaln(2*Dfv + 0.5) - 4*np.real(loggamma(c)))
    return pref*(-8*np.pi*G8/(4*Dfv - 1)/S0)*I
for Dfv in (1.2, 2.6):
    for nuv in (0.0, 3.0):
        v = minusC2(Dfv, nuv); t = 1/(nuv**2 + 2.25)
        pr(f"(f) Df={Dfv} nu={nuv}: -C2/(8piG) = {v.real:.10f}, 1/(nu^2+9/4) = {t:.10f}, rel {abs(v.real/t-1):.1e}, imag {abs(v.imag):.1e}")

# ---------------- (g)
sys.path.insert(0, os.path.join(HERE, '..', 'omegau'))
from nf_common import NF
F = NF(); p = sp.symbols('p', positive=True)
PH2 = sum(sp.Rational(cc.numerator, cc.denominator)*p**i for i, cc in enumerate(F.s2))   # normalised: k(0) = 1
x0 = PH2.coeff(p, 1); Lfin = sp.integrate(sp.expand((PH2 - x0*p)/p**2), (p, 0, 1))
A = 16*x0; B = 16*Lfin; BA = B - A*sp.log(sp.Rational(3, 2))
pr(f"(g) NF02 (k(0)=1): A = {float(A):.4f}, B_AdS = {float(BA):.4f}  (RESULT_D: 3022.3, -6552.5)")
# shape sensitivity: L_A(c) = int psi2(t)/(t^2 + c^2/L^2) dt ~ x0 log(L/c) + L_fin  => B_AdS(c) = B_AdS(3/2) - A log(c/(3/2))
psi2 = sp.lambdify(p, PH2, 'mpmath')
Lval = mp.mpf(98.5)   # Lambda/n at Delta_g = 100, n = 1 (Lambda = sqrt(98*99))
for cval in (mp.mpf(1.5), mp.mpf(1.5)*mp.e**0.1):
    LA = mp.quad(lambda t: psi2(t)/(t**2 + cval**2/Lval**2), [0, cval/Lval, 10*cval/Lval, 1])
    pr(f"(g) c = {mp.nstr(cval, 6)}: 16 L_A - A log L = {mp.nstr(16*LA - float(A)*mp.log(Lval), 10)}"
       f"  (asymptotic B_AdS - A log(c/1.5) = {mp.nstr(float(BA) - float(A)*mp.log(cval/1.5), 10)})")
open(os.path.join(HERE, 'cache', 'v1_lightside_checks.txt'), 'w').write(buf.getvalue())
