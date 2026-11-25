"""t0check/chk_eval.py: an independent float evaluator of the hybrid H and N_s (T0_CERT §1.1), for hostile scans [DV].
  flat:  B2 = int (2 - u p^2) psi2 D_h dp, B3 = sigma u int psi3 D_h dp, D_h = P_J(1 - 2u p^2) (++) or P^(0,8)_{J-4}(1 - 2u p^2) (+-),
         three-term recurrence at Gauss-Legendre nodes (n ~ 2 bhat + 300 nodes; the integrand has ~ bhat/pi oscillations).
  curv:  C_k = (2/pi) int_0^{pi/2} [w Phi(tau) - Phi(x)] dphi (Mehler-Dirichlet in Lemma C variables), x = bhat sin(phi),
         tau = 2 Lam asinh(x/(2 Lam)), w = (1 + x^2/(4 Lam^2))^{-1/2}, valid for every beta; Phi(x) = int psi cos(px) dp by
         Gauss-Legendre (x <= 60) or the exact finite integration-by-parts sum (x > 60); integration over x-panels
         (unit panels to 300, geometric beyond) on [0, bhat/2] and phi-panels on [pi/6, pi/2].
Own code; only the exact psi coefficients are read (prl/data/nf02_psi.json via t0_common.PSI).
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from math import factorial
import t0_common as T

PSI = {k: np.array([float(c) for c in T.PSI[k]]) for k in (2, 3)}
SIG = {'++': 1, '+-': -1}
# exact derivatives at 0 and 1 (from Fractions, then float)
D0 = {}; D1 = {}
for k in (2, 3):
    c = T.PSI[k]; d0 = []; d1 = []; cc = list(c)
    for j in range(len(c)):
        d0.append(float(cc[0])); d1.append(float(sum(cc)))
        cc = [i*cc[i] for i in range(1, len(cc))] or [0]
    D0[k] = d0; D1[k] = d1
_GL = {}
def gl(n, a=0.0, b=1.0):
    if n not in _GL:
        _GL[n] = np.polynomial.legendre.leggauss(n)
    x, w = _GL[n]
    return (b - a)/2*x + (a + b)/2, w*(b - a)/2

import json as _json
_d = _json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'nf02_psi.json')))
_K0 = float(eval(str(_d['k0'])) if not isinstance(_d['k0'], str) else float(_d['k0']))
XY = {2: np.array([float(eval(c)) for c in _d['X']]), 3: np.array([float(eval(c)) for c in _d['Y']])}
def psi(k, p):
    # factored form (1-p)^5 p X_k(p)/k0: no cancellation near p = 1 (the +- backward peak ~ J^8 amplifies any round-off there)
    return (1 - p)**5*p*np.polyval(XY[k][::-1], p)/_K0

PG, WG = gl(240)
PSG = {k: psi(k, PG)*WG for k in (2, 3)}

def Phi(k, x):
    x = np.asarray(x, float); out = np.empty_like(x)
    sm = x <= 60
    if sm.any():
        out[sm] = np.cos(np.outer(x[sm], PG)) @ PSG[k]
    if (~sm).any():
        xl = x[~sm]
        # int_0^1 psi e^{ipx} = sum_j (-1)^j [psi^(j)(1) e^{ix} - psi^(j)(0)]/(ix)^(j+1)
        s = np.zeros_like(xl, dtype=complex); e = np.exp(1j*xl)
        for j in range(len(D0[k])):
            s += (-1)**j*(D1[k][j]*e - D0[k][j])/(1j*xl)**(j + 1)
        out[~sm] = s.real
    return out

def curv(k, bh, Lam):
    """C_k(bhat, Lam) = int psi_k [P_{-1/2+i Lam p}(cosh beta) - J0(bhat p)] dp."""
    def integrand_x(x):        # [w Phi(tau) - Phi(x)]
        tau = 2*Lam*np.arcsinh(x/(2*Lam)); w = 1/np.sqrt(1 + (x/(2*Lam))**2)
        return w*Phi(k, tau) - Phi(k, x)
    tot = 0.0
    # part 1: x in [0, bhat/2], dphi = dx / sqrt(bhat^2 - x^2)
    xe = bhat_half = bh/2
    br = list(np.arange(0, min(300.0, xe), 1.0))
    if xe > 300:
        g = 300.0
        while g < xe:
            br.append(g); g *= 1.03
    br.append(xe); br = np.array(sorted(set(br)))
    for a, b in zip(br[:-1], br[1:]):
        xs, ws = gl(10, a, b)
        tot += np.sum(ws*integrand_x(xs)/np.sqrt(bh*bh - xs*xs))
    # part 2: phi in [pi/6, pi/2]
    nph = int(min(4*bh, 40000)) + 64
    npan = max(1, nph//32)
    edges = np.linspace(np.pi/6, np.pi/2, npan + 1)
    for a, b in zip(edges[:-1], edges[1:]):
        ph, wp = gl(32, a, b)
        tot += np.sum(wp*integrand_x(bh*np.sin(ph)))
    return 2/np.pi*tot

def jac_rec(n, a, b, x):
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

def flat(h, J, u):
    bh = (2*J + 1)*np.sqrt(u)
    n = int(2*bh) + 300
    p, w = gl(n)
    c = 1 - 2*u*p*p
    D = jac_rec(J, 0, 0, c) if h == '++' else jac_rec(J - 4, 0, 8, c)
    B2 = np.sum(w*(2 - u*p*p)*psi(2, p)*D); B3 = SIG[h]*u*np.sum(w*psi(3, p)*D)
    return B2, B3

def lam(Dg):
    return np.sqrt((Dg - 2.0)*(Dg - 1.0))

def hyb(Dg, J, mfl, h):
    """mfl = m/Lam."""
    L = lam(Dg); u = 1/mfl**2; bh = (2*J + 1)/mfl
    B2, B3 = flat(h, J, u)
    c2 = 2*curv(2, bh, L); c3 = u*SIG[h]*curv(3, bh, L)
    Ns = abs(B2 + c2) + abs(B3 + c3)
    return dict(base=B2 + B3 + c2 + c3, Ns=Ns, B2=B2, B3=B3, cv2=c2, cv3=c3, bhat=bh, u=u)
