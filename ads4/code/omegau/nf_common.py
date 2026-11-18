"""Common tools for OMEGA_NF.md: Omega^{NF,SC}, the AdS uplift of the forward-free functional NF02 (NOFORWARD_FUNCTIONAL.md) in the B' basis,
with the superconvergent replacements of SUPERCONVERGENT_UPLIFT.md (C^-_{3,5,7} -> Psi_k/N_k, C^+_{4,6} -> Psi_k/N^ff_k, C^+_2 with S_ff).

Omega^{NF,SC} = int_0^1 dp sum_{k=2..7} phi_k(p) eps_k Lam^{2k} Chat_{k,Lam p}      (no forward terms: a = 0)
  phi_k from psi2, psi3 of NF02 exactly as for Omega^U (FLAT_UNSTRIPPED.md item 5).  At the level of the leading heavy action the normalised
  replacements equal the project moments (SUPERCONVERGENT_UPLIFT.md items 2, 3, 6), so the heavy model is that of adsaction.py with a = 0.
Precision: NF02's expanded monomial coefficients reach ~6e10 (strong cancellation; NOFORWARD_FUNCTIONAL.md section 1.3).  Here
  * psi2, psi3 and phi_k are evaluated in factored form (1-p)^5 p X(p) at the nodes;
  * the Mehler-Dirichlet / J_0 transforms use moments and endpoint derivatives computed in exact rational arithmetic (ExactCosT).
"""
import os, sys, json
from fractions import Fraction as Fr
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from common import KS, EPS, SIGMA, CH, glang, nth, Einvsqrt, Dh_all, ek, state    # noqa
from math import factorial

JS = os.path.join(HERE, '..', 'robust', 'functionals_noforward.json')

def _pmul(a, b):
    r = [Fr(0)]*(len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x*y
    return r

def _padd(*ps):
    n = max(len(p) for p in ps); r = [Fr(0)]*n
    for p in ps:
        for i, c in enumerate(p):
            r[i] += c
    return r

def _sc(c, p):
    return [c*x for x in p]

def _p2n(n):
    return [Fr(0)]*(2*n) + [Fr(1)]

class NF:
    """NF02 normalised to k(0) = 1 (exact)."""
    def __init__(self, name='NF02'):
        d = json.load(open(JS))[name]
        self.r = d['r']
        X = [Fr(c) for c in d['X']]; Y = [Fr(c) for c in d['Y']]
        s2 = [Fr(c) for c in d['psi2_exact']]; s3 = [Fr(c) for c in d['psi3_exact']]
        mu = lambda c, j: sum(ci/(i + j + 1) for i, ci in enumerate(c))
        self.k0 = -2*mu(s2, 6) + 2*mu(s3, 4)
        self.X = np.array([float(c/self.k0) for c in X]); self.Y = np.array([float(c/self.k0) for c in Y])
        self.s2 = [c/self.k0 for c in s2]; self.s3 = [c/self.k0 for c in s3]
        self.a = 0.0
        self.name = name
        # exact phi_k (monomial, rational)
        s2, s3 = self.s2, self.s3
        self.PHx = {2: s2, 3: _padd(s3, _sc(-4, _pmul(_p2n(1), s2))), 4: _padd(_sc(2, _pmul(_p2n(2), s2)), _sc(-1, _pmul(_p2n(1), s3))),
                    5: _padd(_sc(-3, _pmul(_p2n(3), s2)), _sc(2, _pmul(_p2n(2), s3))),
                    6: _padd(_sc(Fr(1, 2), _pmul(_p2n(4), s2)), _sc(Fr(-1, 2), _pmul(_p2n(3), s3))),
                    7: _padd(_sc(Fr(-1, 2), _pmul(_p2n(5), s2)), _sc(Fr(1, 2), _pmul(_p2n(4), s3)))}
        self.x0 = float(s2[1]); self.psi2f = np.array([float(c) for c in s2])
    # factored evaluation
    def psi2(self, p):
        p = np.asarray(p, float); return (1 - p)**self.r*p*np.polynomial.polynomial.polyval(p, self.X)
    def psi3(self, p):
        p = np.asarray(p, float); return (1 - p)**self.r*p*np.polynomial.polynomial.polyval(p, self.Y)
    def phi(self, k, p):
        a, b = self.psi2(p), self.psi3(p); q = np.asarray(p, float)**2
        return {2: a, 3: b - 4*q*a, 4: 2*q**2*a - q*b, 5: -3*q**3*a + 2*q**2*b, 6: 0.5*(q**4*a - q**3*b), 7: 0.5*(-q**5*a + q**4*b)}[k]

class ExactCosT:
    """int_0^1 w(p) cos(p tau) dp for an exact rational polynomial w: Taylor (exact moments) for |tau| < 8, endpoint expansion beyond
    (exact derivatives at 0 and 1).  Same algorithm as rcommon.CosT."""
    TAU0, K = 8.0, 48
    def __init__(self, c):
        c = list(c); d = len(c)
        self.d0, self.d1 = [], []
        a = c[:]
        for _ in range(d):
            self.d0.append(float(a[0]) if a else 0.0); self.d1.append(float(sum(a)))
            a = [i*a[i] for i in range(1, len(a))] or [Fr(0)]
        mom = lambda n: sum(ci/(i + n + 1) for i, ci in enumerate(c))
        self.tay = np.array([float((-1)**k*mom(2*k)/factorial(2*k)) for k in range(self.K)])
        self.n = d
    def __call__(self, tau):
        tau = np.atleast_1d(np.asarray(tau, float)); out = np.empty_like(tau)
        sm = np.abs(tau) < self.TAU0; t2 = tau[sm]**2; r = np.zeros_like(t2)
        for a in self.tay[::-1]:
            r = r*t2 + a
        out[sm] = r
        tl = tau[~sm]
        if tl.size:
            z = -1j*tl; ez = np.exp(-z); acc = np.zeros_like(z); zp = 1.0/z
            for j in range(self.n):
                acc += (self.d0[j] - ez*self.d1[j])*zp; zp = zp/z
            out[~sm] = acc.real
        return out

class TransformsNF:
    def __init__(self, F):
        self.ct = {k: ExactCosT(F.PHx[k]) for k in KS}
    def regge(self, Lam, beta):
        b = Lam*beta; th, wt = glang(nth(b)); s = np.sin(th); E = Einvsqrt(beta, s)
        return {k: (2/np.pi)*np.sum(wt*self.ct[k](b*s)*E) for k in KS}
    def flat(self, bh):
        th, wt = glang(nth(bh))
        return {k: (2/np.pi)*np.sum(wt*self.ct[k](bh*np.sin(th))) for k in KS}

def flat_parts_nf(F, h, u, Jmax, p, w):
    """stripped flat parts (B2, B3) and per-k B' parts, factored weights, on quadrature (p, w); rows J."""
    s = SIGMA[h]; y = u*p*p
    Js, D = Dh_all(h, Jmax, y)
    if h == '+-':
        Dt = D/(1 - y)**4            # dtilde44
    else:
        Dt = D
    B2 = Dt @ (w*F.psi2(p)*(2 - y)); B3 = Dt @ (w*s*u*F.psi3(p))
    parts = D @ np.array([u**(k - 2)*w*F.phi(k, p)*ek(k, y, s) for k in KS]).T
    return Js, B2, B3, parts

def envk(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2) if k % 2 == 0 else 1/(1 - x*x)**((k + 1)/2)

def action_nf(F, TR, h, st, J, mfl, Lam, pn, pw, sk, flatk, B2, B3):
    """as adsaction.action with a = 0 (no forward terms)."""
    s = SIGMA[h]; u = 1/mfl**2; beta = st['beta']; bh = (2*J + 1)/mfl
    I = TR.regge(Lam, beta); f = TR.flat(bh)
    curv = np.array([u**(k - 2)*2*EPS[k]*s**k*(I[k] - f[k]) for k in KS])
    parts = flatk + curv
    Fh = B2 + B3 + curv.sum()
    NU = np.abs(parts).sum()
    Ns = abs(B2 + curv[0]) + abs(B3 + curv[1])
    x = pn*np.sqrt(u)
    E = U = S = Es = 0.0
    for k in KS:
        rho, err, ctrl = sk[k]
        wk = u**(k - 2)*F.phi(k, pn)*pw*envk(k, x)
        S += np.sum((s**k)*wk*rho*ctrl); Es += np.sum(np.abs(wk)*err*ctrl)
        E += np.sum(np.abs(wk)*(np.abs(rho) + err)*ctrl); U += np.sum(np.abs(wk)*(~ctrl))
    return dict(Fh=Fh, Ns=Ns, NU=NU, S=S, Es=Es, E=E, U=U, beta=beta, curv=curv)
