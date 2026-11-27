#!/usr/bin/env python3
"""hev.py -- high-precision evaluator of the exact NF02 heavy angular action (uncut ideal family) at any gap.

S := H_ang,b / w^3 = sum_{k=0}^{n} U_k [ a_k + ( gam_k a_k + (k+1)^2 a_{k+1}/2 + sigma_b Lam^2 g3_k ) / T ]

with U_k = prod_{j<k} (C-u_j)/(T-u_j), u_j=(b+j)(b+j+1), gam_k = (-k^2+(8-4b)k+9-2b(b-1))/2,
a_k = 2 int_0^1 psi2(p) g[k,3](Lam^2 p^2) dp, g3_k = int_0^1 psi3(p) g[k,4](Lam^2 p^2) dp,
g[k,p](y) = coefficient of z^k in (1-z)^-p P_{-1/2+i nu}((1+z)/(1-z)), y = nu^2 (three-term recurrence).
This is the research note's Sec. 3 closed form rewritten with u^k t^k = m^{-2k} (identical algebra; validated
against rv_closed.py exact rationals in validate_hev.py).  w = m^2/T, H = w^3 S.

The moments are computed by Gauss-Legendre quadrature in p (exact for the polynomial integrand when the node
number is >= K/1+9) in Arb ball arithmetic; the ladder sum is truncated at K (numerical diagnostic: the last
terms are reported).  Values are rigorous enclosures of the TRUNCATED sum only; truncation is NOT certified.
Label: numerical evidence.
"""
import json, os, sys, math
from fractions import Fraction as Fr
import flint
from flint import arb, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, '..', '..', 'prl', 'data'))
_nf = json.load(open(os.path.join(DATA, 'nf02_psi.json')))
K0 = Fr(_nf['k0'])
PSI2 = [Fr(c)/K0 for c in _nf['psi2_exact']]
PSI3 = [Fr(c)/K0 for c in _nf['psi3_exact']]


def _arbq(q):
    return arb(q.numerator)/arb(q.denominator)


class Coeffs:
    """state-independent spectral coefficients a_k, g3_k at a given Lam^2 (exact rational or float string)."""

    def __init__(self, lam2, K, prec=None, verbose=False):
        lam2 = Fr(lam2)
        self.lam2 = lam2
        self.K = K
        Lam = math.sqrt(float(lam2))
        if prec is None:
            prec = int(2.5*Lam/math.log(2)) + 320
        self.prec = prec
        ctx.prec = prec
        L2 = _arbq(lam2)
        nn = K + 12
        nodes = [arb.legendre_p_root(nn, i, weight=True) for i in range(nn)]
        ps = [(x+1)/2 for x, w in nodes]
        ws = [w/2 for x, w in nodes]
        c2 = [_arbq(c) for c in PSI2]
        c3 = [_arbq(c) for c in PSI3]
        def poly(c, p):
            s = arb(0)
            for x in reversed(c):
                s = s*p + x
            return s
        w2 = [w*poly(c2, p) for p, w in zip(ps, ws)]
        w3 = [w*poly(c3, p) for p, w in zip(ps, ws)]
        ys = [L2*p*p for p in ps]
        a = [arb(0)]*(K+2)
        g3 = [arb(0)]*(K+2)
        for (pp, wt, kk) in ((3, w2, 'a'), (4, w3, 'g')):
            out = a if kk == 'a' else g3
            for i in range(nn):
                y = ys[i]
                gm, g0 = arb(0), arb(1)
                wi = wt[i]
                out[0] += wi*g0
                for k in range(K+1):
                    g1 = ((8*k*k+8*pp*k+4*pp-1 - 4*y)*g0 - 4*(k-1+pp)**2*gm)/(4*(k+1)**2)
                    gm, g0 = g0, g1
                    out[k+1] += wi*g0
            if verbose:
                print('done p=', pp, file=sys.stderr)
        self.a = [2*x for x in a]
        self.g3 = g3
        self.L2 = L2


def S_state(cf, J, b, a_twist, Kmax=None, want_terms=False):
    """S = H/w^3 for spin J, ladder b, twist a (real; Delta = a+J+2). Returns (S, w, info)."""
    ctx.prec = cf.prec
    a_twist = Fr(a_twist)
    Delta = a_twist + J + 2
    T = _arbq((Delta-1)*(Delta-2))
    C = J*(J+1)
    n = J - b
    sig = 1 if b == 0 else -1
    K = cf.K if Kmax is None else min(Kmax, cf.K)
    kmax = min(n, K)
    S = arb(0)
    U = arb(1)
    last = []
    for k in range(kmax+1):
        gam = arb(-k*k+(8-4*b)*k+9-2*b*(b-1))/2
        term = U*(cf.a[k] + (gam*cf.a[k] + arb((k+1)**2)/2*cf.a[k+1] + sig*cf.L2*cf.g3[k])/T)
        S += term
        if k >= kmax-3:
            last.append(term)
        uj = (b+k)*(b+k+1)
        U = U*(C-uj)/(T-uj)
    m2 = T - C
    w = m2/T
    info = {'truncated': n > K, 'last_terms': [float(x.mid()) for x in last], 'kmax': kmax}
    return S, w, info


def H_LR(cf, J, b, a_twist):
    """leading-Regge proxy (1-r)^3 sum r^k [a_k + sigma Lam^2 g3_k / T]  = 2 int psi2 P(eta) + sigma u int psi3 P(eta)."""
    ctx.prec = cf.prec
    a_twist = Fr(a_twist)
    Delta = a_twist + J + 2
    T = _arbq((Delta-1)*(Delta-2))
    C = J*(J+1)
    r = arb(C)/T
    sig = 1 if b == 0 else -1
    s = arb(0)
    rk = arb(1)
    for k in range(cf.K+1):
        s += rk*(cf.a[k] + sig*cf.L2*cf.g3[k]/T)
        rk *= r
    return (1-r)**3*s


def state_params(lam2, J, a_twist):
    lam = math.sqrt(float(lam2))
    Delta = float(a_twist) + J + 2
    T = (Delta-1)*(Delta-2)
    C = J*(J+1)
    m = math.sqrt(T-C)
    bhat = (2*J+1)*lam/m
    ell = math.log(T/(T-C))
    return {'bhat': bhat, 'ell': ell, 'm': m, 'T': T, 'r': C/T}


if __name__ == '__main__':
    cf = Coeffs(997002, 200, verbose=True)
    for J in (10, 40, 100):
        S, w, info = S_state(cf, J, 0, 998)
        print(J, S.str(20), info)
