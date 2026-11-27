#!/usr/bin/env python3
"""nogo_cone_strong.py -- Theorem NG' : F, G not in C_D for any D <= Lambda^2/2000, every Lambda >= 998.5.
Choice rho = rho(Lambda) = (600/Lambda)^2, so x = cos(theta) = (1-rho)/(1+rho), theta = 2 arctan(600/Lambda) and
X := Lambda*theta lies in [X0, 1200], X0 = 2*998.5*arctan(600/998.5) (X increasing in Lambda), Lambda sin(theta/2) <= 600,
rho <= rho0 = (600/998.5)^2, theta <= theta0.  Cone side: log((1+rho)/(1-rho)) <= 2rho/(1-rho0), so for D <= Lambda^2/2000
the factor is <= exp(720/(1-rho0)/... ) = exp(2*3.6e5/(2000(1-rho0))).  Everything else as in nogo_cone.py, with each
quantity replaced by its worst case over X in [X0,1200] (bounds monotone in X are evaluated at the bad end).
Rigorous: exact rationals + Arb."""
import json, os
from fractions import Fraction as Fr
import sympy as sp
from flint import arb, acb, ctx
ctx.prec = 300
HERE = os.path.dirname(os.path.abspath(__file__))
nf = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'nf02_psi.json')))
K0 = sp.Rational(nf['k0'])
p = sp.symbols('p')
def q(x): x = sp.Rational(x); return arb(int(x.p))/arb(int(x.q))
L0 = arb(9985)/10
rho0 = (arb(600)/L0)**2
th0 = 2*(arb(600)/L0).atan()
X0 = L0*th0; X1 = arb(1200)
s0 = (th0/2).sin()
Pstar = 2*acb(s0*s0).elliptic_k().real/arb.pi()      # P_{-1/2}(cos theta) <= P_{-1/2}(cos theta0) (increasing in theta)
cone_exp = 2*arb(360000)/(2000*(1-rho0))               # D log((1+rho)/(1-rho)) <= this
res = {}
for name, key, Xk, rp in (('F', 'psi2_exact', 'X', 3), ('G', 'psi3_exact', 'Y', 4)):
    psi = sum(sp.Rational(c)*p**i for i, c in enumerate(nf[key]))/K0
    Xp = sp.Poly(sum(sp.Rational(c)*p**i for i, c in enumerate(nf[Xk])), p)
    assert Xp.count_roots(0, 1) == 1
    p_lo, p_hi = (sp.Rational(9168, 10000), sp.Rational(9169, 10000)) if name == 'F' else (sp.Rational(9232, 10000), sp.Rational(9233, 10000))
    assert Xp.eval(p_lo) > 0 > Xp.eval(p_hi)
    supder = sum(abs(c)*i for i, c in enumerate(reversed(sp.Poly(psi, p).all_coeffs())))
    supabs = (p_hi-p_lo)*supder
    Iplus = sp.integrate(psi, (p, 0, p_lo)) + (p_hi-p_lo)*supabs
    Iabs = sp.integrate(psi, (p, 0, p_lo)) - sp.integrate(psi, (p, p_lo, 1)) + 2*(p_hi-p_lo)*supabs
    p2, p3 = sp.Rational(95, 100), sp.Rational(99, 100)
    Iminus = -sp.integrate(psi, (p, p2, p3)); assert Iminus > 0
    P2, PH = q(p2), q(p_hi)
    # [Phi(p2 L) I- - Phi(phi L) I+] / sinh(p2 X) >= I-/(pi p2 600) - P* I+ cosh(phi X)/sinh(p2 X)
    # cosh(phi X)/sinh(p2 X) <= exp(-(p2-phi)X) (1+e^{-2 phi X})/(1-e^{-2 p2 X}), decreasing in X -> evaluate at X0
    ratio_up = (-(P2-PH)*X0).exp()*(1+(-2*PH*X0).exp())/(1-(-2*P2*X0).exp())
    lhs = q(Iminus)/(arb.pi()*P2*600) - Pstar*q(Iplus)*ratio_up
    # need (1+rho)^-rp * lhs * sinh(p2 X) > exp(cone_exp) (1-rho)^-rp Iabs ; worst case rho = rho0, X = X0
    lhs_full = lhs*(P2*X0).sinh()/(1+rho0)**rp
    rhs_full = cone_exp.exp()*q(Iabs)/(1-rho0)**rp
    ok = bool(lhs > 0) and bool(lhs_full > rhs_full)
    res[name] = dict(ok=ok, log_lhs=float(lhs_full.log().mid()), log_rhs=float(rhs_full.log().mid()))
    print(name, res[name])
json.dump(res, open(os.path.join(HERE, 'out', 'nogo_cone_strong.json'), 'w'), indent=1)
assert all(r['ok'] for r in res.values()); print("Theorem NG' PASS (D <= Lambda^2/2000, all Lambda >= 998.5)")
