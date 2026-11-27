#!/usr/bin/env python3
"""nogo_cone.py -- rigorous check of Theorem NG (no positive-cone representation of the NF02 generating functions).

Cone C_D := { sum_{k>=0, 0<=d<=D} c_{k,d} z^k (1-z)^d : c_{k,d} >= 0, convergent on |z|<1 }.
Elementary fact (Lemma NG.1): for G in C_D and 0<rho<1,  |G(-rho)| <= ((1+rho)/(1-rho))^D G(rho).
Claim: for the ideal NF02 functions
    F(z) = 2 (1-z)^-3 int_0^1 psi2(p) P_{-1/2+i Lam p}((1+z)/(1-z)) dp          (k=2 generating function)
    G(z) = Lam^2 (1-z)^-4 int_0^1 psi3(p) P_{-1/2+i Lam p}((1+z)/(1-z)) dp     (k=3 generating function)
and every real Lam >= 998.5 (i.e. every gap >= 1000) we have, at rho = 1/4,
    |F(-1/4)| > (5/3)^(Lam+1) |F(1/4)|   and   |G(-1/4)| > (5/3)^(Lam+1) |G(1/4)|,
so F, G are not in C_D for any D <= Lam+1 (in particular D = M+1 = floor(Delta_g-2)+1, or the twist-edge range a+1).

Proof ingredients (all elementary, see README Sec. NG):
 * z=-rho gives x=(1-rho)/(1+rho)=cos(theta)=3/5; Mehler-Dirichlet  Phi(nu):=P_{-1/2+i nu}(cos theta)
   = (2/pi) int_0^theta cosh(nu t) dt / sqrt(2(cos t - cos theta)) > 0, increasing in nu, and
   sinh(nu theta)/(pi nu sin(theta/2)) <= Phi(nu) <= cosh(nu theta) P_{-1/2}(cos theta),  P_{-1/2}(cos th)=(2/pi)K(sin(th/2)).
 * psi_k > 0 on (0,p*_k), psi_k < 0 on (p*_k,1) (one sign change; Sturm count in [0,1] done exactly below).
 * |P_{-1/2+i nu}(x)| <= 1 for x >= 1, so |F(rho)| <= 2(1-rho)^-3 int|psi2|, |G(rho)| <= Lam^2 (1-rho)^-4 int|psi3|.
 * int psi Phi <= Phi(p_hi Lam) int_0^{p_hi}|psi| - Phi(p2 Lam) int_{p2}^{p3}|psi|  (p_hi > p*, [p2,p3] in negative region)
   hence |int psi Phi| >= Phi(p2 Lam) I_-  - Phi(p_hi Lam) I_+.
 * the resulting inequality, divided by exp(p2 Lam theta), has LHS*Lam increasing and RHS*Lam decreasing for
   Lam >= 998.5 (exponent gaps delta=(p2-p_hi)theta > 1/998.5 and eps = p2 theta - ln(5/3) > 1/998.5), so the
   check at Lam = 998.5 (done in Arb below) proves it for every Lam >= 998.5.
"""
import json, os
from fractions import Fraction as Fr
import sympy as sp
from flint import arb, acb, ctx

ctx.prec = 200
HERE = os.path.dirname(os.path.abspath(__file__))
nf = json.load(open(os.path.join(HERE, '..', '..', 'prl', 'data', 'nf02_psi.json')))
K0 = Fr(nf['k0'])
p = sp.symbols('p')
res = {}
for name, key, Xk, lam_pow, rho_pow in (('F', 'psi2_exact', 'X', 0, 3), ('G', 'psi3_exact', 'Y', 2, 4)):
    psi = sum(sp.Rational(c)*p**i for i, c in enumerate(nf[key]))/sp.Rational(str(K0))
    Xp = sp.Poly(sum(sp.Rational(c)*p**i for i, c in enumerate(nf[Xk])), p)
    assert Xp.count_roots(0, 1) == 1          # exactly one sign change of psi in (0,1)
    # bracket the root by rationals
    p_lo, p_hi = (Fr(9168, 10000), Fr(9169, 10000)) if name == 'F' else (Fr(9232, 10000), Fr(9233, 10000))
    assert Xp.eval(sp.Rational(p_lo.numerator, p_lo.denominator)) > 0 > Xp.eval(sp.Rational(p_hi.numerator, p_hi.denominator))
    P_hi = sp.Rational(p_hi.numerator, p_hi.denominator); P_lo = sp.Rational(p_lo.numerator, p_lo.denominator)
    # I_+ >= int_0^{p*} psi = int_0^{p_hi}|psi| bound: int_0^{p_lo} psi (psi>0 there) + (p_hi-p_lo)*sup|psi| on [p_lo,p_hi]
    # on [p_lo,p_hi] |psi(p)| <= (p_hi-p_lo)*sup_[0,1]|psi'| (psi vanishes at p* in that interval)
    supder = sum(abs(c)*i for i, c in enumerate(reversed(sp.Poly(psi, p).all_coeffs())))
    supabs = (P_hi-P_lo)*supder
    Iplus = sp.integrate(psi, (p, 0, P_lo)) + (P_hi-P_lo)*supabs
    p2, p3 = sp.Rational(95, 100), sp.Rational(99, 100)
    assert all(Xp.eval(sp.Rational(k, 1000)) < 0 for k in range(950, 991))  # (sanity) negative region; proof uses 1 root + sign at p_hi
    Iminus = -sp.integrate(psi, (p, p2, p3))
    Iabs = sp.integrate(psi, (p, 0, P_lo)) - sp.integrate(psi, (p, P_lo, 1)) + 2*(P_hi-P_lo)*supabs  # >= int_0^1 |psi|
    assert Iminus > 0
    Lam = arb(9985)/10
    th = arb(3)/5
    th = th.acos()
    s2 = (arb(1)/5).sqrt()                       # sin(theta/2) = sqrt((1-3/5)/2)
    Kell = acb(s2*s2).elliptic_k().real          # K(m) with parameter m = k^2
    Pm12 = 2*Kell/arb.pi()
    def q(x):
        return arb(sp.Rational(x).p)/arb(sp.Rational(x).q)
    nu2 = q(p2)*Lam; nuh = q(P_hi)*Lam
    Phi2_lo = (nu2*th).sinh()/(arb.pi()*nu2*s2)
    Phih_hi = (nuh*th).cosh()*Pm12
    lower_int = Phi2_lo*q(Iminus) - Phih_hi*q(Iplus)
    pref_m = Lam**lam_pow/(arb(5)/4)**rho_pow * (2 if name == 'F' else 1)
    pref_p = Lam**lam_pow/(arb(3)/4)**rho_pow * (2 if name == 'F' else 1)
    lhs = pref_m*lower_int
    rhs = (arb(5)/3)**(Lam+1)*pref_p*q(Iabs)
    log_ratio = (lhs/rhs).log()
    delta = (q(p2)-q(P_hi))*th
    eps = q(p2)*th - (arb(5)/3).log()
    ok = bool(lhs > rhs) and bool(delta*Lam > 1) and bool(eps*Lam > 1) and bool(lower_int > 0)
    res[name] = dict(ok=ok, log_lhs=float(lhs.log().mid()), log_rhs=float(rhs.log().mid()), log_ratio_lower=float((log_ratio).lower()),
                     delta=float(delta.mid()), eps=float(eps.mid()), Iminus=float(Iminus), Iplus=float(Iplus), Iabs=float(Iabs))
    print(name, res[name])
json.dump(res, open(os.path.join(HERE, 'out', 'nogo_cone.json'), 'w'), indent=1)
assert all(r['ok'] for r in res.values())
print('Theorem NG check PASS at Lam=998.5 (monotone extension to all Lam>=998.5: see docstring)')
