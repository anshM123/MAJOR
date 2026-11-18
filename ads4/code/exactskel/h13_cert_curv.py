"""
h13_cert_curv: explicit (derived) remainder bound for the large-bhat envelope of the D-safe hybrid action, with every constant
computed in exact / Arb arithmetic, and the region where positivity of F_hyb follows (H13 follow-up, task 2).

Exact identities [D] (NF02_LARGE_B.md section 7):
  cosh(beta) - 1 = bhat^2/(2 Lam^2) exactly (m^2 = m2 in omegau.common.state), so M2/Lam^3 = -1/bhat^3 and the envelope's
  A2 M2/(2 Lam^3) IS the flat tail -A2/(2 bhat^3).
  F_hyb = F_flat + curv,  curv = sum_k w_k (I_k - f_k),  w_k = u^{k-2} 2 eps_k sigma^k,
  I_k - f_k = (2/pi) int_0^{pi/2} D_k(bhat sin th) dth,  D_k(s) = w(s) Phi_k(tau(s)) - Phi_k(s),
  Phi_k(x) = int_0^1 phi_k(p) cos(px) dp,  tau(s) = 2 Lam asinh(s/(2 Lam)),  w = (1 + s^2/(4 Lam^2))^{-1/2}
  (Mehler-Dirichlet for P_{-1/2+i nu}, J0(z) = (2/pi) int cos(z sin th), substitution s = 2 Lam sinh(t/2)).
  For x > 0 (finite integration by parts, phi_k polynomial, phi_k ~ (1-p)^5):
    Phi_k(x) = sum_j a_{j,k} x^{-(2j+2)} + Osc_k(x),  a_{j,k} = (-1)^{j+1} phi_k^{(2j+1)}(0),  |Osc_k(x)| <= Om_k(x) = sum_{m>=5} |phi_k^{(m)}(1)| x^{-m-1}.
  (2/pi) int_0^{pi/2} [1/s^2 - w/tau^2]_{s = bhat sin th} dth = M1(beta)/Lam^2   (checked to 15 digits, h13_asym / NF02_LARGE_B s7;
  it is the j = 0 term: a_0 = -A1).
Hence curv - A1 M1/Lam^2 = (2/pi) int_{s < s1} (D + A1 e0) dth + (2/pi) int_{s > s1} [sum_{j>=1} a_j e_j + w Osc(tau) - Osc(s)] dth,
  e_j(s) = w tau^{-(2j+2)} - s^{-(2j+2)},  and for s <= bhat <= 2 Lam (rho = s^2/(24 Lam^2) <= 1/6):
  |e_j| <= kappa_j s^{-2j}/Lam^2, kappa_j = max(1/8, ((6/5)^{2j+2} - 1)/4)   [x - x^3/6 <= asinh x <= x, 1 - y/2 <= (1+y)^{-1/2} <= 1, convexity],
  int_{th1}^{pi/2} (bhat sin th)^{-n} dth <= (pi/2)^n s1^{1-n}/((n-1) bhat)   [sin th >= 2 th/pi, th1 = asin(s1/bhat) >= s1/bhat],
  region s < s1: |.| <= (2/pi) (1/bhat) (1 - s1^2/bhat^2)^{-1/2} Z0,  Z0 >= int_0^{s1} |D + A1 e0| ds (Arb cells, Lipschitz margins).
=> |curv - A1 M1(beta)/Lam^2| <= R(bhat) = [ C_I (1 - s1^2/bhat^2)^{-1/2} + C_sm/Lam^2 + C_osc ] / bhat      (all constants >= 0, explicit)
   valid for s1 <= bhat <= 2 Lam, every m >= 1.5 Lam (|w_k| <= its value at u = 1/2.25), both helicities.
Positivity: F_hyb = F_flat + curv >= F_flat + A1_min M1(beta)/Lam^2 - R(bhat) > 0 wherever F_flat >= 0 (NOFORWARD_FUNCTIONAL s4,
  certified: all J < 400, all m >= 2 Lam, ++ any J, +- m >= 1.633 Lam) and A1_min M1 > R Lam^2 (checked here on bhat boxes, with M1
  lower bounds from h13_cert_m1-type Arb lower sums).
Output: cache/h13_cert_curv.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial, pi as PI
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 160
import h13_common as HC
import h13_cert_m1 as CM          # reuses lower_I() (Arb lower sums for I(beta) = pi M1/sqrt2 on beta boxes)

F = HC.F
KS = (2, 3, 4, 5, 6, 7)
EPS = {k: Fr(1) if k % 2 == 0 else Fr(1, 2) for k in KS}
UMAX = Fr(4, 9)
out = []
def pr(s=''):
    print(s, flush=True); out.append(s)

def A(q):
    q = Fr(q); return arb(q.numerator)/q.denominator

class Poly:
    def __init__(self, c):
        self.c = [Fr(x) for x in c]; self.d = len(self.c) - 1
        self.mom = lambda n: sum(ci/(i + n + 1) for i, ci in enumerate(self.c))
        self.der0 = [factorial(m)*self.c[m] if m <= self.d else Fr(0) for m in range(self.d + 1)]
        self.der1 = [sum(ci*Fr(factorial(i), factorial(i - m)) for i, ci in enumerate(self.c) if i >= m) for m in range(self.d + 1)]
        self.absmom = lambda n: sum(abs(ci)/(i + n + 1) for i, ci in enumerate(self.c))   # >= int p^n |phi|
        self.mu = [A(self.mom(n)) for n in range(0, 260)]
    def Phi(self, x):              # int_0^1 phi(p) cos(px) dp, x a point (arb), x >= 0
        if x < 8:
            s = arb(0); t = arb(1)
            for j in range(0, 120):
                s += ((-1)**j)*self.mu[2*j]*t
                t = t*x*x/((2*j + 1)*(2*j + 2))
            return s + arb(0, float((A(self.absmom(0))*t*(x.cosh())).upper()))
        z = acb(0, -1)*acb(x); ez = acb(0, 1)*acb(x); ez = ez.exp(); acc = acb(0); zp = 1/z
        for m in range(self.d + 1):
            acc += (A(self.der0[m]) - ez*A(self.der1[m]))*zp; zp = zp/z
        return acc.real
    def dPhi(self, x):             # Phi'(x) = -int p phi sin(px)
        if x < 8:
            s = arb(0); t = x
            for j in range(0, 120):
                s += ((-1)**j)*self.mu[2*j + 2]*t
                t = t*x*x/((2*j + 2)*(2*j + 3))
            return -s + arb(0, float((A(self.absmom(2))*t*(x.cosh())).upper()))
        chi = Poly([Fr(0)] + self.c)          # p phi
        z = acb(0, -1)*acb(x); ez = (acb(0, 1)*acb(x)).exp(); acc = acb(0); zp = 1/z
        for m in range(chi.d + 1):
            acc += (A(chi.der0[m]) - ez*A(chi.der1[m]))*zp; zp = zp/z
        return -acc.imag

PH = {k: Poly(F.PHx[k]) for k in KS}
for k in KS:
    assert all(PH[k].der1[m] == 0 for m in range(5)), "phi_k must vanish like (1-p)^5"
WMAX = {k: UMAX**(k - 2)*2*EPS[k] for k in KS}                  # max |w_k| over u <= 1/2.25
A1 = {k: PH[k].der0[1] for k in KS}                             # phi_k'(0): only k = 2, 3 nonzero
A1min = min(sum(u**(k - 2)*2*EPS[k]*(sg**k)*A1[k] for k in KS) for u in (Fr(0), UMAX) for sg in (1, -1))
pr(f"# A1 = sum_k w_k phi_k'(0): min over u in [0, 4/9], sigma = +-1: {float(A1min):.4f} (linear in u)")

def kappa(j):
    return max(Fr(1, 8), (Fr(6, 5)**(2*j + 2) - 1)/4)

def constants(Lam, s1, h=Fr(1, 100)):
    """C_sm, C_osc (u-maximised sums over k) and C_I (from Z0) for given Lam, s1."""
    hp = A(Fr(PI).limit_denominator(10**12)) if False else arb.pi()
    Csm = arb(0); Cosc = arb(0)
    for k in KS:
        P = PH[k]
        for j in range(1, P.d//2 + 1):
            aj = abs(P.der0[2*j + 1]) if 2*j + 1 <= P.d else 0
            if aj:
                Csm += A(WMAX[k]*aj*kappa(j))*(2/hp)*(hp/2)**(2*j)*arb(s1)**(1 - 2*j)/(2*j - 1)
        for m in range(5, P.d + 1):
            dm = abs(P.der1[m])
            if dm:
                Cosc += A(WMAX[k]*dm)*((A(Fr(6, 5))**(m + 1)) + 1)*(2/hp)*(hp/2)**(m + 1)*arb(s1)**(-m)/m
    # Z0 = int_0^{s1} sum_k |w_k| |D_k + a_{0,k} e0| ds, cells of width h, Lipschitz margins
    L = arb(Lam); Z0 = arb(0)
    m1 = {k: A(PH[k].absmom(1)) for k in KS}; m2 = {k: A(PH[k].absmom(2)) for k in KS}
    n = int(Fr(s1)/h)
    for i in range(n):
        sa, sb = A(i*h), A((i + 1)*h); sc = (sa + sb)/2; hh = A(h)
        # w - 1 in [w(sb) - 1, 0]; tau - s in [tau(sb) - sb, 0] (both monotone); |e0| <= kappa0/Lam^2
        wm1 = abs(1/(1 + sb*sb/(4*L*L)).sqrt() - 1)
        tsb = abs(2*L*(sb/(2*L)).asinh() - sb)
        e0b = max(Fr(1, 8), Fr(1, 12)*(1 - Fr(1, 6))**-3)
        cell = arb(0)
        for k in KS:
            P = PH[k]
            # |Phi(tau)| <= |Phi(sc)| + m1 (h + |tau - s|),  |Phi'| on [tau, s] <= |Phi'(sc)| + m2 (h + |tau - s|)
            Pb = abs(P.Phi(sc)) + m1[k]*(hh + tsb)
            dPb = abs(P.dPhi(sc)) + m2[k]*(hh + tsb)
            Db = wm1*Pb + tsb*dPb
            cell += A(WMAX[k])*(Db + abs(A(P.der0[1]))*A(e0b)/(L*L))
        Z0 += arb(cell.upper())*hh
    CI = 2/arb.pi()*Z0
    return CI, Csm, Cosc

def M1_lower(ba, bb, nv=1500):
    return (arb(2).sqrt()/arb.pi())*CM.lower_I(ba, bb, nv)

for Dg in (100, 300, 1000):
    Lam = HC.Lam_of(Dg); LamF = Fr(Lam).limit_denominator(10**9)
    best = None
    for s1 in (10, 15, 20, 30, 40, 60):
        if s1 > 2*Lam:
            continue
        CI, Csm, Cosc = constants(LamF, s1)
        # b0: smallest bhat (grid) such that for all bhat in [b0, 2 Lam] the bound holds (check on boxes)
        L2 = arb(LamF)**2
        boxes = []; b = Fr(max(s1 + 1, 10)); ok_from = None; okall = True
        bend = 2*LamF
        while b < bend:
            bb = min(b*Fr(103, 100), bend)
            ba_beta = 2*(A(b)/(2*arb(LamF))).asinh(); bb_beta = 2*(A(bb)/(2*arb(LamF))).asinh()
            Rb = (CI*(1 - arb(s1)**2/A(b)**2).rsqrt() + Csm/L2 + Cosc)/A(b)
            boxes.append((b, bb, ba_beta, bb_beta, Rb))
            b = bb
        best = best or None
        pr(f"Dg={Dg} Lam={float(Lam):.2f} s1={s1}: C_I = {float(CI.upper()):.4e}, C_sm = {float(Csm.upper()):.4e}, C_osc = {float(Cosc.upper()):.4e};"
           f"  R(bhat)*Lam^2*bhat <= {float((CI*L2 + Csm + Cosc*L2).upper()):.4e}")
        # find b0 by scanning boxes from the top (M1 is decreasing in beta, R decreasing in bhat)
        b0 = None
        for (ba, bb, bt_a, bt_b, Rb) in reversed(boxes):
            bl, bh = Fr(float(bt_a.lower())).limit_denominator(10**8), Fr(float(bt_b.upper())).limit_denominator(10**8) + Fr(1, 10**8)
            m1lo = M1_lower(bl, bh, 600)
            lhs = A(A1min)*m1lo/(arb(LamF)**2)
            if lhs > Rb:
                b0 = ba
            else:
                break
        pr(f"    => curv >= A1 M1/Lam^2 - R(bhat) > 0 on bhat in [{float(b0) if b0 else None}, 2 Lam = {float(2*Lam):.1f}]")
        if b0 is not None and (best is None or b0 < best[0]):
            best = (b0, s1)
    pr(f"## Dg={Dg}: best b0 = {float(best[0]) if best else None} (s1 = {best[1] if best else None});"
       f" positivity of F_hyb certified for bhat in [b0, 2 Lam] where the flat certificate applies")
open(os.path.join(HERE, 'cache', 'h13_cert_curv.txt'), 'w').write('\n'.join(out) + '\n')
