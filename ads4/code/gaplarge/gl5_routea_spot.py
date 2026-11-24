"""
gaplarge/gl5_routea_spot.py (GAP_LARGESPIN.md §3.2): VERBATIM COPY of exactskel/h13a_cert_curv.py (ROUTEA §4.2/4.4 certified window
curv > 0 on bhat in [b0, 2 Lam], m >= 1.5 Lam), with three changes only: (i) sys.dont_write_bytecode; (ii) S1S has an entry for every
requested Dg (H13A_DGS) with s1 in H13A_S1 (default 60,80,120); (iii) output to gaplarge/cache/gl5_routea_<Dgs>.txt.
Read-only imports from exactskel/: h13_common, h13_cert_m1.  Run: H13A_SIGNED=1 H13A_DGS=2000,5000 python3 gl5_routea_spot.py
--- original docstring follows ---
h13a_cert_curv: ROUTE A-prime version (only C^_2, C^_3: weights w_2 = 2, w_3 = u sigma on psi_2, psi_3) of the explicit (derived) remainder bound
for the large-bhat envelope of the D-safe hybrid action, with every constant
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
sys.dont_write_bytecode = True
MYDIR = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.abspath(os.path.join(MYDIR, '..', 'exactskel')); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from math import factorial, pi as PI
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 160
import h13_common as HC
import h13_cert_m1 as CM          # reuses lower_I() (Arb lower sums for I(beta) = pi M1/sqrt2 on beta boxes)

F = HC.F
KS = (2, 3)
EPS = {2: Fr(1), 3: Fr(1, 2)}
SIGNED_PRE = bool(int(os.environ.get('H13A_SIGNED', '0')))
UMAX = Fr(4, 9)
DGS = [int(a) for a in os.environ.get('H13A_DGS', '100,300,1000').split(',')]
S1S = {100: (50, 80, 120), 300: (80, 120, 160, 200), 1000: (120, 160, 200, 300, 400)}
if SIGNED_PRE:
    S1S = {100: (20, 30, 40, 50, 80), 300: (30, 40, 60, 80, 120), 1000: (40, 60, 80, 120, 160)}
_S1 = tuple(int(a) for a in os.environ.get('H13A_S1', '60,80,120').split(','))
for _dg in DGS:
    S1S[_dg] = _S1
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

PH = {2: Poly(F.s2), 3: Poly(F.s3)}
XY = {2: F.X, 3: F.Y}
from h13_common import F as _F
def absmoments(k, nmax=3, ncell=4000):
    # rigorous upper bounds of int_0^1 p^n |psi_k(p)| dp, psi_k = (1-p)^r p X(p) (factored, exact rationals), Arb cells
    from fractions import Fraction as _Fr
    import json
    coef = [_Fr(c) for c in json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_noforward.json')))['NF02']['X' if k == 2 else 'Y']]
    coef = [c/_F.k0 for c in coef]
    r = _F.r
    res = [arb(0)]*(nmax + 1)
    for i in range(ncell):
        p = CM.iv(_Fr(i, ncell), _Fr(i + 1, ncell))
        X = arb(0)
        for c in reversed(coef):
            X = X*p + A(c)
        v = abs((1 - p)**r*p*X)
        for n in range(nmax + 1):
            res[n] = res[n] + arb((v*p**n).upper())/ncell
    return [arb(x.upper()) for x in res]
ABSM = {k: absmoments(k) for k in KS}
for k in KS:
    assert all(PH[k].der1[m] == 0 for m in range(5)), "phi_k must vanish like (1-p)^5"
WMAX = {k: UMAX**(k - 2)*2*EPS[k] for k in KS}                  # max |w_k| over u <= 1/2.25
A1 = {k: PH[k].der0[1] for k in KS}                             # phi_k'(0): only k = 2, 3 nonzero
A1min = min(sum(u**(k - 2)*2*EPS[k]*(sg**k)*A1[k] for k in KS) for u in (Fr(0), UMAX) for sg in (1, -1))
pr(f"# rigorous int p^n |psi_k| (n=0..3): " + "; ".join(f"k={k}: " + ", ".join(f"{float(x.upper()):.3g}" for x in ABSM[k]) for k in KS))
pr(f"# A1 = sum_k w_k phi_k'(0): min over u in [0, 4/9], sigma = +-1: {float(A1min):.4f} (linear in u)")

def kappa(j):
    return max(Fr(1, 8), (Fr(6, 5)**(2*j + 2) - 1)/4)

DIAG = None
SIGNED = bool(int(os.environ.get('H13A_SIGNED', '0')))
def constants(Lam, s1, h=Fr(1, 25)):
    """C_sm, C_osc (u-maximised sums over k) and C_I (from Z0) for given Lam, s1."""
    hp = A(Fr(PI).limit_denominator(10**12)) if False else arb.pi()
    Csm = arb(0); Cosc = arb(0); Cosc2 = arb(0)
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
                # sharper (H13 follow-up): |w Osc(tau) - Osc(s)| <= |1-w| Om(tau) + |s - tau| Om1(tau), |1-w| <= s^2/(8 Lam^2),
                # |s - tau| <= s^3/(24 Lam^2), tau >= (5/6) s, Om1 = sum |d1_m| (1 + (m+1)/tau) tau^{-m-1}; then the same theta integral
                if m >= 4:
                    t1 = A(Fr(1, 8))*(A(Fr(6, 5))**(m + 1))*(2/hp)*(hp/2)**(m - 1)*arb(s1)**(2 - m)/(m - 2)
                    t2 = A(Fr(1, 24))*(A(Fr(6, 5))**(m + 1))*(1 + A(Fr(6, 5))*(m + 1)/arb(s1))*(2/hp)*(hp/2)**(m - 2)*arb(s1)**(3 - m)/(m - 3)
                    Cosc2 += A(WMAX[k]*dm)*(t1 + t2)/A(Lam)**2
    # Z0 >= int_0^{s1} sum_k |w_k| |G_k|,  G_k = D_k + phi_k'(0) e0  (= D_k - a_{0,k} e0: the part of D not in the M1 term).
    # s < 4: |G_k| <= |D_k| + |phi_k'(0)| |e0| (crude, Lipschitz margins with the rigorous m_n);
    # s >= 4: second-order Taylor about the cell centre, |G| <= |G(sc)| + (h/2)|G'(sc)| + (h^2/8) sup|G''|  (G, G' at the point sc
    # in Arb; sup|G''| <= |w''| m0 + 3|w w'| m1 + |1 - w^3| m2 + |s - tau| m3 + |phi'(0)| sup|e0''|, e0'' by ball evaluation).
    L = A(Lam); Z0 = arb(0); L2_ = L*L
    P0 = {k: arb(0) for k in KS}; Q0 = {k: arb(0) for k in KS}      # signed int G_k (enclosure) and int |G_k| s^2 (upper)
    m = {k: ABSM[k] for k in KS}
    n = int(Fr(s1)/h); hh = A(h)
    e0b = max(Fr(1, 8), Fr(1, 12)*(1 - Fr(1, 6))**-3)
    for i in range(n):
        sa, sb = A(i*h), A((i + 1)*h); sc = (sa + sb)/2
        w_ = lambda x: 1/(1 + x*x/(4*L2_)).sqrt()
        tau_ = lambda x: 2*L*(x/(2*L)).asinh()
        cell = arb(0)
        if i*h < 4:
            wm1 = abs(w_(sb) - 1); tsb = abs(tau_(sb) - sb)
            for k in KS:
                P = PH[k]
                Pb = abs(P.Phi(sc)) + m[k][1]*(hh + tsb); dPb = abs(P.dPhi(sc)) + m[k][2]*(hh + tsb)
                gb = wm1*Pb + tsb*dPb + abs(A(P.der0[1]))*A(e0b)/L2_
                cell += A(WMAX[k])*gb
                P0[k] += arb(0, float((gb*hh).upper())); Q0[k] += arb((gb*hh*sb*sb).upper())
        else:
            wc = w_(sc); tc = tau_(sc); wpc = -sc*wc**3/(4*L2_)
            e0 = wc/tc**2 - 1/sc**2
            e0p = wpc/tc**2 - 2*wc*wc/tc**3 + 2/sc**3
            S = CM.iv(Fr(i)*h, Fr(i + 1)*h)
            wS = w_(S); tS = tau_(S); wpS = -S*wS**3/(4*L2_); wppS = -wS**3/(4*L2_) + 3*S*S*wS**5/(16*L2_*L2_)
            e0pp = abs(wppS/tS**2 - 6*wS*wpS/tS**3 + 6*wS**3/tS**4 - 6/S**4)
            Wpp = 1/(4*L2_) + 3*sb*sb/(16*L2_*L2_); Wp = sb/(4*L2_); W3 = 3*sb*sb/(8*L2_); Ts = sb**3/(24*L2_)
            for k in KS:
                P = PH[k]; a1 = A(P.der0[1])
                Dc = wc*P.Phi(tc) - P.Phi(sc)
                Dpc = wpc*P.Phi(tc) + wc*wc*P.dPhi(tc) - P.dPhi(sc)
                Gc = Dc + a1*e0; Gpc = Dpc + a1*e0p
                Gpp = Wpp*m[k][0] + 3*Wp*m[k][1] + W3*m[k][2] + Ts*m[k][3] + abs(a1)*e0pp
                cell += A(WMAX[k])*(abs(Gc) + hh/2*abs(Gpc) + hh*hh/8*Gpp)
                # midpoint rule on the cell: int_cell G = h G(sc) + e, |e| <= h^3/24 sup|G''|
                P0[k] += hh*Gc + arb(0, float((hh**3/24*Gpp).upper()))
                Q0[k] += arb(((abs(Gc) + hh/2*abs(Gpc) + hh*hh/8*Gpp)*hh*sb*sb).upper())
        Z0 += arb(cell.upper())*hh
        if DIAG is not None:
            DIAG.append((float(i*h), float((arb(cell.upper())*hh).upper())))
    CI = 2/arb.pi()*Z0
    CIP = 2/arb.pi()*sum((A(WMAX[k])*arb(abs(P0[k]).upper()) for k in KS), arb(0))
    CIQ = 2/arb.pi()*sum((A(WMAX[k])*Q0[k] for k in KS), arb(0))/2
    if SIGNED:
        return (CIP, CIQ), Csm, (Cosc2 if Cosc2 < Cosc else Cosc)
    return CI, Csm, Cosc

def M1_lower(ba, bb, nv=1500):
    return (arb(2).sqrt()/arb.pi())*CM.lower_I(ba, bb, nv)

if __name__ != '__main__':
    raise SystemExit
for Dg in DGS:
    Lam = HC.Lam_of(Dg); LamF = Fr(Lam).limit_denominator(10**9)
    best = None
    for s1 in S1S[Dg]:
        if s1 > 2*Lam:
            continue
        CI, Csm, Cosc = constants(LamF, s1)
        # b0: smallest bhat (grid) such that for all bhat in [b0, 2 Lam] the bound holds (check on boxes)
        L2 = A(LamF)**2
        boxes = []; b = Fr(max(s1 + 1, 10)); ok_from = None; okall = True
        bend = 2*LamF
        while b < bend:
            bb = min(b*Fr(103, 100), bend)
            ba_beta = 2*(A(b)/(2*A(LamF))).asinh(); bb_beta = 2*(A(bb)/(2*A(LamF))).asinh()
            if SIGNED:
                CIP, CIQ = CI
                Rb = (CIP + CIQ*(1 - arb(s1)**2/A(b)**2)**arb(-1.5)/A(b)**2 + Csm/L2 + Cosc)/A(b)
            else:
                Rb = (CI*(1 - arb(s1)**2/A(b)**2).rsqrt() + Csm/L2 + Cosc)/A(b)
            boxes.append((b, bb, ba_beta, bb_beta, Rb))
            b = bb
        best = best or None
        CIr = CI[0] if SIGNED else CI
        pr(f"Dg={Dg} Lam={float(Lam):.2f} s1={s1}: C_I = {float(CIr.upper()):.4e}" + (f" (signed), C_Q = {float(CI[1].upper()):.4e}" if SIGNED else "")
           + f", C_sm = {float(Csm.upper()):.4e}, C_osc = {float(Cosc.upper()):.4e};  (C_I Lam^2 + C_sm + C_osc Lam^2) <= {float((CIr*L2 + Csm + Cosc*L2).upper()):.4e}")
        # find b0 by scanning boxes from the top (M1 is decreasing in beta, R decreasing in bhat)
        b0 = None
        for (ba, bb, bt_a, bt_b, Rb) in reversed(boxes):
            bl, bh = Fr(float(bt_a.lower())).limit_denominator(10**8), Fr(float(bt_b.upper())).limit_denominator(10**8) + Fr(1, 10**8)
            m1lo = M1_lower(bl, bh, 600)
            lhs = A(A1min)*m1lo/(A(LamF)**2)
            if lhs > Rb:
                b0 = ba
            else:
                break
        pr(f"    => curv >= A1 M1/Lam^2 - R(bhat) > 0 on bhat in [{float(b0) if b0 else None}, 2 Lam = {float(2*Lam):.1f}]")
        if b0 is not None and (best is None or b0 < best[0]):
            best = (b0, s1)
    pr(f"## Dg={Dg}: best b0 = {float(best[0]) if best else None} (s1 = {best[1] if best else None});"
       f" positivity of F_hyb certified for bhat in [b0, 2 Lam] where the flat certificate applies")
open(os.path.join(MYDIR, 'cache', 'gl5_routea_%s%s.txt' % ('_'.join(map(str, DGS)), '_signed' if SIGNED else '')), 'w').write('\n'.join(out) + '\n')
