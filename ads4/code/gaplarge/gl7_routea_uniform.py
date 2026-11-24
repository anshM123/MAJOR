"""
gaplarge/gl7_routea_uniform.py (GAP_LARGESPIN.md §3.2): ROUTEA's certified curvature window (ROUTEA_LARGE_B §4.2/§4.4) made UNIFORM in
Lam >= Lam0 = Lam(1000), i.e. for every Delta_g >= 1000.  [D + certified]

ROUTEA's bound (m >= 1.5 Lam [env GL7_UMAX: any u = Lam^2/m^2 <= UMAX], both helicities, s1 <= bhat <= 2 Lam):
   |curv - A1 M1(beta)/Lam^2| <= R = [C_P + C_Q (1 - s1^2/bhat^2)^{-3/2}/bhat^2 + C_sm/Lam^2 + C_osc2]/bhat,
   C_P = (2/pi) sum_k |w_k| |int_0^{s1} G_k ds|,  C_Q = (1/pi) sum_k |w_k| int_0^{s1} |G_k| s^2 ds,  G_k(s) = w r_k(tau) - r_k(s),
   r_k = Phi_k + a1_k/x^2, w = (1 + s^2/4Lam^2)^{-1/2}, tau = 2 Lam asinh(s/2Lam);  C_sm and Lam^2 C_osc2 are Lam-INDEPENDENT by
   construction (|e_j| <= kappa_j s^{-2j}/Lam^2, |1-w| <= s^2/8Lam^2, |s - tau| <= s^3/24Lam^2).
New here (GAP_LARGESPIN §3.2, [D]): with eps = 1/Lam^2,
   Lam^2 G = g1 + eps h,   g1 = -(s^3/24) r'(s) - (s^2/8) r(s) = -d/ds[(s^3/24) r(s)],
   |h| <= (3 s^5/640)|r'(s)| + (3 s^4/128)|r(s)| + (s^6/1152) sup_{[5s/6,s]}|r''| + (s^5/192) sup_{[5s/6,s]}|r'|
   (alternating Taylor bounds of asinh and (1+y)^{-1/2} on [0, 1]; tau >= 5s/6), hence for every Lam >= Lam0:
   |int_0^{s1} Lam^2 G| <= s1^3 |r(s1)|/24 + eps0 int_0^{s1} hbar,     int_0^{s1} |Lam^2 G| s^2 <= int_0^{s1} (gbar + eps0 hbar) s^2,
   with the decreasing majorants |r^(n)(x)| <= R_n(x) = min(asymptotic: exact finite IBP sum with |.| ; crude: n! a1 (n+1)/... below).
Positivity for all Lam >= Lam0 on a bhat-box [b, 1.03 b] (U decreasing in bhat; beta = 2 asinh(bhat/2Lam) <= beta_max = 2 asinh(1.03b/2Lam0),
and beta <= 2 asinh(1) since bhat <= 2 Lam):  A1min_{u<=4/9} * min_{beta <= beta_max} M1(beta) > U(b) := Lam^2 R(b) (Lam-free upper bound).
M1 lower bounds: lrtight/cache/lrt2_m1.json cells on [0.1, 24] and, on (0, 0.1], the [P] bound
   M1(beta) >= (1/24 - 7 beta^2/1920)/cosh(beta/2) - 2 e^{-2 pi nu1}(nu1/(2 pi) + 1/(4 pi^2)),  nu1 = sqrt2/beta.
Output: cache/gl7_routea_uniform.{txt,json}
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, ctx
import gl_common as G
from gl_common import A, up, lo, W

ctx.prec = 128
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)

PI = arb.pi()
UMAX = Fr(os.environ.get('GL7_UMAX', '4/9'))   # u = Lam^2/m^2 <= UMAX: 4/9 (m >= 1.5 Lam, ROUTEA), 100/169 (m >= 1.3 Lam), 1 (m >= Lam)
WMAX = {2: Fr(2), 3: UMAX}               # max |w_k|: w_2 = 2, w_3 = u sigma
EPS0 = 1/G.LAM0**2


def kappa(j):
    return max(Fr(1, 8), (Fr(6, 5)**(2*j + 2) - 1)/4)


class RB:
    """majorants R_n(x) >= |r^(n)(x)|, n = 0, 1, 2, decreasing in x > 0."""

    def __init__(self, w):
        self.w = w
        n = len(w.c)
        self.nonosc = [(j, abs(w.d0[j])) for j in range(3, n, 2) if w.d0[j] != 0]       # term c/x^{j+1}
        self.osc = [(j, abs(w.d1[j])) for j in range(5, n) if w.d1[j] != 0]              # term c e^{ix}/x^{j+1}
        self.a1 = A(w.a1); self.I = A(w.absint_up)                                       # int p^n |psi| <= int |psi|

    def asym(self, x, n):
        s = arb(0)
        for j, c in self.nonosc:
            f = 1
            for i in range(1, n + 1):
                f *= (j + i)
            s += A(c)*f/x**(j + 1 + n)
        for j, c in self.osc:
            a = j + 1
            fac = arb(1) if n == 0 else (1 + a/x if n == 1 else 1 + 2*a/x + a*(a + 1)/x**2)
            s += A(c)*fac/x**a
        return s

    def crude(self, x, n):
        return self.I + [1, 2, 6][n]*self.a1/x**(2 + n)

    def R(self, x, n):
        a, b = self.asym(x, n), self.crude(x, n)
        return a if a.upper() <= b.upper() else b


def gbar(rb, s_lo, s_hi):
    """sup over the cell of (s^3/24)|r'| + (s^2/8)|r|  (R_n decreasing, powers increasing)."""
    return s_hi**3/24*rb.R(s_lo, 1) + s_hi**2/8*rb.R(s_lo, 0)


def hbar(rb, s_lo, s_hi):
    q = A(Fr(5, 6))*s_lo
    return (3*s_hi**5/640*rb.R(s_lo, 1) + 3*s_hi**4/128*rb.R(s_lo, 0) + s_hi**6/1152*rb.R(q, 2) + s_hi**5/192*rb.R(q, 1))


def first_cell(rb, s0):
    """int_0^{s0} of gbar s^2, hbar, hbar s^2 using the crude (non-singular) forms, bounded by s0 x value at s0 (integrands increasing)."""
    a1, I = rb.a1, rb.I; c = A(Fr(6, 5))
    g = a1/12 + a1/8 + I*(s0**3/24 + s0**2/8)
    h = s0**2*(6*a1/640 + 3*a1/128 + 6*a1*c**4/1152 + 2*a1*c**3/192) + I*(3*s0**5/640 + 3*s0**4/128 + s0**6/1152 + s0**5/192)
    return up(s0*g*s0**2), up(s0*h), up(s0*h*s0**2)


def region_I(k, s1):
    rb = RB(W[k])
    s0 = A(Fr(1, 1000))
    Ig2, Ih, Ih2 = first_cell(rb, s0)
    # cells: geometric 1e-3 .. 1 (ratio 1.01), then width 1/50 to s1
    edges = [Fr(1, 1000)]
    while edges[-1] < 1:
        edges.append(min(edges[-1]*Fr(101, 100), Fr(1)))
    x = Fr(1)
    while x < s1:
        x = min(x + Fr(1, 50), Fr(s1)); edges.append(x)
    for a, b in zip(edges[:-1], edges[1:]):
        sa, sb = A(a), A(b); hh = sb - sa
        gb = gbar(rb, sa, sb); hb = hbar(rb, sa, sb)
        Ig2 += up(gb*sb**2*hh); Ih += up(hb*hh); Ih2 += up(hb*sb**2*hh)
    P = up(A(s1)**3*rb.R(A(s1), 0)/24 + EPS0*Ih)
    Q = up(Ig2 + EPS0*Ih2)
    return P, Q, dict(lead=float(up(A(s1)**3*rb.R(A(s1), 0)/24).mid()), epsh=float(up(EPS0*Ih).mid()), Ig2=float(Ig2.mid()))


def lam_free_consts(s1):
    hp = PI
    Csm = arb(0); Cosc2 = arb(0)
    for k in (2, 3):
        w = W[k]; d = len(w.c) - 1
        for j in range(1, d//2 + 1):
            aj = abs(w.d0[2*j + 1]) if 2*j + 1 <= d else 0
            if aj:
                Csm += A(WMAX[k]*aj*kappa(j))*(2/hp)*(hp/2)**(2*j)*arb(s1)**(1 - 2*j)/(2*j - 1)
        for m in range(5, d + 1):
            dm = abs(w.d1[m])
            if dm:
                t1 = A(Fr(1, 8))*(A(Fr(6, 5))**(m + 1))*(2/hp)*(hp/2)**(m - 1)*arb(s1)**(2 - m)/(m - 2)
                t2 = A(Fr(1, 24))*(A(Fr(6, 5))**(m + 1))*(1 + A(Fr(6, 5))*(m + 1)/arb(s1))*(2/hp)*(hp/2)**(m - 2)*arb(s1)**(3 - m)/(m - 3)
                Cosc2 += A(WMAX[k]*dm)*(t1 + t2)              # = Lam^2 C_osc2 (ROUTEA's Cosc2 without the 1/Lam^2)
    return up(Csm), up(Cosc2)


def M1_small(beta):
    """[P] lower bound of M1 on (0, beta] for beta <= 0.1 (the bound is decreasing in beta)."""
    nu1 = arb(2).sqrt()/beta
    return lo((arb(1)/24 - 7*beta**2/1920)/(beta/2).cosh() - 2*(-2*PI*nu1).exp()*(nu1/(2*PI) + 1/(4*PI**2)))


CELLS = G.m1_cells()


def M1min_upto(bmax):
    """certified lower bound of min_{0 < beta <= bmax} M1(beta)."""
    v = M1_small(A(Fr(1, 10)))
    for a, b, m in CELLS:
        if A(a) < bmax and m.lower() < v.lower():
            v = m
    return v


def main():
    pr(__doc__.split('Output')[0])
    A1min = 2*W[2].a1 - UMAX*W[3].a1
    pr(f'UMAX = {UMAX} (m >= {float(UMAX)**-0.5:.3f} Lam); A1min (u <= UMAX, both sigma) = {float(A1min):.4f};  eps0 = 1/Lam0^2 = {float(EPS0.mid()):.6e}')
    pr(f'M1 >= {float(M1_small(A(Fr(1, 10))).mid()):.6f} on (0, 0.1] [P];  min M1 on (0, 2 asinh 1] >= '
       f'{float(M1min_upto(2*arb(1).asinh()).mid()):.6f}')
    out = {}
    for s1 in (60, 80, 100, 120):
        PQ = {k: region_I(k, s1) for k in (2, 3)}
        CP = up(2/PI*sum((A(WMAX[k])*PQ[k][0] for k in (2, 3)), arb(0)))
        CQ = up(1/PI*sum((A(WMAX[k])*PQ[k][1] for k in (2, 3)), arb(0)))
        Csm, Co = lam_free_consts(s1)
        pr(f'\ns1 = {s1}: Lam^2 C_P <= {float(CP.mid()):.4f}  [k=2: lead {PQ[2][2]["lead"]:.3f}, eps0 h {PQ[2][2]["epsh"]:.2e}; '
           f'k=3: lead {PQ[3][2]["lead"]:.3f}, eps0 h {PQ[3][2]["epsh"]:.2e}];  Lam^2 C_Q <= {float(CQ.mid()):.4e};  C_sm = '
           f'{float(Csm.mid()):.4f};  Lam^2 C_osc2 = {float(Co.mid()):.4f}   (all valid for every Lam >= Lam0)')
        def U(b):
            b = A(b)
            return up((CP + CQ*(1 - arb(s1)**2/b**2)**arb(-1.5)/b**2 + Csm + Co)/b)
        # boxes from bhat = s1 * 1.03 upward to 2 Lam0; then the half-line bhat >= 2 Lam0 (beta <= 2 asinh 1)
        boxes = []
        b = Fr(s1)*Fr(103, 100)
        while b < 2*Fr(998499874812, 10**9):
            bb = b*Fr(103, 100); boxes.append((b, bb)); b = bb
        L0u = up(G.LAM0)
        tailU = U(Fr(1997))         # U decreasing; for bhat >= 2 Lam0 > 1997 every beta <= 2 asinh 1
        tail_ok = A(A1min)*M1min_upto(2*arb(1).asinh()) > tailU
        b0 = None; rows = []
        for (ba, bb) in reversed(boxes):
            bmax = 2*(A(bb)/(2*lo(G.LAM0))).asinh()
            lhs = A(A1min)*M1min_upto(bmax); rhs = U(ba)
            rows.append((float(ba), float(lhs.mid()), float(rhs.mid())))
            if lhs > rhs:
                b0 = ba
            else:
                break
        pr(f'   tail bhat >= 2 Lam0: A1min min M1 = {float((A(A1min)*M1min_upto(2*arb(1).asinh())).mid()):.4f} > U(1997) = '
           f'{float(tailU.mid()):.4f}: {tail_ok}')
        pr(f'   => certified for EVERY Dg >= 1000: curv > 0 on bhat in [{float(b0) if b0 else None:.2f}, 2 Lam]  (m >= {float(UMAX)**-0.5:.3f} Lam, both helicities)'
           if (b0 and tail_ok) else '   => no window')
        for ba, l, r in sorted(rows)[:4]:
            pr(f'      box start {ba:8.2f}: A1min min M1 = {l:.4f} vs U = {r:.4f}')
        out[s1] = dict(CP=float(CP.mid()), CQ=float(CQ.mid()), Csm=float(Csm.mid()), Cosc2=float(Co.mid()),
                       b0=float(b0) if b0 else None, tail_ok=bool(tail_ok))
    best = min((v['b0'], s) for s, v in out.items() if v['b0'] and v['tail_ok'])
    pr(f'\nbest uniform window start: bhat0 = {best[0]:.2f} (s1 = {best[1]})')
    json.dump(out, open(os.path.join(G.CACHE, 'gl7_routea_uniform_u%s.json' % str(UMAX).replace('/', '_')), 'w'), indent=0)
    open(os.path.join(G.CACHE, 'gl7_routea_uniform_u%s.txt' % str(UMAX).replace('/', '_')), 'w').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
