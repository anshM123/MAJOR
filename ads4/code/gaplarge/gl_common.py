"""
gaplarge/gl_common.py (GAP_LARGESPIN.md §1): Lemma LR-inf bound pieces with Lam as an explicit PARAMETER and the s^4-monotone
envelope R^ of the small-x majorant, so that Lam^4 x (every bound piece) is NON-INCREASING in Lam (GAP_LARGESPIN §1.2, [P]).

Read-only imports: lrtight/lrt_common.py (exact psi_k, R_inf, d1, sumc, P_{-1/2} in Arb), lrtight/lrt2_m1.py (cells(), region-(c)
M1 bound).  Read-only data: lrtight/cache/lrt1_constants.json (certified int|psi_k|, A1min, A2'), lrtight/cache/lrt2_m1.json
(certified M1 lower bounds on 830 beta-cells of [0.1, 24]).  Nothing outside gaplarge/cache is written; no bytecode is written.

Envelope.  |r(s)| <= R(s) = min(R_inf(s), a1/s^2 + I), I = int|psi| (LRINF_TIGHT).  s^4 R_inf(s) = sum_j c_j s^{3-j} (j >= 3) is
non-increasing, a1 s^2 + I s^4 is increasing (a1 > 0 for k = 2, 3), so for ANY s1 > 0
      V := max(s1^4 R_inf(s1), a1 s1^2 + I s1^4)  >=  sup_s s^4 R(s),
and R^(s) := min(R_inf(s), V/s^4) >= R(s) >= |r(s)|, with s^4 R^(s) = min(s^4 R_inf(s), V) non-increasing in s.
"""
import os, sys, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
LRT = os.path.abspath(os.path.join(HERE, '..', 'lrtight'))
sys.path.insert(0, LRT)
from fractions import Fraction as Fr
from flint import arb, ctx
import lrt_common as LC                      # read-only
from lrt_common import Weight, A, up, lo, ball, P_mhalf
from lrt2_m1 import M1_lower_region_c, cells  # read-only (definitions only; its main() is not called)

ctx.prec = 128
CACHE = os.path.join(HERE, 'cache')
LRT_CACHE = os.path.join(LRT, 'cache')

CONST = json.load(open(os.path.join(LRT_CACHE, 'lrt1_constants.json')))
A1MIN = Fr(CONST['A1min']); A2P = Fr(CONST['A2p'])
C3 = arb(3)/8*arb(2)**(arb(5)/2)
SQ2PI = arb(2).sqrt()/arb.pi()


def lam_of(Dg):
    """Lam = sqrt((Dg-2)(Dg-1)) as an Arb ball (Dg an integer or Fraction)."""
    return A(Fr(Dg - 2)*Fr(Dg - 1)).sqrt()


LAM0 = lam_of(1000)                            # = lrt_common.LAM
assert abs(float((LAM0 - LC.LAM).mid())) < 1e-30


class HatWeight(Weight):
    """Weight with the envelope R^ (see module doc)."""

    def __init__(self, k):
        super().__init__(k, absint_up=Fr(CONST[str(k)]['absint_upper_rational']))
        a1 = float(self.a1); I = float(self.absint_up)
        f = lambda s: s**4*float(self.Rinf(arb(s)).mid()) - (a1*s*s + I*s**4)
        lo_, hi_ = 0.1, 1000.0
        assert f(lo_) > 0 and f(hi_) < 0
        for _ in range(200):
            mid = (lo_*hi_)**0.5
            if f(mid) > 0:
                lo_ = mid
            else:
                hi_ = mid
        s1 = A(Fr(lo_).limit_denominator(2**40))
        v1 = s1**4*self.Rinf(s1); v2 = A(self.a1)*s1**2 + A(self.absint_up)*s1**4
        self.s1 = s1
        self.V = up(v1) if v1.upper() >= v2.upper() else up(v2)

    def Rhat(self, x):
        b1 = self.Rinf(x)
        b2 = self.V/x**4
        return b1 if b1.upper() <= b2.upper() else b2


W = {k: HatWeight(k) for k in (2, 3)}
WK = ((2, 2), (3, 1))                          # weights 2|E2| + |E3|   (u <= 1)

T1 = 1e-6; QG = 1.0005; TMAX = 32.5


def build_I(w, LAM, hat=True):
    """cumulative UPPER sums of I1(T) = int_0^T Rx(Lam t) t^4 cosh t dt and I2(T) = int_0^T Rx(Lam t)(cosh t - 1)^2 dt on the
    geometric t-grid of lrt3 (Rx decreasing in t: left end; t^4 cosh t, (cosh t - 1)^2 increasing: right end), Rx = R^ (hat=True)
    or LRINF_TIGHT's R (hat=False).  [0, t1]: R^(Lam t) <= V/(Lam t)^4 -> V t1 cosh t1/Lam^4 and V t1 cosh^2 t1/(4 Lam^4);
    R(Lam t) <= a1/(Lam t)^2 + I -> lrt3's closed forms."""
    ts = [T1]
    while ts[-1] < TMAX:
        ts.append(ts[-1]*QG)
    ts = [arb(t) for t in ts]
    t1 = ts[0]; CC = t1.cosh()
    if hat:
        I1 = w.V*t1*CC/LAM**4
        I2 = w.V*t1*CC**2/(4*LAM**4)
        Rx = w.Rhat
    else:
        a1 = A(w.a1); ai = A(w.absint_up)
        I1 = a1*t1**3/(3*LAM**2)*CC + ai*t1**5/5*CC
        I2 = (a1*t1**3/(3*LAM**2) + ai*t1**5/5)*CC**2/4
        Rx = w.R
    cum1 = [up(I1)]; cum2 = [up(I2)]
    for i in range(len(ts) - 1):
        a, b = ts[i], ts[i + 1]
        Ra = Rx(LAM*a)
        h = b - a
        I1 = I1 + up(Ra*b**4*b.cosh()*h)
        I2 = I2 + up(Ra*(b.cosh() - 1)**2*h)
        cum1.append(up(I1)); cum2.append(up(I2))
    return ts, cum1, cum2


def lookup(ts, cum, T):
    lo_, hi_ = 0, len(ts) - 1
    assert ts[-1] >= T
    while lo_ < hi_:
        mid = (lo_ + hi_)//2
        if ts[mid] >= T:
            hi_ = mid
        else:
            lo_ = mid + 1
    return cum[lo_]


def E_up(k, ba, bb, TAB, LAM):
    """certified upper bound of |E_k| on the cell [ba, bb] at Lam (LRINF_TIGHT §1.4 with R -> R^)."""
    w = W[k]; ts, c1, c2 = TAB[k]
    Ba = A(ba); Bb = A(bb)
    D = Ba.cosh() - 1; K0 = 1/D.sqrt(); X = LAM*Ba
    Ta = K0/LAM*w.d1f(X)
    Tb = K0/(4*D*LAM**3)*w.sumc(X)
    Tc1 = K0*(lookup(ts, c1, Bb/2)/(48*D) + C3*lookup(ts, c2, Bb/2)/D**2)
    Tc2 = w.Rinf(LAM*Ba/2)*(arb.pi()/arb(2).sqrt())*P_mhalf(Ba)
    E = SQ2PI*(Ta + Tb + Tc1 + Tc2)
    return up(E), dict(Ta=SQ2PI*Ta, Tb=SQ2PI*Tb, Tc1=SQ2PI*Tc1, Tc2=SQ2PI*Tc2)


def tail_majorant(TAB, LAM, BT=64):
    """for every beta >= BT: Lam^2 (2|E2| + |E3|)/(A1min M1) <= returned value (LRINF_TIGHT §1.5; TAB built with R^ or R)."""
    B = arb(BT); T0 = arb(32)
    ln2 = arb(2).log()
    ell = 1 - B/(2*(B/2).sinh()) - arb(2).sqrt()*(B*(-B/2).exp()/A(Fr(47, 10)) + 4/B) - 2*(2*ln2).sqrt()*B/(B - ln2)**2
    assert ell > 0
    D = B.cosh() - 1
    C1 = 2*arb(2).sqrt()*(arb(1)/4).exp(); C2 = (2/(1 - (-arb(1)).exp())).sqrt()
    tot = arb(0)
    for k, wk in WK:
        w = W[k]; ts, c1, c2 = TAB[k]
        I1T0 = lookup(ts, c1, T0); I2T0 = lookup(ts, c2, T0)
        rinf0 = w.Rinf(LAM*T0)                 # R^, R <= R_inf; t^4 R_inf(Lam t) decreasing => <= rho0 on t >= 32
        rho0 = rinf0*T0**4
        f_a = B*w.d1f(LAM*B)/LAM
        f_b = B*w.sumc(LAM*B)/(4*D*LAM**3)
        f_c1 = B*I1T0/(48*D) + rho0*B/(96*(B/2).sinh())
        f_c2 = C3*(B*I2T0/D**2 + rinf0*(B*B/(16*(B/2).sinh()**4) + B*(B/2).cosh()/(8*(B/2).sinh()**3)))
        f_d = B*w.Rinf(LAM*B/2)*(C1 + C2*B)/arb(2).sqrt()
        tot += wk*(f_a + f_b + f_c1 + f_c2 + f_d)
    return up(LAM**2/(A(A1MIN)*ell)*tot), ell


def m1_cells():
    """certified M1 lower bounds: lrt2 cells on [0.1, 24] (floats stored by lrt2; shaved by a relative 2^-40 to absorb the
    float rounding of the stored lower endpoints) and region-(c) cells of width 0.5 on (24, 64]."""
    M1c = json.load(open(os.path.join(LRT_CACHE, 'lrt2_m1.json')))
    cl = [(Fr(a), Fr(b), arb(m)*(1 - arb(2)**-40)) for a, b, m in M1c]
    for (a, b) in cells([('24', '64', '0.5')]):
        cl.append((a, b, lo(M1_lower_region_c(ball(a, b)))))
    return cl
