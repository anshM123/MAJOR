"""
lemU_cert/lu_bounds.py [P + Arb]: rigorous envelopes of the per-frequency jets for complex frequency.

Lemma J (jet envelope).  For real c in (0,1), d = 1 - c, p in {3,4}, 0 < r < 1/c, and complex nu = nu_r + i nu_i:
    |y_j(nu)| = |d^p [t^j] F_{p,nu}(c + c d t)|  <=  r^{-j} (d / min|1-z|)^p exp((1/2 + |nu_i|) Lw(r) + |nu_r| Th(r)),
where on the circle z = c + c d r e^{i alpha} and for phi in [0, pi], with w = (1 + z + 2 sqrt(z) cos phi)/(1 - z):
    Lw(r) >= max |log|w||,   Th(r) >= max |arg w|.
Proof.  Cauchy in t on |t| = r.  |z| <= c + c d r < 1, so Re X > 0 for X = (1+z)/(1-z), and Laplace's first integral
P_mu(X) = (1/pi) int_0^pi (X + sqrt(X^2-1) cos phi)^mu dphi holds (Whittaker-Watson 15.23; X + sqrt(X^2-1) cos phi = w).
Both factors (1 + sqrt(z) e^{+-i phi}) have positive real part (|sqrt z| < 1); the enclosure below checks |arg w| < pi on every box,
so the principal power is analytic in z on the disc and equals P_mu on the real segment (identity theorem).  With
mu = -1/2 + i nu:  |w^mu| = |w|^{-1/2 - nu_i} e^{-nu_r arg w} <= exp((1/2 + |nu_i|)|log|w|| + |nu_r||arg w|).  []
Enclosures: Arb boxes over alpha in [0, pi] (conjugate symmetry z -> conj z covers [pi, 2pi]); phi is handled exactly, since
w(phi) = A + B cos(phi) runs over a segment (see _box_eval).
"""
import math
from flint import arb, acb, ctx

PI = None


def _box_eval(c, d, r, a0, a1):
    """enclosures on the alpha-box [a0, a1]*pi of the circle.  With q = sqrt(z):  w(phi) = (1 + q^2 + 2 q cos phi)/(1 - q^2)
    runs over the segment between w_- = (1-q)/(1+q) = e^{-xi} and w_+ = (1+q)/(1-q) = e^{xi}.  If 0 is not on the segment,
    arg w is monotone along it, so max|arg w| = max(|arg w_+|, |arg w_-|) = |Im xi|; and min |w| = min_s |(1-q)^2 + 2q(1+s)|/|1-q^2|
    over s in [-1,1], attained at clamp(s*), s* = -Re((1+q^2) conj q)/(2|q|^2) (convexity)."""
    al = arb((a0 + a1)/2, (a1 - a0)/2)*PI
    t = acb(r*al.cos(), r*al.sin())
    z = c + c*d*t
    q = z.sqrt()
    one_z = 1 - z
    iv = 1/abs(one_z)
    wp = (1 + q)/(1 - q); wm = (1 - q)/(1 + q)
    ap, am = abs(wp), abs(wm)
    if not (ap.lower() > 0 and am.lower() > 0) or not iv.is_finite():
        return None
    th = abs(wp.arg())
    th2 = abs(wm.arg())
    if th2.upper() > th.upper():
        th = th2
    if not (th.upper() < PI*0.999):
        return None
    sstar = -((1 + q*q)*q.conjugate()).real/(2*abs(q)**2)
    lo = max(-1.0, min(1.0, float(sstar.lower()))); hi = max(-1.0, min(1.0, float(sstar.upper())))
    S1 = arb((lo + hi)/2 + 1, (hi - lo)/2 + 1e-30)            # 1 + s, s in [lo, hi]
    num = (1 - q)**2 + 2*q*S1
    dist = arb((abs(num)*iv).lower()) if False else arb((abs(num)/abs(1 - q*q)).lower())
    if not dist > 0:
        return None
    lw = abs(ap.log())
    for v in (abs(am.log()), abs(dist.log())):
        if v.upper() > lw.upper():
            lw = v
    if not lw.is_finite():
        return None
    return lw, th, iv


def _box_stats(c, r, na=400, nphi=None):
    """upper bounds (Lw, Th, inv), inv >= max 1/|1-z|, on the circle of radius r (in t); adaptive splitting in alpha."""
    global PI
    PI = arb.pi()
    d = 1 - c
    Lw = arb(0); Th = arb(0); inv = arb(0)
    stack = [(i/na, (i + 1)/na) for i in range(na)]
    while stack:
        a0, a1 = stack.pop()
        res = _box_eval(c, d, r, a0, a1)
        if res is None:
            if a1 - a0 < 1e-9:
                raise ValueError('box splitting failed (r too close to 1/c?)')
            am = (a0 + a1)/2
            stack += [(a0, am), (am, a1)]
            continue
        lw, th, iv = res
        if lw.upper() > Lw.upper():
            Lw = arb(lw.upper())
        if th.upper() > Th.upper():
            Th = arb(th.upper())
        if iv.upper() > inv.upper():
            inv = arb(iv.upper())
    return Lw, Th, inv


class Envelope:
    """precomputed (r_k, Lw_k, Th_k, inv_k) on a radius grid; log-envelope of the jets for frequency sets."""

    def __init__(self, c, nr=60, rmin=0.01, frac=0.995, na=400, nphi=None):
        old = ctx.prec
        ctx.prec = 128
        try:
            self.c = arb(c.mid()) if c.rad() == 0 else c   # c must be exact here (it is: a dyadic)
            assert c.rad() == 0
            rmax = frac/float(c.mid())
            self.rows = []
            for k in range(nr):
                r = rmin*(rmax/rmin)**(k/(nr - 1))
                r = arb(r)
                Lw, Th, inv = _box_stats(self.c, r, na)
                self.rows.append((r, r.log(), Lw, Th, (1 - self.c)*inv))
        finally:
            ctx.prec = old

    def log_ybar(self, J, p, nur, nui):
        """list of arb upper bounds of log max|y_j| (j = 0..J) over complex nu with |nu_r| <= nur, |nu_i| <= nui."""
        old = ctx.prec
        ctx.prec = 64
        try:
            nur = arb(nur); nui = arb(nui)
            base = []
            for (r, lr, Lw, Th, dinv) in self.rows:
                base.append((lr, p*dinv.log() + (arb(0.5) + nui)*Lw + nur*Th))
            out = []
            for j in range(J + 1):
                best = None
                for lr, b0 in base:
                    v = b0 - j*lr
                    if best is None or v.upper() < best.upper():
                        best = v
                out.append(arb(best.upper()))
            return out
        finally:
            ctx.prec = old
