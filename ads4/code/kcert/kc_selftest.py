"""kc_selftest.py -- consistency checks of the ball evaluators used by kc_dk0.py (not part of the certificate's logic).
 (1) chi_re (windowed, with tail ball) vs a full, non-windowed sum and vs the published convenience table (double precision);
 (2) chi/z near 0 (series branch) vs chi(z)/z away from 0;
 (3) acb.integral with the windowed evaluator on a test integrand with a closed form: int_0^b Q(nu) g(nu - c) dnu is not elementary,
     so instead int_{-oo}^{oo} exp(-(nu-c)^2/(2 s^2)) dnu = s sqrt(2 pi) is checked through Chi's machinery with a one-centre weight;
 (4) dK0 vs the published double-precision values (CERT_CLOSE Sec. 3.2).
"""
import sys, os
import numpy as np
from flint import arb, acb, ctx
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kc_common import *  # noqa

ctx.prec = 128
ok = True
for Dg in (300, 1000):
    C = Chi(Dg)
    tab = np.loadtxt(os.path.join(DATA, 'k2', 'k2_continuum_D_%d_table.txt' % Dg))
    idx = np.linspace(0, len(tab) - 1, 41).astype(int)
    worst_tab = 0.0; worst_full = 0.0
    for i in idx:
        nu = tab[i, 0]; z = acb(nu)
        v = C(z)
        # full sum, no window
        full = acb(0)
        for c, s, x in zip(C.c, C.s, C.x):
            g = (-(z - c) ** 2 / (2 * arb(s) ** 2)).exp()
            if c > 0:
                g += (-(z + c) ** 2 / (2 * arb(s) ** 2)).exp()
            full += arb(x) * g
        full = C.Q(z) * full + C.sing(z, False)
        worst_full = max(worst_full, abs(float((v - full).real.mid())), float(v.real.rad()))
        worst_tab = max(worst_tab, abs(float(v.real.mid()) - tab[i, 1]) / max(1e-300, np.abs(tab[:, 1]).max()))
        if not (v.real.overlaps(full.real)):
            ok = False
    print('Dg=%d (1): |windowed - full| <= %.2e (ball radii incl.); max |chi - table|/max|chi| = %.2e' % (Dg, worst_full, worst_tab))
    ok &= worst_tab < 1e-11
    for nu in (0.3, 0.45, 0.6):
        a = C.over_z(acb(nu)); b = C(acb(nu)) / nu
        print('Dg=%d (2): chi/z at %.2f: series %s vs direct %s' % (Dg, nu, a.real.str(15), b.real.str(15)))
        ok &= a.real.overlaps(b.real)
# (3) one-Gaussian test
class One(Chi):
    def __init__(self):
        self.c = np.array([50.0]); self.s = np.array([0.9]); self.x = np.array([1.0]); self.y = [0.0, 0.0, 0.0]; self.xmax = 1.0
        self._prec = None
o = One()
f = lambda z, an: o.B(z) if z.rad() < 4 else NAN
I1 = sum((acb.integral(f, arb(a), arb(a + 5)) for a in range(30, 70, 5)), acb(0))
exact = arb('0.9') * (2 * arb.pi()).sqrt()
print('(3) int g = %s vs %s (tails beyond [30,70] < 1e-100)' % (I1.real.str(20, radius=True), exact.str(20)))
ok &= abs(float(I1.real.mid()) - float(exact.mid())) < 1e-15
# (4) published values
pub = {300: -7.666324482e-4, 1000: -6.830149984e-5}
import kc_dk0
for Dg in (300, 1000):
    _, r = kc_dk0.main(Dg, quick=True)
    rel = abs(float(r['dK0'].mid()) - pub[Dg]) / abs(pub[Dg])
    print('(4) Dg=%d dK0 ball %s vs published double %.9e (rel %.1e)' % (Dg, r['dK0'].str(12, radius=True), pub[Dg], rel))
    ok &= rel < 1e-8
print('SELFTEST', 'PASS' if ok else 'FAIL')
