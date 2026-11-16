"""
Referee probe of the 'outer region' of identity (*) (SPINNING_LOWK.md sect. 6): the exact shadow integrand along rays of growing
zeta = -Phat.X, from the typical region (zeta ~ 1) up to zeta ~ m (Cartesian components of x of order 1), compared with the
leading heavy integrand.  Prints log10 of  |x|^{k-1} max_pol |T_exact|  (the (*) density per d log|x|, k = 2) relative to its value at zeta = 1,
and log10(|T_exact|/|T_lead|).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import mpmath as mp
from vertex import integrand_exact, integrand_leading_weyl
mp.mp.dps = 40
def hyp(a, phi): return [mp.cosh(a), mp.sinh(a)*mp.cos(phi), mp.sinh(a)*mp.sin(phi)]
def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def boost(v, a): return [mp.cosh(a)*v[0] + mp.sinh(a)*v[1], mp.sinh(a)*v[0] + mp.cosh(a)*v[1], v[2]]
l = int(sys.argv[1]) if len(sys.argv) > 1 else 1
k = 2
for (D, Jv) in ((mp.mpf('50.25'), 20), (mp.mpf('200.25'), 80)):
    Ap, Am = D - mp.mpf(3)/2, Jv + mp.mpf(1)/2
    m = mp.sqrt(Ap*Ap - Am*Am); b2 = mp.atanh(Am/Ap)
    ph = hyp(b2, 0); pbh = hyp(-b2, 0)
    print(f"# l={l} Delta={D} J={Jv} m={float(m):.1f}")
    for (ax, phx, ay, phy, which) in ((0, 0, 0, 0, 'x'), (3, 1.0, 0, 0, 'x'), (0, 0, 0, 0, 'y'), (0, 0, 3, 2.0, 'y'), (2, 2.5, 2, 1.0, 'x')):
        xh = boost(hyp(ax, phx), b2); yh = boost(hyp(ay, phy), -b2)
        dx = mp.acosh(max(1, -mdot(xh, ph))); dy = mp.acosh(max(1, -mdot(yh, pbh)))
        pols = [[mp.mpf(1), mp.cos(t), mp.sin(t)] for t in (0.3, 1.9, 3.5, 5.1)]
        ref = None; row = []
        dex = mp.acosh(max(1, xh[0])); dey = mp.acosh(max(1, yh[0]))
        zxmax = m*mp.exp(-dex)*mp.cosh(dx); zymax = m*mp.exp(-dey)*mp.cosh(dy)
        zmax = zxmax if which == 'x' else zymax
        print(f"   support: zeta_max = {float(zmax):.2f}")
        for zeta in (mp.mpf(1), mp.mpf(4), mp.mpf(16)) + tuple(f*zmax for f in (0.25, 0.5, 0.7, 0.85, 0.95, 0.99)):
            zx = zeta if which == 'x' else mp.mpf(1); zy = zeta if which == 'y' else mp.mpf(1)
            x = [v/mp.cosh(dx)*zx/m for v in xh]; y = [v/mp.cosh(dy)*zy/m for v in yh]
            best = None; bestl = None
            for z1 in pols:
                for z3 in pols:
                    try:
                        le, ve = integrand_exact(l, D, Jv, 1.75, x, y, z1, z3)
                        ll, vl = integrand_leading_weyl(l, D, Jv, 1.75, x, y, z1, z3)
                    except Exception as e:
                        le, ve, ll, vl = mp.mpf('-inf'), 0, mp.mpf(0), 1
                    v = le + mp.log(abs(ve)) if ve != 0 else mp.mpf('-inf')
                    w = ll + mp.log(abs(vl)) if vl != 0 else mp.mpf('-inf')
                    best = v if best is None or mp.re(v) > mp.re(best) else best
                    bestl = w if bestl is None or w > bestl else bestl
            r = mp.sqrt(-mdot(x, x)) if which == 'x' else mp.sqrt(-mdot(y, y))
            dens = mp.re(best) + (k - 1)*mp.log(r)
            if ref is None: ref = dens
            row.append((float(zeta), dens, float((mp.re(best) - bestl)/mp.log(10))))
        print(f"  ray in {which}: d(x,p)={float(dx):.1f} d(y,pb)={float(dy):.1f}: " +
              "  ".join(f"z={z:.3g}: dens {float((d - ref)/mp.log(10)):+.1f}, ex/lead {e:+.2f}" for z, d, e in row), flush=True)
