"""
Uniformity of the heavy expansion of the exact shadow-representation integrand T(z_1,z_3;x,y) (ads4/code/e2e/vertex.py) in the
directions of x and y.  The moment is  2 int d^3x d^3y (|x||y|)^{k-4} K(xhat,yhat) T(x,y)  (identity (*) of SPINNING_E2E.md sect. 3.1).
The (r, eta) corner u -> 0 corresponds to a large hyperbolic distance d between xhat and yhat.  If the coefficient c_2 of the
pointwise 1/m^2 correction grew like e^{p d}, term-by-term integration would fail at the order where the growth beats the decay
of the leading integrand (~ e^{-(a_x - 1) d}).  We measure m^2 max|T_sym - T_lead|/max|T_lead| (over 5x5 polarisation pairs) at fixed X = m x, Y = m y with
|X| = 1/cosh d(xhat, phat), |Y| = 1/cosh d(yhat, pbarhat) (the typical moduli at fixed direction) for growing d, in the Weyl-odd
leading form, with the polarisation-symmetrised tensor (sect. 6.1).  Output: out_u1_uniformity.txt
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import mpmath as mp
from vertex import integrand_exact, integrand_leading_weyl
mp.mp.dps = 50
out = []
def pr(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
def hyp(a, phi):  # unit future timelike vector, rapidity a, azimuth phi
    return [mp.cosh(a), mp.sinh(a)*mp.cos(phi), mp.sinh(a)*mp.sin(phi)]
def mdot(u, v): return -u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
def nullv(th): return [mp.mpf(1), mp.cos(th), mp.sin(th)]
def boost_frame(yh):
    """null polarisations attached to yhat: z(th) = yhat + cos th e1 + sin th e2 with e1, e2 orthonormal spacelike, orthogonal to yhat."""
    # e1: in the plane of yhat and the time axis
    t = [mp.mpf(1), mp.mpf(0), mp.mpf(0)]
    c = -mdot(t, yh)
    e1 = [t[i] - c*yh[i] for i in range(3)]; n1 = mp.sqrt(mdot(e1, e1)) if mdot(e1, e1) > 0 else None
    if n1 is None or n1 < mp.mpf('1e-30'):
        e1 = [mp.mpf(0), mp.mpf(1), mp.mpf(0)]; e2 = [mp.mpf(0), mp.mpf(0), mp.mpf(1)]
    else:
        e1 = [v/n1 for v in e1]
        # e2 = yhat x e1 (Lorentzian cross product, up to sign)
        a, b = yh, e1
        e2 = [-(a[1]*b[2] - a[2]*b[1]), -(a[2]*b[0] - a[0]*b[2]), (a[0]*b[1] - a[1]*b[0])]
        e2 = [e2[0], e2[1], e2[2]]
        n2 = mp.sqrt(abs(mdot(e2, e2))); e2 = [v/n2 for v in e2]
    return lambda th: [yh[i] + mp.cos(th)*e1[i] + mp.sin(th)*e2[i] for i in range(3)]
if __name__ == "__main__":
    pr("# m^2 max_pol|T_sym - T_lead|/max_pol|T_lead| (5x5 polarisation pairs on the circle attached to yhat) at fixed X = m x, Y = m y, |X| = 1/cosh d(xhat,phat), |Y| = 1/cosh d(yhat,pbarhat); Delta/J fixed;")
    pr("# columns: m = 4.6e2, 4.6e3, 4.6e4 (Delta = 500.25, 5000.25, 50000.25, J = 0.4 Delta); polarisations z(th) attached to yhat.")
    pr("# Uniformity <=> the converged value does not grow exponentially with the hyperbolic distances.")
    for l in (1, 2):
        for (ax, phx, ay, phy) in ((0.0, 0.0, 0.0, 0.0), (2.0, 1.0, 0.0, 0.0), (4.0, 1.0, 0.0, 0.0), (6.0, 1.0, 0.0, 0.0), (8.0, 1.0, 0.0, 0.0),
                                   (4.0, 2.5, 0.0, 0.0), (8.0, 2.5, 0.0, 0.0), (0.0, 0.0, 4.0, 2.0), (0.0, 0.0, 8.0, 2.0), (6.0, 1.0, 6.0, 4.0)):
            row = []
            for (D, Jv) in ((500.25, 200), (5000.25, 2000), (50000.25, 20000)):
                D = mp.mpf(D)
                Ap, Am = D - mp.mpf(3)/2, Jv + mp.mpf(1)/2
                m = mp.sqrt(Ap*Ap - Am*Am)
                beta2 = mp.atanh(Am/Ap)                      # phat at rapidity +beta/2 along axis 1, pbarhat at -beta/2
                ph = hyp(beta2, 0); pbh = hyp(-beta2, 0)
                xh = hyp(ax, phx); yh = hyp(ay, phy)
                # move xhat, yhat to be measured from phat, pbarhat: use Lorentz boosts along axis 1 (exact, simple)
                def boost(v, a):
                    return [mp.cosh(a)*v[0] + mp.sinh(a)*v[1], mp.sinh(a)*v[0] + mp.cosh(a)*v[1], v[2]]
                xh = boost(xh, beta2); yh = boost(yh, -beta2)
                dx = mp.acosh(max(mp.mpf(1), mp.re(-mdot(xh, ph)))); dy = mp.acosh(max(mp.mpf(1), mp.re(-mdot(yh, pbh))))
                X = [v/mp.cosh(dx) for v in xh]; Y = [v/mp.cosh(dy) for v in yh]
                x = [v/m for v in X]; y = [v/m for v in Y]
                zf = boost_frame(yh)
                ths = [mp.mpf(0.3) + 2*mp.pi*j/5 for j in range(5)]
                num = mp.mpf(0); den = mp.mpf(0)
                for i1 in range(5):
                    for i3 in range(i1, 5):
                        z1, z3 = zf(ths[i1]), zf(ths[i3])
                        le, ve = integrand_exact(l, D, Jv, 1.75, x, y, z1, z3); ll, vl = integrand_leading_weyl(l, D, Jv, 1.75, x, y, z1, z3)
                        le2, ve2 = integrand_exact(l, D, Jv, 1.75, x, y, z3, z1); ll2, vl2 = integrand_leading_weyl(l, D, Jv, 1.75, x, y, z3, z1)
                        ex = mp.exp(le - ll)*ve + mp.exp(le2 - ll2)*ve2; ld = vl + vl2
                        num = max(num, abs(ex - ld)); den = max(den, abs(ld))
                rs = 1 + num/den
                row.append(mp.nstr(m*m*(rs - 1), 6))
                dxy = mp.acosh(max(mp.mpf(1), -mdot(xh, yh)))
            pr(f"  l={l}  d(xhat,phat)={float(dx):4.1f}  d(yhat,pbarhat)={float(dy):4.1f}  d(xhat,yhat)={float(dxy):5.2f}:  m^2 * rel.dev = {row}")
    open(os.path.join(HERE, 'out_u1_uniformity.txt'), 'w').write("\n".join(out) + "\n")
