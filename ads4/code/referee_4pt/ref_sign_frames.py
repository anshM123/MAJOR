"""Referee M2 (round 2): with the ACTUAL polarisation vectors and positions fed to the blocks (e2e/moments.frame_positions,
fourpt/moments4.pol_grid4), evaluate D5 = det[z_in, z_out, n], n = (x_in - x_out)/|.|, at both pairs, and D1 = z_in.z_out + (z_in.n)(z_out.n),
and compare with sin/cos(th_out - th_in) in the mode parametrisation.  Any eps convention gives c*det, so the product's sign is convention-free."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fourpt')); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
from moments import frame_positions
from moments4 import pol_grid4
eta = np.diag([-1., 1., 1.])
dot = lambda a, b: a @ eta @ b
rng = np.random.default_rng(5)
worst = {}
for trial in range(6):
    ru = np.array([rng.uniform(0.05, 0.9)]); xi = np.array([rng.uniform(0.05, 2.5)])
    x1, x2, x3, x4 = [v[0] for v in frame_positions(ru, xi)]
    # generic angles instead of the 3-point grid
    l = 1
    th = rng.uniform(0, 2*np.pi, 4)       # th1, th3, th2, th4
    ch, sh = np.cosh(xi[0]), np.sinh(xi[0])
    zy = lambda t: np.array([1, np.cos(t), np.sin(t)])
    zx = lambda t: np.array([-ch - np.cos(t)*sh, -sh - np.cos(t)*ch, -np.sin(t)])
    # cross-check the generic formulas against pol_grid4 on its own grid
    Z1, Z3, Z2, Z4 = pol_grid4(xi, 1)
    g = 2*np.pi*np.arange(3)/3
    assert np.allclose(Z2[0, 0, 0, 1, 0], zx(g[1])) and np.allclose(Z1[0, 2, 0, 0, 0], zy(g[2]))
    z1, z3, z2, z4 = zy(th[0]), zy(th[1]), zx(th[2]), zx(th[3])
    # sanity: polarisations null and transverse? (null only)
    assert abs(dot(z1, z1)) < 1e-12 and abs(dot(z2, z2)) < 1e-12
    def pair(zi, zo, xi_, xo_):
        n = xi_ - xo_; n = n/np.sqrt(abs(dot(n, n)))
        return np.linalg.det(np.array([zi, zo, n])), dot(zi, zo) + dot(zi, n)*dot(zo, n), dot(n, n)
    d13, e13, nn13 = pair(z1, z3, x1, x3); d24, e24, nn24 = pair(z2, z4, x2, x4)
    s13, s24 = np.sin(th[1]-th[0]), np.sin(th[3]-th[2])
    c13, c24 = np.cos(th[1]-th[0]), np.cos(th[3]-th[2])
    print(f"trial {trial}: n.n = {nn13:+.2f},{nn24:+.2f};  D5_13/sin = {d13/s13:+.6f}, D5_24/sin = {d24/s24:+.6f};  "
          f"product ratio {d13*d24/(s13*s24):+.6f};  D1_13/cos = {e13/c13:+.6f}, D1_24/cos = {e24/c24:+.6f}")
