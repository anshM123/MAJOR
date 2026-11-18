"""
g3 (M2c): x-independence of the relative crystal normalisation of the CLS-D (2,4) building blocks, which makes Q^2 and O*Q images of FIXED
covariant graviton structures.
CLS-D (sec. 4.3, App. D): j'=0 basis for conserved O2, O4: (-2V_{0,42})^{J_F - J2 - J4 + 2n} H20^{J2-n} H40^{J4-n} H24^n / (common X powers).
Photons (J2 = J4 = 1): s1 = (-2V)^{J_F-2} H20 H40, s2 = (-2V)^{J_F} H24.  Crystal images (eqs. H20H40_kernel_toflat, H24_kernel_toflat at D = 4):
  H20 H40 = X20 X40 * PP/(4 nu^2 (m^2-nu^2)),   H24 = -X24 * (VV + WW)/(nu^2 (m^2-nu^2)^2).
So s2/s1 |crystal = rho(x) * (VV + WW)/VV with rho = -(2V)^2 X24/(X20 X40) * 4/(m^2 - nu^2).  CLS-D's printed photon kernels (fixed dual
combinations of s1, s2) have x-independent coefficients in (VV, WW) only if rho is x-independent; then the graviton structures s_n, n = 0,1,2,
map to (common) * P^{2-n} (rho H)^n with the same constant rho, and Q^2 = (H - 2P)^2 is the image of a fixed combination.
Here rho(x) is evaluated from the scattering-crystal coordinates (CLS-D eq. scattering_crystal, d = 3, lightcone x = (u, v, y),
x^2 = -uv + y^2), with V_{0,42} = [(z0.x40) x20^2 - (z0.x20) x40^2]/x24^2 (CLS-D eq. bulkpoint_integrand_factor up to normalisation).
Usage: python3 g3_homogeneity.py -> cache/g3_homogeneity.txt
"""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
def sq(w): return -w[0]*w[1] + w[2]**2
def dot(a, b): return -(a[0]*b[1] + a[1]*b[0])/2 + a[2]*b[2]
out = []
for xx in (0.05, 0.2, 0.4, 0.6, 0.8, 0.95):
    m, nu = 1.0, xx
    a = (m - np.sqrt(m*m - nu*nu))/nu
    X = np.array([-1j*a, 1j*a, 0.0])
    x4 = X/(-sq(X)); x2 = -x4
    vals = []
    for sgn in (+1, -1):                                   # the two points of the S^0 locus
        x0 = np.array([0.0, 0.0, sgn*a]); z0 = np.array([1.0, 0.0, 0.0])
        x40, x20, x24 = x4 - x0, x2 - x0, x2 - x4
        V = (dot(z0, x40)*sq(x20) - dot(z0, x20)*sq(x40))/sq(x24)
        rho = -(2*V)**2*sq(x24)/(sq(x20)*sq(x40))*4/(m*m - nu*nu)
        vals.append(rho)
    out.append(f"x = {xx}: rho = {vals[0]:.12f} (point +), {vals[1]:.12f} (point -)")
    print(out[-1])
open(os.path.join(HERE, 'cache', 'g3_homogeneity.txt'), 'w').write("\n".join(out) + "\n")
