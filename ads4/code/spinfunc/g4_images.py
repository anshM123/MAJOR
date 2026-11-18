"""
g4 (M2c): independent evaluation of the crystal images of the CLS-D photon dual-structure building blocks H20 H40/(X20 X40) and H24/X24
(D = 4, d = 3), using CLS-D's definitions (eq. CFT_building_blocks_VH: X_ij = -2 X_i.X_j = x_ij^2, V_{i,jk} = (Z_i.X_j X_i.X_k - Z_i.X_k X_i.X_j)/
(X_j.X_k), H_ij = -2(Z_i.Z_j X_i.X_j - Z_i.X_j Z_j.X_i); embedding X = (1, x^2, x), Z = (0, 2x.z, z)), the scattering crystal (eq.
scattering_crystal: x4 = -x2 = x/(-x^2), x = (-i a, i a, 0), a = (m - sqrt(m^2-nu^2))/nu, x0 = (0, 0, +-a), z0 = (1, 0, 0) lightcone) and the
polarisation maps (eqs. polarvec_125_CFTtoflat, polarvec_346_CFTtoflat: z2 = (t=0, x = -e2^n, y = e2_perp), z4 = (0, -e4^{n'}, e4_perp)), with
flat e2 = e2^n n_perp + e2_perp x3-hat, n_perp = (0, -nu/m, s/m, 0); e4 = e4^{n'} n'_perp + e4_perp x3-hat, n'_perp = (0, nu/m, s/m, 0).
Compared with CLS-D's printed images (eqs. H20H40_kernel_toflat, H24_kernel_toflat at D = 4):
   H20H40/(X20X40) = (e2.V2)(e4.V4)/(4 nu^2 (m^2-nu^2)),      H24/X24 = -[(e2.V2)(e4.V4) + (e2.W2)(e4.W4)]/(nu^2 (m^2-nu^2)^2),
and with the combination implied by CLS-D's printed photon kernel K^(1)_{j'=0} (eq. JJO_kernels_flatspace) through the printed dual
coefficients (eq. JJO_kernel_dual).  Result printed as 2x2 bilinear forms on the basis e2 in {n_perp, x3}, e4 in {n'_perp, x3}.
Usage: python3 g4_images.py -> cache/g4_images.txt
"""
import os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))

def ldot(a, b):                      # lightcone (u, v, y): a.b = -(a_u b_v + a_v b_u)/2 + a_y b_y
    return -(a[0]*b[1] + a[1]*b[0])/2 + a[2]*b[2]
def cart_to_lc(t, x, y):            # u = t - x, v = t + x  (so -uv = -t^2 + x^2)
    return np.array([t - x, t + x, y], dtype=complex)
def emb(x):  return ('X', x)
def EX(x):   return (1.0 + 0j, ldot(x, x), x)          # embedding X = (X+, X-, x)
def EZ(x, z): return (0.0 + 0j, 2*ldot(x, z), z)
def edot(A, B):                      # X.Y = -(A+ B- + A- B+)/2 + a.b
    return -(A[0]*B[1] + A[1]*B[0])/2 + ldot(A[2], B[2])
def Xij(xi, xj): return -2*edot(EX(xi), EX(xj))
def H(xi, zi, xj, zj):
    Xi, Zi, Xj, Zj = EX(xi), EZ(xi, zi), EX(xj), EZ(xj, zj)
    return -2*(edot(Zi, Zj)*edot(Xi, Xj) - edot(Zi, Xj)*edot(Zj, Xi))
def V(xi, zi, xj, xk):
    Xi, Zi, Xj, Xk = EX(xi), EZ(xi, zi), EX(xj), EX(xk)
    return (edot(Zi, Xj)*edot(Xi, Xk) - edot(Zi, Xk)*edot(Xi, Xj))/edot(Xj, Xk)

def flat_vectors(x, m=1.0):
    nu = x*m; s = np.sqrt(m*m - nu*nu)
    v1 = np.array([1j*nu/m, 0, 0, s/m]); v2 = np.array([0, 0, 1.0, 0])
    q1 = -(1j*nu/2)*v1 + (nu/2)*v2; q3 = (1j*nu/2)*v1 + (nu/2)*v2
    V2 = m*m*q1 + (m*m - nu*nu)*q3; V4 = (m*m - nu*nu)*q1 + m*m*q3
    W2 = m*m*q1 - (m*m - nu*nu)*q3; W4 = (m*m - nu*nu)*q1 - m*m*q3
    eta = np.diag([-1.0, 1, 1, 1]); fd = lambda a, b: a @ eta @ b
    nperp = np.array([0, -nu/m, s/m, 0]); npperp = np.array([0, nu/m, s/m, 0]); x3 = np.array([0, 0, 0, 1.0])
    return nu, s, (V2, W2, V4, W4), fd, (nperp, x3), (npperp, x3)

for x in (0.1, 0.3, 0.5, 0.7, 0.9):
    m = 1.0
    nu, s, (V2, W2, V4, W4), fd, B2, B4 = flat_vectors(x, m)
    a = (m - s)/nu
    Xc = np.array([-1j*a, 1j*a, 0.0])
    x4 = Xc/(-ldot(Xc, Xc)); x2 = -x4
    z0 = np.array([1.0, 0.0, 0.0], dtype=complex)
    # CFT polarisations from flat basis vectors: e2 = n_perp -> z2 = (0, -1, 0); e2 = x3 -> z2 = (0, 0, 1)   (Cartesian t, x, y)
    Z2 = [cart_to_lc(0, -1.0, 0.0), cart_to_lc(0, 0.0, 1.0)]
    Z4 = [cart_to_lc(0, -1.0, 0.0), cart_to_lc(0, 0.0, 1.0)]
    pr(f"x = {x}:")
    for sgn in (+1, -1):
        x0 = np.array([0.0, 0.0, sgn*a], dtype=complex)
        cft_P = np.array([[H(x2, Z2[i], x0, z0)*H(x4, Z4[j], x0, z0)/(Xij(x2, x0)*Xij(x4, x0)) for j in range(2)] for i in range(2)])
        cft_H = np.array([[H(x2, Z2[i], x4, Z4[j])/Xij(x2, x4) for j in range(2)] for i in range(2)])
        VV = np.array([[fd(B2[i], V2)*fd(B4[j], V4) for j in range(2)] for i in range(2)])
        WW = np.array([[fd(B2[i], W2)*fd(B4[j], W4) for j in range(2)] for i in range(2)])
        pr_P = VV/(4*nu**2*(m*m - nu*nu)); pr_H = -(VV + WW)/(nu**2*(m*m - nu*nu)**2)
        v042 = V(x0, z0, x4, x2)
        lam = 4*v042**2*Xij(x2, x4)/(Xij(x2, x0)*Xij(x4, x0))
        pr(f"  S0 point {sgn:+d}: (2V042)^2 X24/(X20 X40) = {lam:.10f}")
        pr(f"    H20H40/(X20X40) CFT  = {np.round(cft_P.flatten(), 10)}")
        pr(f"    printed image        = {np.round(pr_P.flatten(), 10)}   ratio CFT/printed = {np.round((cft_P/np.where(abs(pr_P) > 1e-14, pr_P, np.nan)).flatten(), 8)}")
        pr(f"    H24/X24 CFT          = {np.round(cft_H.flatten(), 10)}")
        pr(f"    printed image        = {np.round(pr_H.flatten(), 10)}   ratio CFT/printed = {np.round((cft_H/np.where(abs(pr_H) > 1e-14, pr_H, np.nan)).flatten(), 8)}")
open(os.path.join(HERE, 'cache', 'g4_images.txt'), 'w').write("\n".join(out) + "\n")
