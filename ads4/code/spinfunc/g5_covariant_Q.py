"""
g5 (M2c): the covariant flip-free (2,4) structure, verified at the D = 4 crystal from CFT definitions (g4 helpers).
Claim: q := (2 H20 H40 + V_{0,42}^2 H24)/(X20 X40)  maps, under CLS-D's polarisation maps, to  -(W2W4 - V2V4)... i.e. to Q/(4 nu^2 (m^2-nu^2)^2) up to
sign, for every x, at both S^0 points (diagonal entries; the n_perp x x3 off-diagonal entries are odd under the S^0 reflection and cancel in
the j' = 0 sum).  Then the graviton (2,4) structure of C^_2 is q^2 (fixed covariant), i.e. in CLS-D's basis
S_n = (-2V_{0,42})^{J_F-4+2n} H20^{2-n} H40^{2-n} H24^n:  q^2 = [4 S_0 + S_1 + S_2/16] x (common prefactor).
Also re-derives that CLS-D's printed dual coefficients (16(d-1)/(d-2), 4/(d-2)) with the CORRECTED H20H40 image reproduce their printed K^(1).
Usage: python3 g5_covariant_Q.py -> cache/g5_covariant_Q.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
src = open(os.path.join(HERE, 'g4_images.py')).read().split("for x in (0.1")[0]
G = {'__file__': os.path.join(HERE, 'g4_images.py')}; exec(compile(src, 'g4h', 'exec'), G)
H, V, Xij, ldot, cart_to_lc, flat_vectors = G['H'], G['V'], G['Xij'], G['ldot'], G['cart_to_lc'], G['flat_vectors']
out = []
def pr(s=""):
    print(s, flush=True); out.append(str(s))
pr("# g5: point-by-point check: the g1 shock frame is one S0 point; the other is its x3-reflection (flips the shock helicity too)")
pr("# g5: covariant flip-free photon structure q = (2 H20 H40 + V042^2 H24)/(X20 X40) at the crystal vs Q = W W - V V")
for x in (0.1, 0.3, 0.5, 0.7, 0.9, 0.97):
    m = 1.0
    nu, s, (V2, W2, V4, W4), fd, B2, B4 = flat_vectors(x, m)
    a = (m - s)/nu; Xc = np.array([-1j*a, 1j*a, 0.0]); x4 = Xc/(-ldot(Xc, Xc)); x2 = -x4
    z0 = np.array([1.0, 0.0, 0.0], dtype=complex)
    Z = [cart_to_lc(0, -1.0, 0.0), cart_to_lc(0, 0.0, 1.0)]
    qs = []
    for sgn in (+1, -1):
        x0 = np.array([0.0, 0.0, sgn*a], dtype=complex)
        v = V(x0, z0, x4, x2)
        q = np.array([[(2*H(x2, Z[i], x0, z0)*H(x4, Z[j], x0, z0) + v**2*H(x2, Z[i], x4, Z[j]))/(Xij(x2, x0)*Xij(x4, x0))
                       for j in range(2)] for i in range(2)])
        qs.append(q)
    qavg = 0.5*(qs[0] + qs[1])
    QQ0 = np.array([[fd(B2[i], W2)*fd(B4[j], W4) - fd(B2[i], V2)*fd(B4[j], V4) for j in range(2)] for i in range(2)])/(4*nu**2*(m*m - nu*nu)**2)
    for sgn, qq in zip((+1, -1), qs):
        dev = np.max(np.abs(qq + QQ0))/np.max(np.abs(QQ0))
        pr(f"   x = {x}, S0 point {sgn:+d}: max|q + Q|/max|Q| (all four entries incl. off-diagonal) = {dev:.1e}")
    QQ = np.array([[fd(B2[i], W2)*fd(B4[j], W4) - fd(B2[i], V2)*fd(B4[j], V4) for j in range(2)] for i in range(2)])/(4*nu**2*(m*m - nu*nu)**2)
    ratio = qavg[np.abs(QQ) > 1e-12]/QQ[np.abs(QQ) > 1e-12]
    pr(f"x = {x}: q (S0-avg) diag = {np.round(np.diag(qavg), 10)}, offdiag = {np.round([qavg[0, 1], qavg[1, 0]], 12)}; "
       f"Q/(4nu^2(m^2-nu^2)^2) = {np.round(QQ.flatten(), 10)}; ratio q/Q on nonzero entries = {np.round(ratio, 10)}")
    # dual-coefficient consistency: K1 ~ 16(d-1)/(d-2) s1 + 4/(d-2) s2, s1 -> P-hat, s2 -> lam * H-hat (corrected image)
    d = 3
    Pc = np.array([[H(x2, Z[i], x0, z0)*H(x4, Z[j], x0, z0)/(Xij(x2, x0)*Xij(x4, x0)) for j in range(2)] for i in range(2)])
    Hc = np.array([[H(x2, Z[i], x4, Z[j])/Xij(x2, x4) for j in range(2)] for i in range(2)])
    lam = 4*v**2*Xij(x2, x4)/(Xij(x2, x0)*Xij(x4, x0))
    K1 = 16*(d - 1)/(d - 2)*Pc + 4/(d - 2)*lam*Hc/4      # (-2V)^2 H24 = 4 V^2 H24 and lam = 4V^2 X24/(X20X40)
    VVm = np.array([[fd(B2[i], V2)*fd(B4[j], V4) for j in range(2)] for i in range(2)]); WWm = np.array([[fd(B2[i], W2)*fd(B4[j], W4) for j in range(2)] for i in range(2)])
    Kp = (VVm + WWm - (4 - 2)*VVm)                       # printed K^(1) at D = 4 (up to overall factor)
    r = np.diag(K1)/np.diag(Kp)
    pr(f"        dual check: K1(from CLS-D duals + corrected images)/printed K^(1), diagonal entries: {np.round(r, 10)} (x-independent => consistent)")
open(os.path.join(HERE, 'cache', 'g5_covariant_Q.txt'), 'w').write("\n".join(out) + "\n")
