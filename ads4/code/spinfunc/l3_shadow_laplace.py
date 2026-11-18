"""
l3 (M6 step 2, first piece): the shadow coefficient by Laplace expansion of the shadow integral, validated on the scalar beta_{Delta,J}.

For spinning externals the gauge-fixed CLS-D formula needs the shadow coefficient MATRIX of <T T O>, not known in closed form.
CLS-D (sec. shadow_coeff) obtain it at leading order from the saddle of the shadow integral at x5 = e, with x3 = y*, x4 = x*/(-x*^2) at
the crystal and x6 = (-e)^+:
     (1/beta) int d^3x5 <O^dag(x6,z6) O(x5)> |<O~^dag(x5) phi3 phi4>| = |<O^dag(x6,z6) phi3 phi4>| .
Here the x5 integral is expanded to relative O(1/Delta^2) with the l1 engine (3 live variables + 4 inert Gaussian dummies), and
|beta_Laplace| / |beta_exact| is tested for Delta-independence (J-only constant from the index-contraction normalisation) and for
independence of the crystal point (nu/Delta) used for x3, x4.  "+"-shifted separations enter as -x56^2 -> +|.| (X6^+ = -X6).
Usage: python3 l3_shadow_laplace.py  -> cache/l3_shadow_laplace.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
from l1_crystal import dot, sub, scal, cpt, log_beta

DF = 3.0


def build(zv, Delta, J, x3, x4):
    w0, w1, w2 = zv[0], zv[1], zv[2]
    e = cpt(1.0, 1.0, 0.0)
    x5 = [e[0] + w0, e[1] + w1, e[2] + w2]
    x6 = cpt(-1.0, -1.0, 0.0)
    X3 = [T.const(c) for c in x3]; X4 = [T.const(c) for c in x4]
    def msq(p, q):
        r = sub(p, q); return -dot(r, r)
    d56 = msq(x5, x6)                                   # = -x56^2 (= 4 at the saddle; spacelike to 6^+)
    m35, m45, m34 = msq(X3, x5), msq(X4, x5), msq(X3, X4)
    lamF = Delta*(-T.log(d56) + 0.5*T.log(m35) + 0.5*T.log(m45))
    for i in range(3, T.NV):                            # inert dummies
        lamF = lamF - 0.5*T.mul(zv[i], zv[i])
    G = -((3 + J)/2)*(T.log(m35) + T.log(m45)) - ((2*DF - 3 + Delta + J)/2)*T.log(m34)
    if J > 0:
        z6 = cpt(1.0, 0.0, 0.0)
        x56 = sub(x5, x6)
        r = T.div(dot(x56, z6), dot(x56, x56))
        a = sub(z6, scal(2*r, x56))                     # I(x56).z6  (null)
        W = sub(scal(m45, sub(x5, X3)), scal(m35, sub(x5, X4)))
        G = G + J*T.log(dot(a, W))
    return lamF, G


def rhs_log(Delta, J, x3, x4):
    """log |<O^dag(x6, z6) phi3 phi4>| at x6 = -e (spacelike-shifted)."""
    x3 = np.array(x3); x4 = np.array(x4); x6 = np.array([-1.0, -1.0, 0.0]); z6 = np.array([1.0, 0.0, 0.0])
    d = lambda p, q: -0.5*(p[0]*q[1] + p[1]*q[0]) + p[2]*q[2]
    m36 = d(x3 - x6, x3 - x6); m46 = d(x4 - x6, x4 - x6)          # '+': X6 -> -X6 flips the sign of x_{i6}^2
    m34 = -d(x3 - x4, x3 - x4)
    val = -((2*DF - Delta + J)/2)*np.log(m34 + 0j) - ((Delta + J)/2)*(np.log(m36 + 0j) + np.log(m46 + 0j))
    if J > 0:
        num = 2*d(z6, x4 - x6)*m36 - 2*d(z6, x3 - x6)*m46
        val += J*np.log(num + 0j)
    return val


def laplace_beta(Delta, J, xi):
    a = (1 - np.sqrt(1 - xi**2))/xi                    # crystal for nu/Delta = xi
    y = np.array([1j*a, -1j*a, 0.0]); x = np.array([-1j*a, 1j*a, 0.0])
    xsq = -x[0]*x[1]
    x3 = y; x4 = x/(-xsq)
    zs = np.zeros(T.NV, complex)
    zv = [T.var(i, 0.0) for i in range(T.NV)]
    lamF, G = build(zv, Delta, J, x3, x4)
    g, H = T.grad_hess(lamF)
    L = np.linalg.inv(T.csym_cholesky(-H)).T
    for _ in range(2):
        lamF2, G2 = build(T.linsub(zs, L), Delta, J, x3, x4)
        _, H2 = T.grad_hess(lamF2)
        L = L @ np.linalg.inv(T.csym_cholesky(-H2)).T
    h = [T.homog(lamF2, k) for k in range(T.DEG + 1)]; gg = [T.homog(G2, k) for k in range(T.DEG + 1)]
    A1 = h[3] + gg[1]; A2 = h[4] + gg[2]; A3 = h[5] + gg[3]; A4 = h[6] + gg[4]
    A1s = T.mul(A1, A1)
    K1 = T.expect(A2) + 0.5*T.expect(A1s)
    K2 = T.expect(A4) + 0.5*T.expect2(A2, A2) + T.expect2(A1, A3) + 0.5*T.expect2(A1s, A2) + T.expect2(A1s, A1s)/24.0
    logI0 = lamF[0] + G[0] + 1.5*np.log(2*np.pi) + np.log(np.linalg.det(L[:3, :3]))
    return dict(grad=np.abs(g[:3]).max(), K1=K1, K2=K2, logI0=logI0, rhs=rhs_log(Delta, J, x3, x4))


if __name__ == '__main__':
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr("# l3: |beta_Laplace|/|beta_exact| (shadow integral expanded around x5 = e to O(1/Delta^2)); should be a J-only constant")
    for J in (0, 2, 4):
        for xi in (0.4, 0.6):
            row = []
            for Dl in (100.0, 300.0, 1000.0, 3000.0):
                r = laplace_beta(Dl, J, xi)
                lb = [r['logI0'].real + np.log(abs(1 + c)) - r['rhs'].real for c in (0, r['K1'], r['K1'] + r['K2'])]
                row.append((Dl, r['grad'], [l - log_beta(Dl, J) for l in lb], r['K1']*Dl, r['K2']*Dl**2))
            ref = row[-1][2]
            pr(f"  J={J} nu/Delta={xi}: |grad| at x5=e max {max(rw[1] for rw in row):.1e}; K1*Delta = {row[-1][3].real:.4f}, "
               f"K2*Delta^2 = {row[-1][4].real:.3f}")
            for (Dl, gr, l, k1, k2) in row:
                pr(f"     Delta={Dl:6.0f}: log(|beta_L|/|beta|) - (Delta=3000 value), orders 0/1/2: "
                   f"{l[0]-ref[0]:+.2e} / {l[1]-ref[1]:+.2e} / {l[2]-ref[2]:+.2e};  constant (order 2) = {l[2]:.8f}")
    open(os.path.join(HERE, 'cache', 'l3_shadow_laplace.txt'), 'w').write("\n".join(out) + "\n")
