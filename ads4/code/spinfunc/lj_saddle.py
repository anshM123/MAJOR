"""
lj_saddle (HLJ_ANALYTIC.md §1, §3): travelling-wave saddles of the exact CLS-D scalar crystal integrand at large J.

Exact structure used [P, lj_struct.py]: the scalar integrand of l1_crystal.build is
    A^{i nu/2} B^{Delta/2} e^{G_0(z)} P_J(t(z)) / beta_{Delta,J}   (times z-independent factors),
with G_0 the J = 0 slow exponent (the J-dependent 3-pt powers cancel against (Z^2 W^2)^{J/2} identically: Z^2 = -m15 m25 m12,
W^2 = -m35 m45 m34), t = Z.W/sqrt(Z^2 W^2).  Split P_J = W_+ + W_- (travelling waves, W_+- ~ e^{+-i(J+1/2)theta}/sqrt(2 pi J sin theta)
up to constant phases, theta = arccos t).  Each wave has the exponent
    Phi_+-(z; nu) = (i nu/2) log A + (Delta/2) log B + G_0 +- i(J+1/2) theta(z) - (1/2) log sin theta(z).
Saddle z_+-(nu) of Phi_+- (Newton, continued in J from the J = 0 crystal, then in nu).  Envelope theorem:
    d Phi_+-(z_+-(nu), nu)/d nu = (i/2) log A(z_+-)  =: i kappa_+-(nu),
so Re kappa_+- is BOTH the local nu-frequency of the wave and its exponential decay rate in Im nu (Cauchy-Riemann):
the wave is in the analyticity/decay class AG of HLJ_ANALYTIC §2 to all orders of its saddle expansion wherever z_+-(nu) is
non-degenerate and Re kappa_+ >= beta_* > 0 (Re kappa_- <= -beta_*).
Usage: python3 lj_saddle.py Delta J nu1,nu2,...   (prints kappa, det Hess, min |eig|, t, comparison with flat/LR frequencies)
"""
import os, sys
os.environ.setdefault('L1_DEG', '2')
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import l1_tps as T
import l1_crystal as Cr

assert T.DEG >= 2


def tvar(zv):
    """t = Z.W/sqrt(Z^2 W^2) (l1_crystal.build geometry) as TPS; also Z^2, W^2 and the six msq for the identity check."""
    xt0, ux, vx, uy, vy, xp, yp = zv
    x = [ux, vx, xp]; y = [uy, vy, yp]
    xsq = Cr.dot(x, x); rx = T.recip(-xsq)
    x4 = Cr.scal(rx, x); x2 = Cr.neg(x4); x1 = Cr.neg(y); x3 = y
    e = Cr.cpt(1.0, 1.0, 0.0)
    def msq(p, q):
        r = Cr.sub(p, q); return -Cr.dot(r, r)
    m = dict(m12=msq(x1, x2), m34=msq(x3, x4), m15=msq(x1, e), m25=msq(x2, e), m35=msq(x3, e), m45=msq(x4, e))
    x15, x25, x53, x54 = Cr.sub(x1, e), Cr.sub(x2, e), Cr.sub(e, x3), Cr.sub(e, x4)
    Z = Cr.sub(Cr.scal(m['m25'], x15), Cr.scal(m['m15'], x25))
    W = Cr.sub(Cr.scal(m['m45'], x53), Cr.scal(m['m35'], x54))
    ZZ, WW, ZW = Cr.dot(Z, Z), Cr.dot(W, W), Cr.dot(Z, W)
    iq = T.powr(T.mul(ZZ, WW), -0.5)
    if QREF[0] is not None and abs(-iq[0] - QREF[0]) < abs(iq[0] - QREF[0]):
        iq = -iq                                           # branch of (Z^2 W^2)^{-1/2} continued along a path (QREF)
    QLAST[0] = iq[0]
    t = T.mul(ZW, iq)
    return t, ZZ, WW, m


TREF = [None]          # (unused, kept for callers)
QREF = [None]          # reference value of (Z^2 W^2)^{-1/2} for continuation along paths
QLAST = [None]
SREF = [None]          # reference value of theta for continuation along complex-nu paths (None: principal branch)


def phi(zv, nu, Delta, J, sgn, slow=True):
    """Phi_sgn as TPS; returns (Phi, lA, t, theta)."""
    lamF, G, parts = Cr.build(zv, nu, Delta, 0, parts=True)
    t, _, _, _ = tvar(zv)
    one = T.const(1.0)
    s = T.powr(one - T.mul(t, t), 0.5)                     # sin theta (principal branch: theta in (0, pi) on the real crystal)
    th = -1j*T.log(t + 1j*s)                               # theta = arccos t
    if SREF[0] is not None:                                # analytic continuation along a path: the branch of theta (theta ->
        cands = []                                         # +-theta + 2 pi k, sin -> +-sin) closest to the reference SREF = theta_ref
        for sg in (1, -1):
            for k in (-1, 0, 1):
                cands.append((abs(sg*th[0] + 2*np.pi*k - SREF[0]), sg, k))
        _, sg, k = min(cands)
        th = sg*th; th[0] += 2*np.pi*k; s = sg*s
    P = lamF + G
    if J is not None:
        P = P + sgn*1j*(J + 0.5)*th
        if slow:
            P = P - 0.5*T.log(s)
    return P, parts['lA'], t, th


def newton(z, nu, Delta, J, sgn, it=60, tol=1e-11):
    for k in range(it):
        zv = [T.var(i, z[i]) for i in range(T.NV)]
        P, _, _, _ = phi(zv, nu, Delta, J, sgn)
        g, H = T.grad_hess(P)
        dz = np.linalg.solve(H, g)
        # damped step
        lam = 1.0
        if np.abs(dz).max() > 0.2*max(1e-3, np.abs(z).max()):
            lam = 0.2*max(1e-3, np.abs(z).max())/np.abs(dz).max()
        z = z - lam*dz
        if np.abs(dz).max() < tol*max(1.0, np.abs(z).max()):
            break
    else:
        return z, False
    return z, True


def info(z, nu, Delta, J, sgn):
    zv = [T.var(i, z[i]) for i in range(T.NV)]
    P, lA, t, th = phi(zv, nu, Delta, J, sgn)
    g, H = T.grad_hess(P)
    ev = np.linalg.eigvals(H)
    return dict(kappa=0.5*lA[0], Phi=P[0], t=t[0], theta=th[0], gnorm=np.abs(g).max(), detH=np.linalg.det(H),
                mineig=np.min(np.abs(ev)), eig=np.sort(np.abs(ev)))


def track_J(nu, Delta, J, sgn, nstep=None):
    """continue the saddle from the J = 0 crystal (exponent without the spin term) to spin J at fixed (nu, Delta)."""
    zs, a = Cr.saddle(nu, Delta)
    z = zs.copy()
    nstep = nstep or max(8, int(J/4))
    ok = True
    for Jc in np.linspace(0, J, nstep + 1)[1:]:
        z, ok = newton(z, nu, Delta, Jc, sgn)
        if not ok:
            return z, False, Jc
    return z, ok, J


def track_nu(z, nus, Delta, J, sgn):
    out = []
    for nu in nus:
        z, ok = newton(z, nu, Delta, J, sgn)
        out.append((nu, z.copy(), ok))
        if not ok:
            break
    return out


def flat_freq(nu, Delta, J):
    """d/dnu of (J+1/2) theta_flat, theta_flat = 2 arcsin(nu/m); and the LR frequency beta."""
    alpha = Delta - J - 2.0; abar = alpha + 2*J + 1
    m = np.sqrt(alpha*abar)
    beta = np.arccosh(1 + (2*J + 1)**2/(2*m*m))
    return (2*J + 1)/np.sqrt(m*m - nu*nu), beta, m


if __name__ == '__main__':
    Delta = float(sys.argv[1]); J = int(sys.argv[2])
    nus = [float(v) for v in sys.argv[3].split(',')]
    sgn = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    z, ok, Jr = track_J(nus[0], Delta, J, sgn)
    print('J-continuation at nu=%g: ok=%s (reached J=%g)' % (nus[0], ok, Jr))
    for nu, z, ok in track_nu(z, nus, Delta, J, sgn):
        I = info(z, nu, Delta, J, sgn)
        ff, beta, m = flat_freq(nu, Delta, J)
        print('nu=%8.3f ok=%d  kappa=%+.6f%+.6fi  flat=%.6f beta=%.6f  t=%.6f%+.2ei  |g|=%.1e  min|eig|=%.3g  eigs=%s' % (
            nu, ok, I['kappa'].real, I['kappa'].imag, ff, beta, I['t'].real, I['t'].imag, I['gnorm'], I['mineig'],
            np.array2string(I['eig'][:4], precision=3)))
