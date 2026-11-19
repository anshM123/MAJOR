"""
q1b (S4 fix, S1_S4_FIXES.md §1.3): end-to-end check of the q1 C^_3 Regge kernel against the independent crystal integrand (l9_points C3H vs
C2 at the heavy reference Delta_ref = 2e4), the analogue of H7_I14 §18.2 check X.

Crystal side: at Delta_ref (m_ref = 2e4, so O(1/m^2) effects ~ 1e-9), the C3H/C2 ratio of the ++ element (l9 helicity basis, K2 linear
chain) is  R_cr(nu) = [N3(nu)/N2(nu)] * (-2i x^2/(2 - x^2)) (1 + O(x^2/m^2)),  x = nu/m_ref (H0 §23.3).  So
    Q_cr(nu) := R_cr (2 - x^2)/(-2i x^2) = N3/N2 (nu-only).
Kernel side (q1 + v3 engines, identical nu-only rescalings on both functionals): forward (beta = 0) images
    Q_k(nu) := L3(0) rho_B(nu) / (L2(0) nu^2),  rho_B = B(1/2, 1 - delta/2)/B(1/2, 1/2 - delta/2)
(rho_B restores the theta-moment normalisation of the j' = 1 side relative to v3's j' = 0 one), with
  L2 = - eta^{mu nu} dA/dp dB/dpb of the q^2 / j' = 0 eta-kernel (b = 6, v3) and L3 = the o q / j' = 1 / f = -1 image (b = 8, q1).
Both are nu-only up to a nu-independent constant, so log Q = c + p log nu + a1/nu + a2/nu^2 + ... must agree in p, a1, a2 (c differs by the
dropped constants).  Fit window nu in [40, 1000] (crystal: K2 truncation small there; H0 §22.1 used the same window for C2).
Output: cache/q1b_c3_crystal_check.txt
"""
import os, sys, io
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
from scipy.special import loggamma
import q1_c3_regge_kernel as Q
import v3_regge_kernel as V

DREF = 20000.0


def crystal(J=4, K=2):
    import l9_analysis as A
    L = A.L
    R2, R3 = A.Recs('C2', J, K).R, A.Recs('C3H', J, K).R
    Rb = dict(L.records(J)); Rb.update(A._load('C2', J, '_k2')); Rb.update(A._load('C2', J, ''))
    B = L.basis(J, DREF, Rb)
    mr = A.mphys(DREF, J)
    out = []
    for k, r3 in R3.items():
        if k[0] != DREF:
            continue
        c = [kk for kk in R2 if kk[0] == DREF and abs(kk[1] - k[1]) < 1e-7*k[1]]
        if not c:
            continue
        r2 = R2[c[0]]
        M3, M2 = A.helmat(r3, B, K), A.helmat(r2, B, K)
        nu = k[1]; x = nu/mr
        v = np.exp(r3['logC'] - r2['logC'])*A.pred([Mo[0, 0] for Mo in M3], [Mo[0, 0] for Mo in M2], K)
        out.append((nu, v*(2 - x*x)/(-2j*x*x)))
    out.sort()
    return np.array([o[0] for o in out]), np.array([o[1] for o in out])


def kernel_forward(nu, b2=6, b3=8):
    delta = 0.5 + 1j*nu
    beta = 0.0
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(beta)
    nphi = 256
    phi = 2*np.pi*np.arange(nphi)/nphi
    z = np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)
    h1, h2 = 1, -1
    E = {(+1, 'p'): ep_p, (-1, 'p'): em_p, (+1, 'pb'): ep_pb, (-1, 'pb'): em_pb}
    e2, e4, e1, e3 = E[(h2, 'p')], E[(-h2, 'p')], E[(h1, 'pb')], E[(-h1, 'pb')]
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    E2, E4, E1, E3 = (np.broadcast_to(v, z.shape) for v in (e2, e4, e1, e3))
    cA = (V.mdot(z, E2), V.mdot(z, E4), np.full(nphi, V.mdot(e2, e4)))
    cB = (V.mdot(z, E1), V.mdot(z, E3), np.full(nphi, V.mdot(e1, e3)))
    lsA = V.logKc(6 - delta, 1 - delta).real; lsB = V.logKc(5 + delta, delta).real      # common to both functionals
    Ah = V.side_transform(P, z, E2, E4, V.A_monomials(True), b2, delta, cA, True, lsA)
    Bh = V.side_transform(PB, z, E1, E3, V.B_monomials(True), b2, delta, cB, True, lsB)
    L2 = (-(Ah*Bh*V.ETA).sum(-1)).mean()*2*np.pi
    A3 = Q.side_transform_w(P, z, E2, E4, Q.A3_monomials('C3'), b3, delta, cA, {'WE': Q.wvec(z, e2), 'WF': Q.wvec(z, e4)}, lsA)
    B3 = Q.side_transform_w(PB, z, E1, E3, Q.B_monomials_j(1, True), b3, delta, cB, {'W1': Q.wvec(z, e1), 'W3': Q.wvec(z, e3)}, lsB)
    L3 = (A3*B3).mean()*2*np.pi
    lrho = (loggamma(1 - delta/2) - loggamma(1.5 - delta/2)) - (loggamma(0.5 - delta/2) - loggamma(1 - delta/2))
    return L3*np.exp(lrho)/(L2*nu*nu)


def fit(nus, Qv, npow=3, free_p=False):
    lq = np.unwrap(np.angle(Qv))
    y = np.log(np.abs(Qv)) + 1j*lq
    cols = [np.ones_like(nus)] + ([np.log(nus)] if free_p else []) + [nus**-k for k in range(1, npow + 1)]
    A_ = np.array(cols).T
    c, *_ = np.linalg.lstsq(A_.astype(complex), y, rcond=None)
    res = np.abs(A_ @ c - y).max()
    return c, res


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(s=""):
        print(s, flush=True); buf.write(s + "\n")
    nus_c, Qc = crystal(4)
    sel = (nus_c >= 40) & (nus_c <= 1000)
    pr(f"# q1b: crystal N3/N2 (J = 4, Delta_ref = 2e4, K2) at {sel.sum()} nodes nu in [40, 1000] vs the q1 Regge-kernel forward ratio")
    nk = np.geomspace(40, 1000, 24)
    Qk = np.array([kernel_forward(v) for v in nk])
    for lab, nus, Qv in (('crystal', nus_c[sel], Qc[sel]), ('kernel ', nk, Qk)):
        for fp in (True, False):
            c, res = fit(nus, Qv, 3, fp)
            co = list(c[1:])
            pstr = (f"p = {co.pop(0).real:+.4f}{0:+.0e}i, " if fp else "p := 0, ")
            pr(f"  {lab}: log Q = c + {pstr}a1 = {co[0]:.3f}, a2 = {co[1]:.1f}, a3 = {co[2]:.0f}   (max fit resid {res:.1e})")
    # direct comparison of the nu-dependence (constant removed at nu = 200)
    pr("  Q(nu)/Q(200) crystal vs kernel (kernel evaluated at the crystal nodes):")
    i0 = np.argmin(np.abs(nus_c - 200))
    k0 = kernel_forward(nus_c[i0])
    for i in np.where(sel)[0][::max(1, sel.sum()//12)]:
        kv = kernel_forward(nus_c[i])
        pr(f"    nu={nus_c[i]:8.2f}: crystal {Qc[i]/Qc[i0]:.6f}   kernel {kv/k0:.6f}   rel diff {abs((Qc[i]/Qc[i0])/(kv/k0) - 1):.1e}")
    open(os.path.join(HERE, 'cache', 'q1b_c3_crystal_check.txt'), 'w').write(buf.getvalue())
