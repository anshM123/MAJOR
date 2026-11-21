"""Result Q, Stage 5 (RESULT_Q_STAGE5.md §1.3): INDEPENDENT numerical check of the exact rank-2 R_j (stage5_hj) for the raw spinning
k = 5 image: q1's numerical split-representation engine (side_transform_w, the actual boosted frames v3.frames(beta), 256-point z0
circle), with the extra x_mu x_nu realised by central finite differences in p (Richardson h, 2h) at fixed polarisations / w-vectors,
contracted with eta x eta, x 16.  Compared with N_q1(nu) sum_j R_j(nu) P_{l+j}(cosh beta) (N_q1 = 2 pi K_A^ref K_B^ref e^{-lsA-lsB},
verified exactly at rank 0 against q1's c3).  Also: free 5-harmonic fit (j = -2..2) at fixed nu: |h_{+-1}| / |h_0| (must vanish),
both helicity classes, and the light-ray-line flip.
Output: cache/stage5_check.txt
"""
import os, sys, io, time
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14')); sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
import numpy as np
import v3_regge_kernel as V
import q1_c3_regge_kernel as Q1
import stage5_hj as S5

mdot = V.mdot
ETA = np.array([-1.0, 1.0, 1.0])


def side_val(p, z, e, f, side, delta, b=8):
    nphi = z.shape[0]
    P = np.broadcast_to(p, z.shape); E = np.broadcast_to(e, z.shape); F = np.broadcast_to(f, z.shape)
    c = (mdot(z, E), mdot(z, F), np.full(nphi, mdot(e, f)))
    if side == 'A':
        ls = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real
        return Q1.side_transform_w(P, z, E, F, Q1.A3_monomials('C3'), b, delta, c, {'WE': Q1.wvec(z, e), 'WF': Q1.wvec(z, f)}, ls)
    ls = V.logKc(b - 3 + delta, delta).real
    return Q1.side_transform_w(P, z, E, F, Q1.B_monomials_j(1, True), b, delta, c, {'W1': Q1.wvec(z, e), 'W3': Q1.wvec(z, f)}, ls)


def hess(p, z, e, f, side, delta, h):
    """lower-index Hessian d^2/dp^mu dp^nu at fixed (z, e, f, w): array (nphi, 3, 3)."""
    E3 = np.eye(3)
    H = np.zeros((z.shape[0], 3, 3), complex)
    f0 = side_val(p, z, e, f, side, delta)
    for mu in range(3):
        fp = side_val(p + h*E3[mu], z, e, f, side, delta); fm = side_val(p - h*E3[mu], z, e, f, side, delta)
        H[:, mu, mu] = (fp - 2*f0 + fm)/h**2
        for nu_ in range(mu + 1, 3):
            a = side_val(p + h*E3[mu] + h*E3[nu_], z, e, f, side, delta)
            bq = side_val(p + h*E3[mu] - h*E3[nu_], z, e, f, side, delta)
            c = side_val(p - h*E3[mu] + h*E3[nu_], z, e, f, side, delta)
            dd = side_val(p - h*E3[mu] - h*E3[nu_], z, e, f, side, delta)
            H[:, mu, nu_] = H[:, nu_, mu] = (a - bq - c + dd)/(4*h**2)
    return H


def L5_num(nu, beta, h1, h2, flipB=False, nphi=256, h=2e-3):
    delta = 0.5 + 1j*nu
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(beta)
    Ep = {+1: ep_p, -1: em_p}; Epb = {+1: ep_pb, -1: em_pb}
    e2, e4 = Ep[h2], Ep[-h2]
    e1, e3 = Epb[h1], (Epb[h1] if flipB else Epb[-h1])
    phi = 2*np.pi*np.arange(nphi)/nphi
    z = np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)
    out = []
    for hh in (h, 2*h):
        HA = hess(ph, z, e2, e4, 'A', delta, hh); HB = hess(pbh, z, e1, e3, 'B', delta, hh)
        integ = 16*np.einsum('kmn,m,n,kmn->k', HA, ETA, ETA, HB)
        out.append(integ.mean()*2*np.pi)
    return (4*out[0] - out[1])/3, abs(out[0] - out[1])/abs(out[0])


def Nq1(nu, b=8):
    delta = 0.5 + 1j*nu
    lsA = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real; lsB = V.logKc(b - 3 + delta, delta).real
    return complex(S5.Nref(nu))*np.exp(-lsA - lsB)


if __name__ == '__main__':
    buf = io.StringIO()
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); buf.write(s + "\n")
    t0 = time.time()
    F = {cls: S5.numfuncs(f'C5_{cls}') for cls in ('pp', 'pm')}
    pr("# stage5_check: numerical (finite-difference, q1 engine) vs exact rank-2 R_j; ratio num/exact - 1 (Richardson spread in brackets)")
    for nu in (0.4, 1.3, 3.0, 7.5):
        for beta in (0.3, 0.9, 1.8):
            row = []
            for cls, (h1, h2) in (('pp', (1, -1)), ('pm', (1, 1))):
                num, spr = L5_num(nu, beta, h1, h2)
                ex = Nq1(nu)*sum(complex(F[cls][j](nu + 0j))*V.Pcon(nu - 1j*j, beta) for j in F[cls])
                row.append(f"{cls}: {abs(num/ex - 1):.1e} [{spr:.0e}]")
            pr(f"  nu = {nu:4.1f}, beta = {beta:3.1f}: " + "; ".join(row))
    pr("# free fit to j = -2..2 harmonics (++), 9 betas in [0.2, 2.4]: |h_j|/|h_0| for j = +-1 and fitted/exact for j = 0, +-2")
    betas = np.linspace(0.2, 2.4, 9)
    for nu in (0.7, 2.2):
        vals = np.array([L5_num(nu, bb, 1, -1)[0] for bb in betas])
        M = np.array([[V.Pcon(nu - 1j*j, bb) for j in (-2, -1, 0, 1, 2)] for bb in betas])
        sol, *_ = np.linalg.lstsq(M, vals, rcond=None)
        res = np.linalg.norm(M @ sol - vals)/np.linalg.norm(vals)
        exv = {j: Nq1(nu)*complex(F['pp'][j](nu + 0j)) for j in (-2, 0, 2)}
        pr(f"  nu = {nu}: resid {res:.1e}; |h_-1|/|h_0| {abs(sol[1])/abs(sol[2]):.1e}, |h_+1|/|h_0| {abs(sol[3])/abs(sol[2]):.1e}; "
           f"fit/exact - 1: j=-2 {abs(sol[0]/exv[-2] - 1):.1e}, j=0 {abs(sol[2]/exv[0] - 1):.1e}, j=2 {abs(sol[4]/exv[2] - 1):.1e}")
    pr("# light-ray-line flip (z1 = z3) at k = 5: |L_flip| / |L_preserving|")
    for nu, beta in ((1.3, 0.9), (3.0, 1.8)):
        lf, _ = L5_num(nu, beta, 1, -1, flipB=True); lp, _ = L5_num(nu, beta, 1, -1)
        pr(f"  nu = {nu}, beta = {beta}: {abs(lf)/abs(lp):.1e}")
    pr(f"[{time.time() - t0:.0f}s]")
    open(os.path.join(HERE, 'cache', 'stage5_check.txt'), 'w').write(buf.getvalue())
