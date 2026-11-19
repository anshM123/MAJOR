"""
s1 (OFFDIAG_KERNEL.md §1-2): EXACT leading-Regge heavy image of the helicity-CHANGING (++ <-> +-) matrix elements of C^_2' (q^2, j' = 0,
eta-kernel) and C^_3 (o_H q, j' = 1, f = -1), by the split representation of H7 §18.1 / S1_S4 §1.1 (v3_regge_kernel, q1_c3_regge_kernel
imported read-only).  v3/q1 computed only PRESERVING tensors (tau13 = eps_{h1} x eps_{-h1}, tau24 = eps_{h2} x eps_{-h2}); here the four
light-ray / dual polarisations are set independently:
    z1 = eps_{a1}(pbh), z3 = eps_{a3}(pbh)   [pair 1,3 at pbh],      z2' = eps_{a2}(ph), z4' = eps_{a4}(ph)   [pair 2,4 at ph]
(parallel-transported frame of v3.frames).  Preserving: a3 = -a1, a4 = -a2.  A helicity flip on a line is a3 = a1 (or a4 = a2).
Heavy 2x2 (parity-even, even J): ++ := {(h1,h2) = (+,-), (-,+)},  +- := {(+,+), (-,-)};  OPE vectors (lam_{(+,-)} + lam_{(-,+)})/sqrt2 etc.
F_{ab} = 1/2 sum_{(h1,h2) in a} sum_{(k1,k2) in b} L(a1=h1, a2=h2, a3=-k1, a4=-k2).  The off-diagonal therefore contains single flips only
(one line preserved, the other flipped), the diagonal contains preserved terms + double flips.
Normalisation (l9 units): divide by the forward (beta = 0) preserving image of the SAME functional (C2: h_- + h_+; C3: c3 = h0(++)), i.e.
the diagonal is P_l(cosh beta) (+ double flips), and Gamma_off/m := Re[F_off]/norm  (to be compared with r4's beta^2 g(nu beta)).
Usage: python3 s1_offdiag_kernel.py explore
"""
import os, sys, io, time, pickle, functools
os.environ.setdefault('OMP_NUM_THREADS', '2')
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..', 'h7i14'))
import numpy as np
import v3_regge_kernel as V
import q1_c3_regge_kernel as Q1

mdot = V.mdot
CACHE = os.path.join(HERE, 'cache')


def _pols(beta, pols):
    ph, pbh, ep_p, em_p, ep_pb, em_pb = V.frames(beta)
    E = {(+1, 'p'): ep_p, (-1, 'p'): em_p, (+1, 'pb'): ep_pb, (-1, 'pb'): em_pb}
    a1, a3, a2, a4 = pols
    return ph, pbh, E[(a1, 'pb')], E[(a3, 'pb')], E[(a2, 'p')], E[(a4, 'p')]


def _zgrid(nphi):
    phi = 2*np.pi*np.arange(nphi)/nphi
    return phi, np.stack([np.ones_like(phi), np.cos(phi), np.sin(phi)], -1)


def img2(nu, beta, pols, nphi=256, with_eta=True, a=7, Amonos=None):
    """C2-type (j' = 0) image with general polarisations; Amonos: (2,4) monomial list (default v3's q^2)."""
    delta = 0.5 + 1j*nu
    ph, pbh, e1, e3, e2, e4 = _pols(beta, pols)
    phi, z = _zgrid(nphi)
    b = a - 1 if with_eta else a
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    E2 = np.broadcast_to(e2, z.shape); E4 = np.broadcast_to(e4, z.shape)
    E1 = np.broadcast_to(e1, z.shape); E3 = np.broadcast_to(e3, z.shape)
    cA = (mdot(z, E2), mdot(z, E4), np.full(nphi, mdot(e2, e4)))
    cB = (mdot(z, E1), mdot(z, E3), np.full(nphi, mdot(e1, e3)))
    lsA = V.logKc(6 - 2 - delta, 1 - delta).real; lsB = V.logKc(6 - 3 + delta, delta).real
    Am = V.A_monomials(True) if Amonos is None else Amonos
    Ah = V.side_transform(P, z, E2, E4, Am, b, delta, cA, with_eta, lsA)
    Bh = V.side_transform(PB, z, E1, E3, V.B_monomials(True), b, delta, cB, with_eta, lsB)
    if with_eta:
        integ = -(Ah*Bh*V.ETA).sum(-1)
    else:
        integ = Ah*Bh
    return integ.mean()*2*np.pi


def img3(nu, beta, pols, nphi=256, b=8, kindA='C3', jp=1):
    delta = 0.5 + 1j*nu
    ph, pbh, e1, e3, e2, e4 = _pols(beta, pols)
    phi, z = _zgrid(nphi)
    P = np.broadcast_to(ph, z.shape); PB = np.broadcast_to(pbh, z.shape)
    E2 = np.broadcast_to(e2, z.shape); E4 = np.broadcast_to(e4, z.shape)
    E1 = np.broadcast_to(e1, z.shape); E3 = np.broadcast_to(e3, z.shape)
    cA = (mdot(z, E2), mdot(z, E4), np.full(nphi, mdot(e2, e4)))
    cB = (mdot(z, E1), mdot(z, E3), np.full(nphi, mdot(e1, e3)))
    wA = {'WE': Q1.wvec(z, e2), 'WF': Q1.wvec(z, e4)}
    wB = {'W1': Q1.wvec(z, e1), 'W3': Q1.wvec(z, e3)}
    lsA = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real; lsB = V.logKc(b - 3 + delta, delta).real
    Ah = Q1.side_transform_w(P, z, E2, E4, Q1.A3_monomials(kindA), b, delta, cA, wA, lsA)
    Bh = Q1.side_transform_w(PB, z, E1, E3, Q1.B_monomials_j(jp, True), b, delta, cB, wB, lsB)
    return (Ah*Bh).mean()*2*np.pi


PP = [(+1, -1), (-1, +1)]        # ++ class (h1, h2)
PM = [(+1, +1), (-1, -1)]        # +- class


def mat2x2(f, nu, beta, **kw):
    """F_{ab} = 1/2 sum_{(h1,h2) in a} sum_{(k1,k2) in b} L(a1=h1, a3=-k1, a2=h2, a4=-k2); also the 16 raw images."""
    raw = {}
    for (h1, h2) in PP + PM:
        for (k1, k2) in PP + PM:
            raw[(h1, h2, k1, k2)] = f(nu, beta, (h1, -k1, h2, -k2), **kw)
    F = np.zeros((2, 2), complex)
    for ia, A in enumerate((PP, PM)):
        for ib, B in enumerate((PP, PM)):
            F[ia, ib] = 0.5*sum(raw[(h1, h2, k1, k2)] for (h1, h2) in A for (k1, k2) in B)
    return F, raw


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'explore'
    if mode == 'explore':
        for nu in (20.0, 80.0):
            n2 = img2(nu, 0.0, (1, -1, -1, 1)); n3 = img3(nu, 0.0, (1, -1, -1, 1))
            for t in (1.0, 2.0, 3.7, 6.0):
                be = t/nu
                F2, r2 = mat2x2(img2, nu, be); F3, r3 = mat2x2(img3, nu, be)
                print(f"nu={nu} t={t} beta={be:.4f}")
                print("  C2 F/n2/beta^2:", np.round(F2/n2/be**2, 4).tolist())
                print("  C3 F/n3/beta^2:", np.round(F3/n3/be**2, 4).tolist())
                for k in sorted(r2):
                    print("   ", k, f"C2 {r2[k]/n2/be**2:.4f}  C3 {r3[k]/n3/be**2:.4f}")
