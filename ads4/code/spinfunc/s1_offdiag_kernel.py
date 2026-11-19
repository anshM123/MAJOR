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


if __name__ == '__main__' and (len(sys.argv) < 2 or sys.argv[1] == 'explore'):
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


# ------------------------------------------------------------------------------------------------ vectorised over beta (one nu)
def _frames_vec(betas):
    b = np.asarray(betas, float)[:, None]
    c, s = np.cosh(b/2), np.sinh(b/2)
    z0 = np.zeros_like(c)
    ph = np.concatenate([c, s, z0], -1); pbh = np.concatenate([c, -s, z0], -1)
    npp = np.concatenate([s, c, z0], -1); npb = np.concatenate([-s, c, z0], -1)
    ep = np.array([0.0, 0.0, 1.0])
    eps = lambda npar, h: (npar + 1j*h*ep)/np.sqrt(2)
    return ph, pbh, {(+1, 'p'): eps(npp, 1), (-1, 'p'): eps(npp, -1), (+1, 'pb'): eps(npb, 1), (-1, 'pb'): eps(npb, -1)}


def img_vec(nu, betas, pols, kind='C2', nphi=256, Amonos=None, Bmonos=None, jp=None):
    """images for all betas at one nu (same conventions as img2 / img3).  kind 'C2': j' = 0, eta-kernel (b = 6, gradient);
    kind 'C3': o_H q, j' = 1, f = -1 (b = 8).  Amonos: optional (2,4) monomial list in q1's format (with wkey) or v3's (without)."""
    delta = 0.5 + 1j*nu
    ph, pbh, E = _frames_vec(betas)
    a1, a3, a2, a4 = pols
    e1, e3, e2, e4 = E[(a1, 'pb')], E[(a3, 'pb')], E[(a2, 'p')], E[(a4, 'p')]
    phi, z1 = _zgrid(nphi)
    nb = len(betas)
    z = np.broadcast_to(z1[None], (nb, nphi, 3))
    bc = lambda v: np.broadcast_to(v[:, None, :], (nb, nphi, 3))
    P, PB, E1, E3, E2, E4 = bc(ph), bc(pbh), bc(e1), bc(e3), bc(e2), bc(e4)
    cA = (mdot(z, E2), mdot(z, E4), mdot(E2, E4))
    cB = (mdot(z, E1), mdot(z, E3), mdot(E1, E3))
    if kind == 'C2':
        lsA = V.logKc(6 - 2 - delta, 1 - delta).real; lsB = V.logKc(6 - 3 + delta, delta).real
        Am = V.A_monomials(True) if Amonos is None else Amonos
        if Am and len(Am[0]) == 6:                                    # q1-format monomials (w-vectors): no gradient support -> error
            raise ValueError('C2 with w-monomials not supported')
        Ah = V.side_transform(P, z, E2, E4, Am, 6, delta, cA, True, lsA)
        Bh = V.side_transform(PB, z, E1, E3, V.B_monomials(True), 6, delta, cB, True, lsB)
        integ = -(Ah*Bh*V.ETA).sum(-1)
    else:
        b = 8
        lsA = V.logKc(b - 3 + 1 - delta + 1, 1 - delta).real; lsB = V.logKc(b - 3 + delta, delta).real
        wA = {'WE': Q1.wvec(z, E2), 'WF': Q1.wvec(z, E4)}
        wB = {'W1': Q1.wvec(z, E1), 'W3': Q1.wvec(z, E3)}
        Am = Q1.A3_monomials('C3') if Amonos is None else Amonos
        Ah = Q1.side_transform_w(P, z, E2, E4, Am, b, delta, cA, wA, lsA)
        Bh = Q1.side_transform_w(PB, z, E1, E3, Q1.B_monomials_j(1 if jp is None else jp, True), b, delta, cB, wB, lsB)
        integ = Ah*Bh
    return integ.mean(-1)*2*np.pi


PRES = (1, -1, -1, 1)          # ++ preserving: z1 = eps_+, z3 = eps_-, z2' = eps_-, z4' = eps_+
FLIP24 = (1, -1, -1, -1)       # single flip on the (2,4) line: the only non-vanishing off-diagonal configuration (C2)


BETAS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.13, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.75, 0.9, 1.1, 1.35, 1.7, 2.1]


def nu_grid(beta, numax=3000.0, tmax=100.0):
    a = np.arange(0.01, 8.0, 0.04)
    t = np.arange(8.0*beta, min(tmax, numax*beta), 0.12)
    return np.concatenate([a, t[t/beta > 8.0]/beta])


def table(betas, tag, Amonos=None, kind='C2', pols=FLIP24, nphi=192):
    """Y(nu, beta) := L_off(nu, beta)/L_pres(nu, 0) (forward-normalised off-diagonal; Re = Gamma_off/m in the geometric basis)."""
    out = {}
    t0 = time.time()
    for be in betas:
        nus = nu_grid(be)
        Y = np.zeros(len(nus), complex)
        for k, nu in enumerate(nus):
            v = img_vec(nu, np.array([0.0, be]), pols, kind=kind, nphi=nphi, Amonos=Amonos)
            n0 = img_vec(nu, np.array([0.0]), PRES, kind=kind, nphi=nphi, Amonos=Amonos)[0] if pols != PRES else v[0]
            Y[k] = v[1]/n0
        out[be] = (nus, Y)
        print(f"  [{tag}] beta={be}: {len(nus)} nodes  [{time.time() - t0:.0f}s]", flush=True)
        pickle.dump(out, open(os.path.join(CACHE, f's1_table_{tag}.pkl'), 'wb'))
    return out


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'table':
    idx = [int(v) for v in sys.argv[3].split(',')] if len(sys.argv) > 3 else range(len(BETAS))
    table([BETAS[i] for i in idx], sys.argv[2])


# ------------------------------------------------------------------------------------------------ other (2,4) structures at k = 2 (j' = 0)
import sympy as sp
_X, _s, _ex, _fx = V._X, V._s, V._ex, V._fx
_a2, _a4, _e24, _dl = V._a2, V._a4, V._e24, V._dl
# building blocks at x0 = 0 (H7 §18.1 reduction; z' inversion-referred):  H20 = a2, H40 = a4, V0 = s/X = -(z0.xh), V2 = ex/X, V4 = fx/X,
# H24^ = z2'.I(xh).z4' = e24 + 2 V2 V4.  q = 2 H20 H40 - V0^2 H24^ (v3).  P := H20 H40,  Hh := V0^2 H24^.
_P = _a2*_a4
_Hh = _s**2*_X**-2*(_e24 + 2*_ex*_fx*_X**-2)
_V0, _V2, _V4 = _s/_X, _ex/_X, _fx/_X
STRUCTS = {
    'q2': (2*_P - _Hh)**2,                        # route A' (C^_2')
    'PP': _P**2, 'PH': _P*_Hh, 'HH': _Hh**2,      # symmetric j' = 0 span (crystal: a P^2 + b P H + c H^2, flip-free iff a + 2b + 4c = 0)
    'Nq': _Hh*(2*_P - _Hh),                       # = 2 PH - HH: flip-free AND null at the bulk point (a = 0)
    'A1': _a2**2*_V0**2*_V4**2,                   # 2<->4-asymmetric H20^2 V0^2 V4^2
    'A1m': _a4**2*_V0**2*_V2**2,                  # its mirror
    'A2': _P*_V0**2*_V2*_V4,                      # H20 H40 V0^2 V2 V4
    'A3': _a2**2*_a4*_V0*_V4,                     # H20^2 H40 V0 V4  (v3's "H20^2 H40 V4")
    'A3m': _a4**2*_a2*_V0*_V2,
    'VV': _V0**4*_V2**2*_V4**2,                   # V0^4 V2^2 V4^2
    'EE': _V0**4*_e24**2,                         # V0^4 (z2'.z4')^2 = V0^4 (H24^ - 2 V2 V4)^2: vanishes on z2' = z4' (Regge flip-free)
    'EH': _V0**4*_e24*(_e24 + 2*_ex*_fx*_X**-2),  # V0^4 e24 H24^
}


@functools.lru_cache(maxsize=None)
def A_monos_struct(name):
    expr = sp.expand(16*_X**(1 - _dl + 4)*_s**(-(1 - _dl) - 4)*sp.expand(STRUCTS[name]))
    out = []
    for term in sp.Add.make_args(expr):
        pw = term.as_powers_dict()
        k = sp.simplify(pw.get(_X, 0)); m = sp.simplify(pw.get(_s, 0)); i = int(pw.get(_ex, 0)); j = int(pw.get(_fx, 0))
        assert i <= 2 and j <= 2, (name, term)
        coef = sp.simplify(term/(_X**k*_s**m*_ex**i*_fx**j))
        out.append((sp.lambdify((_a2, _a4, _e24, _dl), coef, 'numpy'), sp.lambdify(_dl, k, 'numpy'),
                    sp.lambdify(_dl, m, 'numpy'), i, j))
    return out
