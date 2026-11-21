"""
pw_grav (HWIN_THIMBLE.md): EXACT per-wave spin-J contraction of the graviton crystal integrand (kinds 'EE' = w0 EE integrand, 'C3H').

The <T T O_J> numerators of l5_spin.tto_samples are, as polynomials of the O polarisation z5 (degree J),
    F_a(z_b; z5) = (v.z5)^{J-4} T_{a}(z_b; z5),      T_a: a quartic form in z5,
with v.z5 = V_{5,ab}(z5) (plain component sums u.z = sum_k u_k z_k), because every monomial V5^{a5} H15^{b15} H25^{b25} has
b15 + b25 <= 4 (degree counting in z_a, z_b).  With z^alpha (v.z)^{J-4} = [(J-4)!/J!] d_U^alpha (U.z)^J |_{U=v} (|alpha| = 4) and the
spin-J pairing (l5_spin: <T, S> = sum_n T_n S_{-n}/kappa_n) of two J-th powers
    G_f(U, U') := < (U.z)^J, (U'.z)^J >  =  c_J (g(U,U) g(U',U'))^{J/2} f(g(U,U')/sqrt(g(U,U) g(U',U'))),   f = P_J,
c_J = 2^J J!^2/(2J)!, g(u, u') = -2(u0 u1' + u1 u0') + u2 u2' (the inverse lightcone metric on plain coefficients), the spin-J pairing is
    <F_a(z_b), G_c(z_d)> = [(J-4)!/J!]^2 sum_{|alpha| = |beta| = 4} T_{a,alpha} T'_{c,beta} d_U^alpha d_U'^beta G_f(v, v').
Everything is linear in f, so f = W_+- gives the two travelling-wave parts exactly (P_J = W_+ + W_-); the derivatives are Taylor
coefficients of a bidegree-(4,4) truncated power series in (dU, dU') around (v, v') (f composed through its t-derivatives to order 8,
pw_wave.wave_derivs).  No Fourier sum over the 2J+1 spin-J modes (so no eps^2 4^J round-off, FOLD_WINDOW §1.1).
The spin-2 indices (z_b, z_d on the 5-point null circles) and the D = E (x) E (EE) / q (x) o_H (C3H) structures are contracted exactly as
in w0_ee_common._contract / l9_points.build_spin.  The large scalar c_J (g g')^{J/2} f(t0) [(J-4)!/J!]^2 goes into the exponent.
The wave label uses t0 = sigma t_scalar with the scalar's continued sqrt (pw_core), sigma = the fixed sign relating the two (checked).
"""
import os, sys, itertools
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import scipy.sparse as sp
from scipy.special import gammaln
import n1_pw as PW
import pw_wave as WV
import w0_ee_common as W0

# ------------------------------------------------------------------ bidegree-(4,4) TPS in (dU, dU') (3 + 3 variables)
DG = 4
MON3 = [e for d in range(DG + 1) for e in itertools.product(range(d + 1), repeat=3) if sum(e) == d]
MON3 = sorted(set(MON3), key=lambda e: (sum(e), e))
NM = len(MON3)                                      # 35
IDX3 = {e: i for i, e in enumerate(MON3)}
DEG3 = np.array([sum(e) for e in MON3])
import math
FACT3 = np.array([math.prod(math.factorial(k) for k in e) for e in MON3], float)
_P3 = [(i, j, IDX3[tuple(a + b for a, b in zip(MON3[i], MON3[j]))]) for i in range(NM) for j in range(NM)
       if DEG3[i] + DEG3[j] <= DG]
P3I = np.array([p[0] for p in _P3]); P3J = np.array([p[1] for p in _P3]); P3K = np.array([p[2] for p in _P3])
NT = NM*NM
_I = (P3I[:, None]*NM + P3I[None, :]).ravel(); _Jx = (P3J[:, None]*NM + P3J[None, :]).ravel()
_K = (P3K[:, None]*NM + P3K[None, :]).ravel()
_PM = sp.csr_matrix((np.ones(len(_K)), (_K, np.arange(len(_K)))), shape=(NT, len(_K)))
_PM3 = sp.csr_matrix((np.ones(len(P3K)), (P3K, np.arange(len(P3K)))), shape=(NM, len(P3K)))
TOP = np.where(DEG3 == DG)[0]                       # the 15 quartic monomials


def mul(a, b):
    """bidegree TPS product; a, b: (NT, n)."""
    return _PM @ (a[_I]*b[_Jx])


def mul3(a, b):
    """single-side (35, n) product."""
    return _PM3 @ (a[P3I]*b[P3J])


def series3(A, coefs):
    """sum_k coefs[k] A^k for a nilpotent single-side TPS A (no constant term), k <= DG."""
    out = np.zeros_like(A); out[0] = coefs[0]
    p = None
    for k in range(1, DG + 1):
        p = A.copy() if p is None else mul3(p, A)
        out = out + coefs[k]*p
    return out


def outer(a, b):
    """U-only (35, n) times U'-only (35, n) -> (NT, n)."""
    return (a[:, None, :]*b[None, :, :]).reshape(NT, -1)


def binom_coefs(p, K=DG):
    c = [1.0]
    for k in range(1, K + 1):
        c.append(c[-1]*(p - k + 1)/k)
    return c


def gmet(u, w):
    return -2*(u[0]*w[1] + u[1]*w[0]) + u[2]*w[2]


def g_tps(v, vp, fder):
    """Taylor coefficients (NT, n) of G_f(v + dU, v' + dU') / [c_J (g g')^{J/2} f(t0)], with fder = (f^{(k)}(t0)/f(t0), k = 0..8)."""
    n = v.shape[1]
    J2 = fder['J']/2.0
    a0 = gmet(v, v); b0 = gmet(vp, vp); c0 = gmet(v, vp)
    # single-side quadratic forms a(dU)/a0 - 1 etc.
    A = np.zeros((NM, n), complex); B = np.zeros((NM, n), complex)
    # g(v + dU, v + dU) = a0 + 2 g(v, dU) + g(dU, dU)
    E = np.eye(3)
    for k in range(3):
        A[IDX3[tuple(E[k].astype(int))]] = 2*gmet(v, E[k][:, None]*np.ones((1, n)))/a0
        B[IDX3[tuple(E[k].astype(int))]] = 2*gmet(vp, E[k][:, None]*np.ones((1, n)))/b0
    for k in range(3):
        for l in range(k, 3):
            e = tuple((E[k] + E[l]).astype(int))
            gk = gmet(E[k][:, None], E[l][:, None])[0]
            fac = 1.0 if k == l else 2.0
            A[IDX3[e]] += fac*gk/a0; B[IDX3[e]] += fac*gk/b0
    # C = g(v + dU, v' + dU') / c0 - 1: linear + bilinear
    Cc = np.zeros((NT, n), complex)
    for k in range(3):
        ik = IDX3[tuple(E[k].astype(int))]
        Cc[ik*NM + 0] += gmet(E[k][:, None]*np.ones((1, n)), vp)/c0
        Cc[0*NM + ik] += gmet(v, E[k][:, None]*np.ones((1, n)))/c0
        for l in range(3):
            il = IDX3[tuple(E[l].astype(int))]
            Cc[ik*NM + il] += gmet(E[k][:, None], E[l][:, None])[0]/c0
    one3 = np.zeros((NM, n), complex); one3[0] = 1.0
    Ah = series3(A, binom_coefs(-0.5)); Bh = series3(B, binom_coefs(-0.5))       # (1 + A)^{-1/2}
    AJ = series3(A, binom_coefs(J2)); BJ = series3(B, binom_coefs(J2))           # (1 + A)^{J/2}
    Cfull = Cc.copy(); Cfull[0] += 1.0
    tt = mul(Cfull, outer(Ah, Bh))                       # t/t0
    t0 = fder['t0']
    dt = tt*t0; dt[0] = 0.0                              # t - t0
    cf = fder['c']                                        # f^{(k)}/(k! f)
    K = len(cf) - 1
    fv = np.zeros((NT, n), complex); fv[0] = cf[K]
    for k in range(K - 1, -1, -1):                       # Horner
        fv = mul(fv, dt); fv[0] += cf[k]
    return mul(fv, outer(AJ, BJ))


# ------------------------------------------------------------------ quartic forms
def lin_poly(u):
    """linear form u.z as a degree-1 coefficient array over MON3 (35, n)."""
    out = np.zeros((NM, u.shape[-1]), complex)
    for k in range(3):
        e = [0, 0, 0]; e[k] = 1
        out[IDX3[tuple(e)]] = u[k]
    return out


def side(S, xa, za, xb, x5, W2, mons, sol):
    """T (nstruct, 5[z_b samples], 15 quartic coefs, n) and the v-vector (3, n) for one <T T O> side."""
    V1 = S.Vfun(xa, za, xb, x5)                                                 # (n)
    V2 = S.lin_samples([S.Vfun(xb, S.BASIS[m], x5, xa) for m in range(3)], W2)  # (5, n)
    H12 = S.lin_samples([S.Hfun(xa, za, xb, S.BASIS[m]) for m in range(3)], W2)
    v = np.array([S.Vfun(x5, S.BASIS[k], xa, xb) for k in range(3)])           # (3, n)
    h1 = np.array([S.Hfun(xa, za, x5, S.BASIS[k]) for k in range(3)])
    H25b = np.array([[S.Hfun(xb, S.BASIS[m], x5, S.BASIS[k]) for k in range(3)] for m in range(3)])   # (m, k, n)
    h2 = np.einsum('jm,mkn->jkn', W2, H25b)                                     # (5, 3, n)
    n = v.shape[-1]
    Lv = lin_poly(v); Lh1 = lin_poly(h1)
    Lh2 = [lin_poly(h2[j]) for j in range(5)]
    powv = [None]*5; powh1 = [None]*3
    one = np.zeros((NM, n), complex); one[0] = 1.0
    powv[0] = one; powh1[0] = one
    for k in range(1, 5):
        powv[k] = mul3(powv[k - 1], Lv)
    for k in range(1, 3):
        powh1[k] = mul3(powh1[k - 1], Lh1)
    powh2 = []
    for j in range(5):
        p = [one, Lh2[j], mul3(Lh2[j], Lh2[j])]
        powh2.append(p)
    ns = sol.shape[1]
    T = np.zeros((ns, 5, NM, n), complex)
    cache = {}
    for i, (a1, a2, a5, b12, b15, b25) in enumerate(mons):
        if np.all(np.abs(sol[i]) < 1e-14):
            continue
        assert b15 + b25 <= 4 and a5 + b15 + b25 == len(sol) * 0 + a5 + b15 + b25
        for j in range(5):
            key = (4 - b15 - b25, b15, b25, j)
            if key not in cache:
                cache[key] = mul3(mul3(powv[4 - b15 - b25], powh1[b15]), powh2[j][b25])
            A = (V1**a1)*(V2[j]**a2)*(H12[j]**b12)
            T[:, j] += sol[i][:, None, None]*(A[None, None, :]*cache[key][None])
    return T[:, :, TOP], v


def pairing(Ta, Tc, v, vp, J, t0, mode):
    """P5[a, c, j, l] (relative to the scalar exp(lognorm)) and lognorm = log{[(J-4)!/J!]^2 c_J (g g')^{J/2} f(t0)}."""
    K = 8
    d = WV.wave_derivs(J, t0, mode, K)
    f0 = d[0]
    cf = [d[k]/(math.factorial(k)*f0) for k in range(K + 1)]
    G = g_tps(v, vp, dict(J=J, t0=t0, c=cf))
    G44 = G.reshape(NM, NM, -1)[np.ix_(TOP, TOP)]                              # (15, 15, n)
    fa = FACT3[TOP]
    G44 = G44*(fa[:, None, None]*fa[None, :, None])                             # d^alpha d^beta G
    P5 = np.einsum('ajxn,clyn,xyn->acjln', Ta, Tc, G44)
    a0 = gmet(v, v); b0 = gmet(vp, vp)
    lcJ = J*np.log(2.0) + 2*gammaln(J + 1) - gammaln(2*J + 1)
    lnorm = 2*(gammaln(J - 3) - gammaln(J + 1)) + lcJ + (J/2)*np.log(a0*b0) + np.log(f0)
    return P5, lnorm


CH = int(os.environ.get('PW_GCH', '96'))


def kind_parts(ev, zl):
    """(lamF, Gs, D) for kind ev.kind (point-wise), and the geometry."""
    S, Cr, A6, P9 = PW.S, PW.Cr, PW.A6, PW.P9
    lamF, G0, parts = Cr.build(zl, ev.nu, ev.Delta, 0, n0=1.0, parts=True)
    Pg = A6.geom(zl, 1.0)
    x0, x1, x2, x3, x4, e, z0 = (Pg[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'z0'))
    W2 = S.null_circle(2)
    if ev.kind == 'EE':
        Gs, _ = W0._gs(PW, parts['g'], x0, x1, x2, x3, x4, e, z0, ev.J, A6, S)
        Eb = W0.e_factor(PW, S, x0, z0, x2, x4, W2)
        Dm = S.bm(Eb, Eb)
    elif ev.kind == 'C3H':
        g = parts['g']; msq = A6.msq
        m12, m34 = msq(x1, x2), msq(x3, x4)
        m15, m25, m35, m45 = msq(x1, e), msq(x2, e), msq(x3, e), msq(x4, e)
        v042 = S.Vfun(x0, z0, x4, x2)
        J = ev.J
        Gs = (g['meas'] + g['xsq'] + g['FP'] + g['dual'] + g['cel'] - 4*PW.log(v042)
              - ((10 - J)/2)*PW.log(m12) - (J/2)*(PW.log(m15) + PW.log(m25))
              - ((7 - J)/2)*PW.log(m34) - ((3 + J)/2)*(PW.log(m35) + PW.log(m45)))
        H20 = S.lin_samples([S.Hfun(x2, S.BASIS[m], x0, z0) for m in range(3)], W2)
        H40 = S.lin_samples([S.Hfun(x4, S.BASIS[m], x0, z0) for m in range(3)], W2)
        H24b = np.array([[S.Hfun(x2, S.BASIS[m], x4, S.BASIS[k]) for k in range(3)] for m in range(3)])
        H24 = np.einsum('im,jk,mkn->ijn', W2, W2, H24b)
        inv = 1.0/(S.Xfun(x2, x0)*S.Xfun(x4, x0))
        q = (2*H20[:, None]*H40[None, :] + (v042*v042)[None, None]*H24)*inv
        Dm = q*P9.odd_structure(x0, z0, x2, x4, W2, 'H')
    else:
        raise ValueError(ev.kind)
    return lamF, Gs, Dm, dict(x0=x0, x1=x1, x2=x2, x3=x3, x4=x4, e=e, z0=z0, W2=W2)


def grav(ev, xi, logs_only=False):
    """(Pg, pref (na, nc, n)) of the mode ev.mode (0: P_J, +-1: waves) at the points xi (n3 Evaluator convention)."""
    import pw_core as C
    S = PW.S
    zl, lj = ev._z(xi)
    lamF, Gs, Dm, geo = kind_parts(ev, zl)
    zw, zz = C.zw_geom(zl)
    ts = zw*PW.powr(zz, -0.5)                                                   # continued (last log of the sequence)
    Pg = lamF + Gs + lj
    if logs_only:
        return Pg, None
    x0, x1, x2, x3, x4, e, W2 = (geo[k] for k in ('x0', 'x1', 'x2', 'x3', 'x4', 'e', 'W2'))
    Ta, v = side(S, x1, PW.Cr.sub(x1, x0), x2, e, W2, ev.mons, ev.sol)
    Tc, vp = side(S, x3, PW.Cr.sub(x3, x0), x4, e, W2, ev.mons, ev.solt)
    tV = gmet(v, vp)/np.sqrt(gmet(v, v)*gmet(vp, vp))
    sig = np.sign(np.real(tV/ts))                                               # tV = sig * ts (sig = +-1 per point, checked)
    mode_eff = ev.mode*sig if ev.mode != 0 else 0*sig
    J = ev.J
    out_P = np.zeros((Ta.shape[0], Tc.shape[0], 5, 5, len(ts)), complex)
    lnorm = np.zeros(len(ts), complex)
    for md in (np.unique(mode_eff) if ev.mode != 0 else [0]):
        sel = (mode_eff == md) if ev.mode != 0 else np.ones(len(ts), bool)
        t0 = sig[sel]*ts[sel]                                                   # = tV up to the branch of the sqrt
        ii = np.where(sel)[0]
        for c0 in range(0, len(ii), CH):
            jj = ii[c0:c0 + CH]
            P5, ln = pairing(Ta[..., jj], Tc[..., jj], v[:, jj], vp[:, jj], J, sig[jj]*ts[jj], int(md))
            out_P[..., jj] = P5; lnorm[jj] = ln
    k2 = S.kappa(2)
    Ph = S.fourier(S.fourier(out_P, 2, 2), 3, 2)[:, :, ::-1, ::-1]             # P5hat[a, c, -n2, -n4]
    Dh = S.fourier(S.fourier(Dm, 0, 2), 1, 2)
    pref = np.einsum('ij,i,j,acijn->acn', np.ones((5, 5)), 1/k2, 1/k2, Ph*Dh[None, None])
    return Pg + lnorm, pref
