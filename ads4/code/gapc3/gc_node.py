"""
gc_node (GAP_C3.md): certified polydisc data of one crystal node (J, nu, Delta) for the Laplace remainder lemma of GAP_C3 §2.

For the Morse coordinates v of the Delta-INDEPENDENT phase Phi_xi = lamF/Delta (z = z* + L_Phi v, L_Phi = sqrt(Delta) L, L the
l1_crystal Morse matrix of lamF) and a polydisc P_rho = {|v_i| <= rho}, returns enclosures of
  c0    = g(z*)                      amplitude at the crystal (centre half, ball arithmetic)
  grad  = max |dPhi/dv_i (z*)|       (must enclose 0: exact crystal)
  E     = || d^2Phi(z*) + I ||_max   (Morse normalisation defect of the float L)
  S3    = sum_{|alpha|=3} sup_{P_rho} |Phi_alpha|        (Taylor coefficients of Phi at every point of P_rho)
  Mg    = sup_{P_rho} |g|            (min of the direct enclosure and the mean-value form |g(z*)| + rho sum_i sup|dg/dv_i|)
  bad   = branch / division failures (the polydisc is rejected if non-empty)
Amplitude g = exp(G) x prefactor, exactly as the l1_crystal / l9_points / w0_ee_common builders define it.
Kinds: 'S' (scalar k = 2 functional, l1_crystal with Q_J as prefactor), 'EE' (w0 even/odd EE, k = 2), 'C3H' (l9 / p3, k = 3),
       'C2' (l9 q^2, k = 2).  Spinning kinds return the (na x nc) prefactor matrix; Mg, c0 per entry.
"""
import os, sys, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
import numpy as np
import gc_ball as B
from flint import acb, arb

import l1_tps as TT                      # real float engine (Morse matrix only)
import l1_crystal as CrR


def morse(nu, Delta):
    """float Morse matrix L of lamF (as l1_crystal.laplace: Cholesky + 2 refinements)."""
    zs, a = CrR.saddle(nu, Delta)
    zv = [TT.var(i, zs[i]) for i in range(TT.NV)]
    lamF = CrR.build(zv, nu, Delta, 0)[0]
    _, HF = TT.grad_hess(lamF)
    L = np.linalg.inv(np.linalg.cholesky(-HF + 0j) if False else TT.csym_cholesky(-HF)).T
    for _ in range(2):
        lamF2 = CrR.build(TT.linsub(zs, L), nu, Delta, 0)[0]
        _, H2 = TT.grad_hess(lamF2)
        L = L @ np.linalg.inv(TT.csym_cholesky(-H2)).T
    return zs, a, L


def saddle_ball(nu, Delta):
    """exact crystal z* = (0, -ia, ia, ia, -ia, 0, 0), a = (Delta - sqrt(Delta^2 - nu^2))/nu, in Arb."""
    n, D = arb(nu), arb(Delta)
    a = (D - (D*D - n*n).sqrt())/n
    z = [acb(0), acb(0, -a), acb(0, a), acb(0, a), acb(0, -a), acb(0), acb(0)]
    return z, a


def _abs_up(z):
    return float(abs(z).upper())


def _abs_lo(z):
    return float(abs(z).lower())


def build(kind, J, nu, Delta, rho, sol_override=None):
    """Taylor models (on P_rho, Morse variables of Phi) of lamF, G and the prefactor matrix (na, nc, N)."""
    B.BAD.clear(); B.WMAX[0] = 0.0; B.RHO[0] = rho
    zs, a, L = morse(nu, Delta)
    LP = L*np.sqrt(Delta)                                  # Morse matrix of Phi = lamF/Delta
    zc, aa = saddle_ball(nu, Delta)
    zv = B.linsub_ball(zc, LP)
    if kind == 'S':
        lamF, G, Q = B.Cr.build(zv, nu, Delta, J, pref=True)
        prefs = B.obj(Q)[None, None, :]
    elif J % 2 == 0:
        mons, sol, _ = B.S.conserved(Delta, J)
        _, solt, _ = B.S.conserved(3.0 - Delta, J)
        if sol_override is not None:
            sol, solt = sol_override
        if kind == 'EE':
            bEE = B.W0.make_build_even_EE(B, B.S, B.Cr, B.A6)
            lamF, G, prefs = bEE(zv, nu, Delta, J, mons, sol, solt)
        else:
            lamF, G, prefs = B.P9.build_spin(zv, nu, Delta, J, mons, sol, solt, kind)
    else:
        P3 = B.load_p3()
        if kind == 'EE':
            P3.build = B.W0.make_build_odd_EE(P3)
        out = P3.build(zv, nu, Delta, J, kind, 1.0)
        lamF, G, prefs = out[0], out[1], out[2]
    return dict(lamF=lamF, G=G, prefs=B.obj(prefs), LP=LP, a=float(aa.mid()))


def node(kind, J, nu, Delta, rho, sol_override=None):
    t0 = time.time()
    R = build(kind, J, nu, Delta, rho, sol_override)
    NM = B.NM
    Phi = R['lamF']/acb(Delta)
    grad = max(_abs_up(Phi[1 + i]) for i in range(7))
    E = 0.0
    for k in np.where(B.DEGS == 2)[0]:
        e = B.EXP[k]
        tgt = -0.5 if e.max() == 2 else 0.0
        E = max(E, _abs_up(Phi[k] - tgt))
    # f3 = Phi - Phi(0) - (linear) - (quadratic): sup on P_rho <= sum_{k>=3} B_k + |rem|  (Schwarz: |f3(v)| <= M3 (|v|_inf/rho)^3)
    tot, Bk = B.bound(Phi)
    Bk = Bk[0] if Bk.ndim > 1 else Bk
    M3 = sum(float(Bk[k].upper()) for k in range(3, B.DEG + 1)) + _abs_up(Phi[NM])
    eG = B.exp(R['G'])
    prefs = R['prefs']
    na, nc = prefs.shape[:2]
    c0 = np.zeros((na, nc), complex); c0lo = np.zeros((na, nc)); Mg = np.zeros((na, nc)); remg = np.zeros((na, nc))
    for i in range(na):
        for j in range(nc):
            g = B.mul(eG, prefs[i, j])
            c0[i, j] = complex(g[0].mid()); c0lo[i, j] = _abs_lo(g[0])
            Mg[i, j] = float(B.bound(g)[0][0].upper()) if np.ndim(B.bound(g)[0]) else float(B.bound(g)[0].upper())
            remg[i, j] = _abs_up(g[NM])
    return dict(kind=kind, J=J, nu=nu, Delta=Delta, xi=nu/Delta, rho=rho, a=R['a'], grad=grad, E=E, M3=M3, S3=M3/rho**3,
                remPhi=_abs_up(Phi[NM]), c0=c0, c0lo=c0lo, Mg=Mg, remg=remg, bad=list(B.BAD), qmax=B.WMAX[0], L=R['LP'],
                t=time.time() - t0, deg=B.DEG)


if __name__ == '__main__':
    kind, J, nu, D, rho = sys.argv[1], int(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])
    r = node(kind, J, nu, D, rho)
    for k in ('xi', 'a', 'grad', 'E', 'M3', 'S3', 'remPhi', 'c0', 'c0lo', 'Mg', 'remg', 'bad', 'qmax', 't'):
        print(k, r[k])
