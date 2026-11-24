"""
gm_cert (GAP_MORSE.md §3): certified Taylor-model data of one crystal node on the HERMITIAN ball B_rho = {v in C^7: |v|_2 <= rho}
of the Morse coordinates v of the Delta-independent phase Phi_xi (z = z* + L_Phi v, L_Phi = sqrt(Delta) L, L = l1_crystal Morse matrix;
the same v as GAP_C3).  Taylor models: gm_ball (Bombieri degree bounds), integrand: gm_scalar (= gapc3/gc_scalar, log-factored).

Output (all upper bounds unless stated; every number is an Arb ball endpoint):
  q, bad        max |f_k/f_k* - 1| over B_rho (< 1 certifies analyticity of every log / power) and the number of failures
  qPhi, badPhi  the same for the compositions needed by the phase alone (Phi certified on B_rho if badPhi = 0)
  grad          max |dPhi/dv_i(0)| (an enclosure; contains 0: exact crystal, hypothesis (b))
  EF            Frobenius norm of E' = -Hess Phi(0) - I
  trE           |tr E'|
  bPhi[k]       Bombieri norm [Phi_k]_B of the degree-k homogeneous Taylor part of Phi - Phi*, k = 3..D (radius-1 normalisation)
  remPhi        sup_{B_rho} |Phi - Phi* - P_D|  (Taylor-model remainder; vanishes to order D+1)
  bTr[k]        Bombieri norms of the homogeneous parts (degree k-2) of p_tr = -sum_{k=3}^D Lap(Phi_k)/k   (index k)
  bh[k], remh   the same for h = g / g_env (g_env = g(z*)/P_J(cos th*), the envelope normalisation of GAP_C3 §1)
  h0, h0lo      |h(0)|
Usage: GC_DEG=6 python3 gm_cert.py J xi rho
"""
import os, sys, time, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gm_guard  # noqa: F401
import numpy as np
import gm_ball as B
import gm_scalar as GS
from flint import acb, arb
sys.path.insert(0, os.path.join(HERE, '..', 'gapc3'))
import gc_node as N                    # read-only: float Morse matrix (morse) and exact crystal ball (saddle_ball)

DREF = 1002.0                          # any Delta: at fixed xi every quantity here is Delta-independent (GAP_UNIFORM Prop. S)


def up(z):
    return float(abs(z).upper()) if isinstance(z, acb) else float(z.upper())


def lap_coeffs(c):
    """coefficients (over B.MONS) of Lap(p) = sum_i d_i^2 p for a polynomial with coefficients c (object array, len NM)."""
    idx = {m: i for i, m in enumerate(B.MONS)}
    out = B.zeros(B.NM)
    for a, m in enumerate(B.MONS):
        if sum(m) < 2:
            continue
        for i in range(7):
            if m[i] >= 2:
                mm = list(m); mm[i] -= 2
                out[idx[tuple(mm)]] = out[idx[tuple(mm)]] + c[a]*(m[i]*(m[i] - 1))
    return out


def cert(J, xi, rho, Delta=DREF):
    t0 = time.time()
    B.BAD.clear(); B.WMAX[0] = 0.0; B.RHO[0] = rho
    nu = xi*Delta
    zs, a, L = N.morse(nu, Delta)
    LP = L*np.sqrt(Delta)
    zc, aa = N.saddle_ball(nu, Delta); xib = acb(arb(nu)/arb(Delta))
    zv = B.linsub_ball(zc, LP)
    K = GS.TMB(B)
    snap = {}
    def mark():
        snap['bad'] = len(B.BAD); snap['q'] = B.WMAX[0]
    r = GS.build(zv, nu, Delta, J, K, acb(aa), xi=xib, mark=mark)
    NM = B.NM
    dPhi = r['dPhi']
    grad = max(up(dPhi[1 + i]) for i in range(7))
    # E' = -Hess - I ; Hess_ii = 2 c_{2e_i}, Hess_ij = c_{e_i + e_j}
    E = [[None]*7 for _ in range(7)]
    for k in np.where(B.DEGS == 2)[0]:
        e = B.EXP[k]
        nz = [i for i in range(7) if e[i]]
        if len(nz) == 1:
            i = nz[0]; E[i][i] = -2*dPhi[k] - 1
        else:
            i, j = nz; E[i][j] = E[j][i] = -dPhi[k]
    EF = sum(abs(E[i][j])**2 for i in range(7) for j in range(7)).sqrt()
    trE = abs(sum(E[i][i] for i in range(7)))
    bPhi = B.bombieri(dPhi[:NM][None])[0]
    pc = dPhi[:NM].copy()
    tr = B.zeros(NM)
    for k in range(3, B.DEG + 1):
        ck = pc.copy(); ck[B.DEGS != k] = B.ZERO
        tr = tr - lap_coeffs(ck)/acb(k)
    bTr0 = B.bombieri(tr[None])[0]                    # index = degree of p_tr part (k - 2)
    bTr = [0.0, 0.0, 0.0] + [up(bTr0[k - 2]) for k in range(3, B.DEG + 1)]   # index k: the part from Phi_k
    g = r['grel']
    if J > 0:                       # h = Q_J g/(g* |ZZWW*|^{J/2}) up to a unimodular constant (gm_scalar, centred Legendre)
        g = r['hJ']
    bh = B.bombieri(g[:NM][None])[0]
    out = dict(J=J, xi=xi, nu=nu, Delta=Delta, rho=rho, deg=B.DEG, norm=B.NORM, a=float(aa.mid()), q=B.WMAX[0], bad=len(B.BAD), qPhi=snap['q'], badPhi=snap['bad'],
               grad=grad, EF=up(EF), trE=up(trE), bPhi=[up(bPhi[k]) for k in range(B.DEG + 1)], remPhi=up(dPhi[NM]),
               bTr=bTr, bh=[up(bh[k]) for k in range(B.DEG + 1)], remh=up(g[NM]), h0=abs(complex(g[0].mid())),
               h0lo=float(abs(g[0]).lower()), t=time.time() - t0)
    return out


if __name__ == '__main__':
    J, xi, rho = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    o = cert(J, xi, rho)
    print(json.dumps(o))
