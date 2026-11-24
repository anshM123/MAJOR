"""
gc_cert (GAP_C3.md §3): certified polydisc data for a crystal node from the log-factored integrand (gc_scalar, Taylor models of
gc_ball on P_rho in the Morse variables v of the Delta-independent phase Phi_xi).

Output per (J, xi, rho):  q (max |f_k/f_k* - 1| over P_rho: < 1 certifies analyticity of every Log / power), E (Morse defect),
grad (|dPhi(z*)|), M3 = sup_P |f3| (f3 = Phi - Phi* - quadratic), S3 = M3/rho^3, Mg = sup_P |g/g*_env| with g*_env = g(z*)/P_J(cos th*)
(the envelope normalisation of Gamma), c0n = |P_J(cos th*)| (= |g(z*)|/g*_env).
Usage: python3 gc_cert.py J nu Delta rho [deg]
"""
import os, sys, time, json
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import gc_guard  # noqa: F401
if len(sys.argv) > 5 and __name__ == '__main__':
    os.environ['GC_DEG'] = sys.argv[5]
import numpy as np
import gc_ball as B
import gc_scalar as GS
from flint import acb, arb
import gc_node as N


def saddle_xi(xlo, xhi):
    """crystal for every xi in [xlo, xhi]: a = xi/(1 + sqrt(1 - xi^2)) (ball), z* = (0, -ia, ia, ia, -ia, 0, 0)."""
    x = arb(xlo).union(arb(xhi))
    a = x/(1 + (1 - x*x).sqrt())
    return [acb(0), acb(0, -a), acb(0, a), acb(0, a), acb(0, -a), acb(0), acb(0)], a, acb(x)


def cert(J, nu, Delta, rho, xint=None):
    """xint = (xlo, xhi): certify for every xi in the interval (Morse matrix at the midpoint; E then carries the xi-variation)."""
    t0 = time.time()
    B.BAD.clear(); B.WMAX[0] = 0.0; B.RHO[0] = rho
    zs, a, L = N.morse(nu, Delta)
    LP = L*np.sqrt(Delta)
    if xint is None:                       # exact xi = nu/Delta (the same number as the crystal's), so grad Phi(z*) = 0 exactly
        zc, aa = N.saddle_ball(nu, Delta); xib = acb(arb(nu)/arb(Delta))
    else:
        zc, aa, xib = saddle_xi(*xint)
    zv = B.linsub_ball(zc, LP)
    K = GS.TMB(B)
    r = GS.build(zv, nu, Delta, J, K, acb(aa), xi=xib)
    NM = B.NM
    dPhi = r['dPhi']
    up = lambda z: float(abs(z).upper())
    grad = max(up(dPhi[1 + i]) for i in range(7))
    E = 0.0
    for k in np.where(B.DEGS == 2)[0]:
        tgt = -0.5 if B.EXP[k].max() == 2 else 0.0
        E = max(E, up(dPhi[k] - tgt))
    _, Bk = B.bound(dPhi)
    Bk = np.asarray(Bk).ravel()
    M3 = sum(float(Bk[k].upper()) for k in range(3, B.DEG + 1)) + up(dPhi[NM])
    g = r['grel']
    if J > 0:
        env = abs(r['ZZWW'][0])**(arb(J)/2)                    # |ZZWW(z*)|^{J/2}:  h(z*) = P_J(cos th*) x phase
        g = B.mul(g, r['Q'])/acb(env)
    tot, Bh = B.bound(g)
    Bh = np.asarray(Bh).ravel()
    Mg = float(np.asarray(tot).ravel()[0].upper())
    out = dict(J=J, nu=nu, Delta=Delta, xi=nu/Delta, xint=xint, rho=rho, deg=B.DEG, a=float(aa.mid()), q=B.WMAX[0], bad=len(B.BAD),
               grad=grad, E=E, M3=M3, S3=M3/rho**3, remPhi=up(dPhi[NM]), Mg=Mg, h0=abs(complex(g[0].mid())),
               h0lo=float(abs(g[0]).lower()),
               Bphi=[float(Bk[k].upper()) for k in range(B.DEG + 1)], Bh=[float(Bh[k].upper()) for k in range(B.DEG + 1)],
               remh=up(g[NM]), t=time.time() - t0)
    return out, r


if __name__ == '__main__':
    J, nu, D, rho = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), float(sys.argv[4])
    o, _ = cert(J, nu, D, rho)
    print(json.dumps(o))
