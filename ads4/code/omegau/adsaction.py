"""AdS heavy action of Omega^U per heavy state, from the exact scalar skeleton (CHMRS order-1/m^8) in the B' dictionary.

Per state (m = mfl Lam, J, h; u = 1/mfl^2, x = p sqrt(u)), model of the k-th moment [I: B' dictionary, reading B, skeleton transfer]:
  eps_k m^{2k} C^h_k(Lam p) = e_k(x^2) D_h(1-2x^2)  +  2 eps_k sigma^k [P_{-1/2+i Lam p}(cosh beta) - J_0(bhat p)]  +  sigma^k env_k(x) rho_k(x)
  env_k = E_k (even k; for ++ this is literally the exact scalar moment),  (1-x^2)^{-(k+1)/2} (odd k: scalar relative remainder transferred
  to the odd B' envelope; the odd moments are the k-continued scalar ones [I]).
Forward: a u^3 [flat(2JJ - c_h) + (2 m^2 W - 2JJ - 1/2) + R4] + (5a/2) sigma u^3 [2 P0 + R5].
Returned quantities (normal form, per |c^h|^2):
  Fh   hybrid action (flat + leading-Regge curvature, all k and forward);  Ns = stripped component norm, NU = sum_k |part_k| + |FW|
  S    signed skeleton remainder on controlled nodes (incl. forward if controlled);  Es = its error estimate
  Sall signed remainder on all nodes (truncated series value also where 'uncontrolled'; used only for the flat-limit check u3)
  E    absolute bound sum |rho| + err on controlled nodes (per k, per node; no cross-k cancellation);  U = weight of uncontrolled nodes
"""
import numpy as np
from common import KS, EPS, SIGMA, CH, peval, ek

def envk(k, x):
    return (2 - x*x)/(1 - x*x)**(1 + k/2) if k % 2 == 0 else 1/(1 - x*x)**((k + 1)/2)

def action(F, PH, TR, h, st, J, mfl, Lam, pn, pw, sk, flatk, flatFW, strip):
    """flatk: per-k flat parts at (h, u, J) (array of 6); flatFW: flat forward; strip: (B2, B3, FW) stripped flat parts."""
    s = SIGMA[h]
    u = 1/mfl**2
    m, beta = st['m'], st['beta']
    bh = (2*J + 1)/mfl
    I, W, P0 = TR.regge(Lam, beta)
    f = TR.flat(bh)
    curv = np.array([u**(k - 2)*2*EPS[k]*s**k*(I[k] - f[k]) for k in KS])
    JJ = J*(J + 1.0)
    FWh = flatFW + F.a*u**3*(2*m*m*W - 2*JJ - 0.5) + 5*F.a*s*u**3*(P0 - 1)
    parts = flatk + curv
    B2, B3, FWs = strip
    # flat part taken from the stripped functional (identical by the exact density identity, u1; avoids the floating-point
    # cancellation among the per-k B' parts at m -> Lam); the per-k parts enter only the norm NU
    Fh = B2 + B3 + FWs + curv.sum() + (FWh - flatFW)
    NU = np.abs(parts).sum() + abs(FWh)
    Ns = abs(B2 + curv[0]) + abs(B3 + curv[1] + 0) + abs(FWh)      # stripped-style component norm (B2, B3 incl. their curvature)
    x = pn*np.sqrt(u)
    E = U = S = Es = Sall = 0.0
    for i, k in enumerate(KS):
        rho, err, ctrl = sk[k]
        wk = u**(k - 2)*peval(PH[k], pn)*pw*envk(k, x)
        S += np.sum((s**k)*wk*rho*ctrl); Es += np.sum(np.abs(wk)*err*ctrl); Sall += np.sum((s**k)*wk*rho)
        E += np.sum(np.abs(wk)*(np.abs(rho) + err)*ctrl)
        U += np.sum(np.abs(wk)*(~ctrl))
    R4, e4, R5, e5, cF = sk['fw']
    if cF:
        S += F.a*u**3*(R4 + 2.5*s*R5); Es += F.a*u**3*(e4 + 2.5*e5)
        E += F.a*u**3*(abs(R4) + e4 + 2.5*(abs(R5) + e5))
    else:
        U += abs(FWh)
    Sall += F.a*u**3*(R4 + 2.5*s*R5)
    return dict(Fh=Fh, Ns=Ns, NU=NU, S=S, Es=Es, E=E, U=U, Sall=Sall, beta=beta, curv=curv)
