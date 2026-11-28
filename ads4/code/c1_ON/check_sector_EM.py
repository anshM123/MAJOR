# Check of the per-sector Euler-Maclaurin/curvature expansion for n != 0:
#   S_n(r) = sum_J (2J+1)[log P_n(u_J) - log Pi_vac(w_n,J)]  vs  r^2 A_n + B_n + O(1/r^2)
# with P_n = Pi0_n - (K/3) u Pi0_n'' (validated against exact 3j sums) and the exact Gamma-function reference.
import numpy as np, sys
from scipy.special import loggamma
import fastflat as ff
from c1_modes import mode, tails
x, wgl = np.polynomial.legendre.leggauss(48)

def logvac(w, J, r):
    z1 = (J + 1 + 1j*w*r)/2; z2 = (J + 2 + 1j*w*r)/2
    return np.log(r/16) + 2*np.real(loggamma(z1) - loggamma(z2))

def Sn(n, r, tfac=200):
    w = 2*np.pi*n; K = 1/r**2
    Jc = int(tfac*r)
    J = np.arange(0, Jc+1, dtype=float); u = J*(J+1)*K
    P0, P1, P2 = ff.Pi0(n, u)
    d = np.log(P0 - (K/3)*u*P2) - logvac(w, J, r)
    head = ((2*J+1)*d).sum()
    # tail J > Jc: midpoint E-M, integrand (2J+1) delta(J) with delta from the K-expanded forms of both factors
    a = Jc + 0.5
    def dl(Jx):
        uu = Jx*(Jx+1)*K
        Q0, Q1, Q2 = ff.Pi0(n, uu); p2 = uu + w*w
        th = ff.pi_th(n, uu)[0]
        eps = -(2/np.pi)*np.arctan(2*ff.M0/np.sqrt(p2)) + 8*np.sqrt(p2)*th
        return (2*Jx+1)*(np.log1p(eps) + np.log1p(-(K/3)*uu*Q2/Q0) - np.log1p(-(K/3)*uu*(3/4)/p2**2))
    ed = [a]
    while ed[-1] < 1e5*r: ed.append(ed[-1]*1.6)
    ed = np.array(ed)
    Jq = np.concatenate([0.5*(b-c)*x+0.5*(b+c) for c, b in zip(ed[:-1], ed[1:])])
    W = np.concatenate([0.5*(b-c)*wgl for c, b in zip(ed[:-1], ed[1:])])
    integ = dl(Jq) @ W
    Um = ed[-1]*(ed[-1]+1)*K
    tA, tD = tails(Um, w*w)
    integ += r*r*(tA + K*tD)
    g = dl(np.array([a-1, a+1]))
    return head + integ + (g[1]-g[0])/2/24

for n in [1, 2]:
    d0, I, B, A, _, _ = mode(n)
    for r in [5., 10., 20., 40.]:
        s = Sn(n, r)
        print(n, r, 'S_n - r^2 A_n - B_n =', s - r*r*A - B, ' (x r^2:', (s - r*r*A - B)*r*r, ')', flush=True)
