# Vacuum limit: the T = 0, m = 0 bubble on R x S^2 (3j sum) vs the conformal Delta=1 formula
#   Pi_vac(w,J) = (r/16) |Gamma((J+1+i w r)/2) / Gamma((J+2+i w r)/2)|^2,
# and its curvature expansion Pi_vac = 1/(8p) - (K/3) u d^2_u (1/(8p)) + O(K^2) at fixed w != 0.
import numpy as np
from s2bubble import log3j2
from mpmath import mp, mpf, mpc, gamma, sqrt, pi, log
mp.dps = 20

def vac_exact(w, J, r):
    z1 = mpc(J + 1, w * r) / 2; z2 = mpc(J + 2, w * r) / 2
    return float(r / 16 * abs(gamma(z1) / gamma(z2)) ** 2)

def vac_3j(w, J, r, L0=20000):
    def part(Ls):
        Lmax = max(Ls); l1 = np.arange(Lmax + 1); tot = np.zeros(Lmax + 1)
        for d in range(-J, J + 1, 2):
            l2 = l1 + d
            ok = (l2 >= 0) & (l2 >= np.abs(l1 - J)) & (l2 <= l1 + J)
            L1, L2 = l1[ok], l2[ok]
            E1 = (L1 + 0.5) / r; E2 = (L2 + 0.5) / r
            tot[ok] += (2*L1+1)*(2*L2+1)*np.exp(log3j2(L1, L2, J)) * 2*(E1+E2)/(w*w+(E1+E2)**2)/(4*E1*E2)
        cs = np.cumsum(tot) / (4*np.pi*r*r)
        return np.array([cs[x] for x in Ls])
    Ls = [int(L0 * 2 ** (i / 2)) for i in range(7)]
    S = part(Ls); x = 1.0/np.array(Ls, float)
    A = np.vander(x[-6:], 6, increasing=True)
    return np.linalg.solve(A, S[-6:])[0]

for (w, J, r) in [(0.0, 0, 1.0), (0.0, 3, 1.0), (1.3, 2, 1.0), (5.0, 7, 2.0)]:
    print('w,J,r', w, J, r, '3j-sum', vac_3j(w, J, r), 'Gamma formula', vac_exact(w, J, r))
print('curvature expansion of the reference at w = 2 pi:')
w = 2*np.pi
for r in [8, 16, 32]:
    for J in [3, 20]:
        K = 1.0/r**2; u = J*(J+1)*K; p = np.sqrt(u + w*w)
        ex = vac_exact(w, J, r); P0 = 1/(8*p); corr = -(K/3)*u*(3/4)*(u+w*w)**-2.5/8
        print(f'  r={r} J={J} (ex-P0)/corr={(ex-P0)/corr:.8f}  resid*r^4/P0={(ex-P0-corr)/P0*r**4:.6e}')
