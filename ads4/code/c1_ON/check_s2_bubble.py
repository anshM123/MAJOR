# Validate the curvature expansion Pi_{n,J} = Pi0(u) - (K/3) u Pi0''(u) + O(K^2), u = J(J+1)K, K = 1/r^2,
# against the exact 3j-sum bubble on S^1 x S^2 (beta = 1).
import numpy as np, time
from s2bubble import bubble, m_of_r, M0
import flat
from mpmath import mpf
for r in [4, 8, 16]:
    m = m_of_r(r)
    print(f'r={r}  m(r)={m:.12f}  (m^2-m0^2) r^2 = {(m*m-M0*M0)*r*r:.8f}')
    for n in [0, 1, 3]:
        for J in [0, 1, 4, 12]:
            t = time.time()
            ex, err = bubble(n, J, r, L0=max(4000, 200*r))
            K = 1.0/r**2; u = J*(J+1)*K
            P0 = float(flat.Pi0(n, mpf(u))); P2 = float(flat.Pi0(n, mpf(u), 2))
            corr = -(K/3)*u*P2
            print(f'  n={n} J={J:2d} exact={ex:.13e} (+-{err:.1e})  Pi0={P0:.13e}  (ex-Pi0)/Pi0={(ex-P0)/P0: .4e}'
                  f'  pred={corr/P0: .4e}  resid*r^4/Pi0={(ex-P0-corr)/P0*r**4: .4e}')
