# Method A ingredient: exact 3j-sum bubble vs curvature-expanded bubble, summed over J in each Matsubara sector:
#   diff_n(r) = sum_{J<=Jc} (2J+1) [log Pi^exact_{n,J} - log P_n(u_J)],  P_n = Pi0 - (K/3) u Pi0''.
# Since everything else in log Z_1 is identical between the exact and the expanded treatment, the method-B constant is
# exact iff  sum_n diff_n(r) -> 0  as r -> infinity (expected O(1/r^2)), with convergence in Jc.
import sys, numpy as np, time
from s2bubble import bubble, m_of_r
import fastflat as ff
r = float(sys.argv[1]); ns = [int(a) for a in sys.argv[2].split(',')]; Jc = int(float(sys.argv[3])*r)
K = 1/r**2
for n in ns:
    t = time.time()
    acc = 0.0; rows = []
    for J in range(0, Jc+1):
        ex, err = bubble(n, J, r, L0=max(3000, int(150*r + 30*J)))
        u = J*(J+1)*K
        P0, P1, P2 = ff.Pi0(n, np.array([u]))
        d = (2*J+1)*(np.log(ex) - np.log(P0[0] - (K/3)*u*P2[0]))
        acc += d
        if J in (0, 1, 2, 5, 10, 20, 40, 80, 160) or J == Jc:
            rows.append((J, d, acc, err/ex))
    print(f'r={r} n={n} Jc={Jc} diff={acc:.6e}  diff*r^2={acc*r*r:.5f}  ({time.time()-t:.0f}s)', flush=True)
    for row in rows:
        print('    J=%d term=%.3e cum=%.6e relerr=%.1e' % row, flush=True)
