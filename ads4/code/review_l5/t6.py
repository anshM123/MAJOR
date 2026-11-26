# finite-mass Mehler kernel: kq(t) = L[(1-z)^{-5/2} (1 - z/z_t)^{-1/2}],  z_t = tanh^2(t/2); infinite twist: (1-c)^{-5/2}(1-c/z_t)^{-1/2} (t>beta_c)
import sys
from an_common import *
ctx.prec = 6000
J,a,b = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
st = state(J,a,b); U = ladder(st); c = st['c']; n = st['n']
bc = float((2*c.sqrt().atanh()).mid())
print('beta_c', bc)
A = [arb(1)]
for i in range(n): A.append(A[-1]*(arb(5)/2+i)/(i+1))
def kq(t):
    zt = (arb(t)/2).tanh()**2
    Bc = [arb(1)]
    for j in range(n): Bc.append(Bc[-1]*(arb(1)/2+j)/((j+1)*zt))
    # coefficient of z^k of product: conv
    tot = arb(0)
    # sum_k U_k sum_{i<=k} A_i B_{k-i} = sum_j B_j sum_{i<=n-j} U_{i+j} A_i
    for j in range(n+1):
        s = arb(0)
        for i in range(n-j+1):
            s += U[i+j]*A[i]
        tot += Bc[j]*s
    return tot
import math
for t in [float(x) for x in sys.argv[4].split(',')]:
    v = kq(t)
    zt = math.tanh(t/2)**2; cf = float(c.mid())
    vi = (1-cf)**-2.5*(1-cf/zt)**-0.5 if zt > cf else float('nan')
    print(f't={t:.3f} k_fm={float(v.mid()):+.6e}  k_inf={vi:.6e}', flush=True)
