"""lemU_cert/lu_hybrid_dv.py [DV]: check of the 'global coefficients -> c-jets' route (coordinator/author sliver_next, IDEAL weights).
A_k = [z^k] G_A from the author's triangular moment recurrence (midpoint arb), then Y_j(c) = d^3 (cd)^j sum_k C(k,j) A_k c^{k-j},
compared with the per-c integrated jets (lu_core.integrated_jets_mid) and the Lemma U test value.
Usage: lu_hybrid_dv.py L K prec Np"""
import sys, os, math, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from fractions import Fraction as Fr
from flint import arb, ctx
import lu_core as C
L = float(sys.argv[1]); K = int(sys.argv[2]); prec = int(sys.argv[3]); Np = int(sys.argv[4])
ctx.prec = prec
LAM2 = 997002
def coeffs(p, psi):
    I = [C.A(sum((a/Fr(2*j + i + 1) for i, a in enumerate(psi)), Fr(0))) for j in range(K + 2)]
    prev = [arb(0)]*(K + 2); out = []
    for k in range(K + 1):
        out.append(I[0])
        A_ = 8*k*k + 8*p*k + 4*p - 1; B_ = 4*(k + p - 1)**2; D_ = 4*(k + 1)**2
        nxt = [((A_*I[j] - 4*LAM2*I[j + 1] - B_*prev[j])/D_).mid() for j in range(K - k)]
        prev, I = I, nxt + [arb(0)]*2
    return out
t = time.time()
Ak = [2*v for v in coeffs(3, C.PSI2)]; Dk = [LAM2*v for v in coeffs(4, C.PSI3)]
print('global coefficients', time.time() - t, 's; max log10|A_k| =', max(float(abs(v).log()) for v in Ak if v != 0)/math.log(10), flush=True)
c = arb((1 - arb(-L).exp()).mid()); d = 1 - c
J = 3003
t = time.time()
def jets_from(A):
    # Y_j = d^3 (cd)^j sum_{k>=j} C(k,j) A_k c^{k-j}; via repeated synthetic division (Taylor shift of the polynomial sum A_k z^k to z = c)
    a = list(A)
    n = len(a)
    for i in range(n):                       # Horner-type shift: after the loops a[j] = sum_k C(k,j) A_k c^{k-j}
        for k in range(n - 2, i - 1, -1):
            a[k] = (a[k] + c*a[k + 1]).mid()
    out = []; f = d**3
    for j in range(J + 1):
        out.append(f*a[j]); f = f*c*d
    return out
YA_h = jets_from(Ak); YD_h = jets_from(Dk)
print('Taylor shift', time.time() - t, 's', flush=True)
t = time.time()
Y = C.integrated_jets_mid(c, J, Np)
print('per-c quadrature jets', time.time() - t, 's', flush=True)
for j in (0, 1, 10, 100, 500, 1000, 2000, 3000):
    print(j, 'hybrid', YA_h[j].str(12), ' quadrature', Y['A'][j].str(12), '| D:', YD_h[j].str(8), Y['D'][j].str(8))
