# [exact/DV] check of the dual-Hahn difference equation in n for S_n(alpha) = 3F2(alpha, -n, n+s; -N, N+s; 1) at real N:
#   alpha S_n = B_n S_{n+1} - (B_n + D_n) S_n + D_n S_{n-1},
#   B_n = (N+n+s)(n+s)(N-n)/((2n+s)(2n+s+1)),  D_n = -n(N+n+s)(N-n)/((2n+s-1)(2n+s)).   (s = 1: B = -D, S_{n+1} - S_{n-1} = 2(2n+1) alpha S_n/((N-n)(N+n+1)))
from flint import arb, acb, ctx
ctx.prec = 300
def S(n, N, s, al):
    C0 = n*(n+s); D = N*(N+s); U = arb(1); t = acb(1); tot = acb(1)
    for k in range(1, n+1):
        U = U*(C0 - (k-1)*(k-1+s))/(D - (k-1)*(k-1+s)); t = t*(al + k - 1)/k; tot += U*t
    return tot
mx = 0
for s in (1, 9):
    for N in (arb('57.3'), arb('1234.75')):
        for al in (acb(2.5, 7.3), acb(-3.1, 40), acb(0.25, 300)):
            for n in (1, 5, 17, 40):
                B = (N+n+s)*(n+s)*(N-n)/((2*n+s)*(2*n+s+1)); Dn = -n*(N+n+s)*(N-n)/((2*n+s-1)*(2*n+s))
                lhs = al*S(n, N, s, al); rhs = B*S(n+1, N, s, al) - (B+Dn)*S(n, N, s, al) + Dn*S(n-1, N, s, al)
                r = float((abs(lhs - rhs)/abs(lhs)).mid()); mx = max(mx, r)
print('max relative defect over 48 cases:', mx)
