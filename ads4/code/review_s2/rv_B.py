"""rv_B.py [exact containment + DV]: Lemma B as used in K_cell.
(1) the state box (bchain.state_box) contains the exact g = 1/(d^2 N(N+s)) of every tested integer state (both channels, cell edges,
    every a-range of rt2_kernel.a_ranges incl. the geometric ones and [2e5, oo));
(2) the K_cell chain xs (chain_bounds_multi started from the certificate's PP_box M0, M1) against exact |S^(5/2 - j + i nu)|,
    j = -1..jmax-1, at integer corner states of the box (fastS: dual-Hahn recurrence in n, complex128 + rescaling)."""
import sys, math, random
sys.path.insert(0, '/home/user/MAJOR/ads4/code/review_s2/ws/code/analytic2')
import common2
import rt2_kernel as RK
import bchain as BC
from fastS import S_vec
from fractions import Fraction as Fr
from flint import arb, acb, ctx
ctx.prec = 300
random.seed(11)
La, Lb = float(sys.argv[1]), float(sys.argv[2])
d_hi, d_lo = math.exp(-La), math.exp(-Lb)


def dval(n, a, s):
    N = n + a
    return a*(2*n + a + s)/(N*(N + s))


def n_edge(a, s, dt, side):
    lo, hi = 1, 10**15
    while hi - lo > 1:
        m = (lo + hi)//2
        if dval(m, a, s) > dt: lo = m
        else: hi = m
    return hi if side == 'hi' else lo


# (1) containment of g
bad = 0; cnt = 0
for (a_lo, a_hi) in RK.a_ranges(1032.1):
    _, g, _ = BC.state_box(a_lo, a_hi, d_lo, d_hi)
    for a in {int(math.ceil(a_lo)), int(min(a_hi, 1e12)) if a_hi < math.inf else 10**12}:
        for b in (0, 4):
            s = 2*b + 1
            for (dt, side) in ((d_hi, 'hi'), (d_lo, 'lo')):
                n = n_edge(a, s, dt, side)
                if n < 1001 - b: continue
                d = Fr(a*(2*n + a + s), (n + a)*(n + a + s))
                if not (d_lo <= float(d) <= d_hi): continue
                gx = 1/(d*d*(n + a)*(n + a + s))
                cnt += 1
                if not (float(g.lower()) <= float(gx) <= float(g.upper())):
                    bad += 1; print('g NOT contained', a_lo, a_hi, a, b, n, float(gx), g)
print(f'(1) g containment: {cnt} states, {bad} failures', flush=True)

# (2) chain vs exact
jm = 110 if La < 3.9 else 60
for (T0, alo) in [(1003.0, 998.0), (1031.0, 998.0), (1031.0, 1031.0), (1020.0, 1010.0), (900.0, 998.0), (500.0, 998.0), (1045.0, 1045.0)]:
    T1 = T0 + 1.1; ahi = alo + 1
    M0 = RK.PP_box(2.5, T0, T1, alo, ahi, d_lo, d_hi); M1 = RK.PP_box(3.5, T0, T1, alo, ahi, d_lo, d_hi)
    xs = BC.chain_bounds_multi(2.5, T0, T1, M0, M1, jm, d_lo, d_hi, alo, ahi)
    worst = -1e9; at = None
    for a in (int(alo), int(ahi)):
        for b in (0, 4):
            s = 2*b + 1
            for n in (n_edge(a, s, d_hi, 'hi'), n_edge(a, s, d_lo, 'lo')):
                dd = dval(n, a, s)
                for nu in (T0, T1):
                    S0, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu)])
                    N = arb(n) + arb(a); D = N*(N + s); dA = 1 - arb(n*(n + s))/D
                    for pert in (1.0, 1.0 + 1e-9):     # stability check of the downward beta-recurrence
                        Sd = {0: acb(complex(S0[0]))*arb(lg[0]).exp()*pert, 1: acb(complex(S0[1]))*arb(lg[1]).exp()}
                        al = acb(2.5, nu)
                        for j in range(0, jm):
                            beta = al - j
                            Sd[-(j + 1)] = ((dA*D + beta*(s - 2*beta - 1))*Sd[-j] + beta*(beta + 1)*Sd[-j + 1])/(D - beta*(beta - s))
                        exs = [float((abs(Sd[1 - k]*dA**(al + 1 - k))).log().mid()) for k in range(jm + 1)]
                        if pert == 1.0: ex0 = exs
                        else: stab = max(abs(x - y) for x, y in zip(ex0, exs))
                    for k in range(jm + 1):
                        r = ex0[k] - math.log(xs[k])
                        if r > worst: worst = r; at = (a, b, n, nu, k - 1, round(stab, 9))
    print(f'(2) L[{La},{Lb}] nu[{T0},{T1}] a[{alo},{ahi}]: max_j log(exact/bound) = {worst:.3f} at (a,b,n,nu,j,stab)={at};'
          f' bound j=6: {xs[7]:.3e}', flush=True)
