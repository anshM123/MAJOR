#!/usr/bin/env python3
"""Referee check: near_sum (bottom/left contour majorant of sum_j |c_j| d^j |dressed X_j - 1| + infinite tail)
against the exact value at physical states, nu = 1/4 + i y (left edge) and nu = 0.1 (bottom)."""
import sys
from fractions import Fraction as F
from flint import arb, acb, ctx
sys.path.insert(0, '/home/user/MAJOR/ads4/code/sliver_ext')
from contour_direct import near_sum
ctx.prec = 400
J = 1500

def state(a, b, dt):
    s = 2*b+1; lo, hi = 1, 10**8
    while hi-lo > 1:
        m = (lo+hi)//2; N = m+a
        if F(a*(N+m+s), N*(N+s)) <= dt: hi = m
        else: lo = m
    n = lo; N = n+a; D = N*(N+s)   # d just above dt
    return n, N, D, s

def Xdirect(n, D, s, beta):
    c = arb(n*(n+s))/D; h = arb(1)/D; d = 1-c
    tot = acb(0); term = acb(1); U = arb(1)
    for k in range(n+1):
        tot += term*U; term = term*(beta+k)/(k+1); U = U*(c-h*k*(k+s))/(1-h*k*(k+s))
    return acb(d)**beta*tot

def exact_near(a, b, p, dt, nu):
    n, N, D, s = state(a, b, dt); q = b*(b+1)
    c = arb(n*(n+s))/D; h = arb(1)/D; d = 1-c
    alpha = acb(arb(2*p-1)/2)+acb(0, 1)*nu
    Xm1 = Xdirect(n, D, s, alpha+1); X0 = Xdirect(n, D, s, alpha)
    X = {-1: Xm1, 0: X0}
    for j in range(0, J):  # X_{j+1} = X(alpha-j-1) from homogeneous beta recurrence at beta=alpha-j
        be = alpha-j
        X[j+1] = ((1+h*be*(s-2*be-1)/d)*X[j]+h*be*(be+1)/(d*d)*X[j-1])/(1-h*be*(be-s))
    tau = h/(1+q*h)
    W = arb(1); tot = arb(0)
    a0 = acb(arb(1)/2)-acb(0, 1)*nu; b0 = 1-acb(0, 2)*nu
    cj = acb(1)
    for j in range(J):
        be = alpha-j
        if p == 3:
            A0 = arb(9)/2-b*(b-1)-be*(be+8-4*b)/2; A1 = be*(be+9-4*b)/2
            term = X[j]+tau*(A0*X[j]+A1*X[j-1]/d)-1
        else:
            term = X[j]-1
        tot += abs(cj)*d**j*abs(term)
        cj = cj*(a0+j)*(a0+j)/((b0+j)*(j+1))
    # infinite twist tail sum_{j>=J}|c_j| d^j  (direct to 4000 more terms + crude remainder ignored: exact is smaller)
    tail = arb(0)
    for j in range(J, J+20000):
        tail += abs(cj)*d**j; cj = cj*(a0+j)*(a0+j)/((b0+j)*(j+1))
    return tot+tail, n, float((1-c).mid())

dl, du = F(68779, 1250000), F(46421, 625000)
worst = 0
for p, b, t, y in [(3, 4, '0.25', '0'), (3, 4, '0.25', '5'), (3, 0, '0.25', '10.25'), (3, 4, '0.25', '10.25'), (4, 4, '0.25', '12'), (4, 0, '0.25', '20.25'), (3, 4, '0.1', '0')]:
    tf, yf = F(t), F(y)
    bound = near_sum(dl, du, p, b, tf, tf, yf, yf)
    for a, dt in [(998, du), (998, dl+F(1, 10**6)), (3000, du)]:
        nu = acb(arb(tf.numerator)/tf.denominator, arb(yf.numerator)/yf.denominator)
        ex, n, d = exact_near(a, b, p, dt, nu)
        r = float(ex.upper())/float(bound); worst = max(worst, r)
        print(p, b, t, y, a, n, 'd=%.6f' % d, 'exact %.4g' % float(ex.mid()), 'bound %.4g' % float(bound), 'ratio %.3g' % r, flush=True)
print('worst', worst, 'PASS' if worst < 1 else 'FAIL')
