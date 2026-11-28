"""Exact certificate for the Lemma-M bound (see README, Milestone 2).
Find a polynomial q of degree <= n with q(0)=1 and -theta_a <= q(z) <= 1 on z=0..Mb, maximising S = q(2)-2q(1)+1.
Basis: integer-valued discrete Chebyshev (Gram) polynomials t_j on {0..Mb}, divided by integer scales nu_j.
LP in floats, then the coefficients are rationalised and q is evaluated EXACTLY (Python integers / Fractions)
on the whole grid.  The check is exact: q(0)=1, max q <= 1, min q >= -theta_a."""
import numpy as np, sys, json, math
from fractions import Fraction as Fr
from scipy.optimize import linprog

def gram_rows(Mb, n):
    """yield (j, t_j as list of python ints on z=0..Mb)"""
    N = Mb+1
    t0 = [1]*N; yield 0, t0
    t1 = [2*z - Mb for z in range(N)]; yield 1, t1
    tm, tc = t0, t1
    for j in range(1, n):
        # (j+1) t_{j+1} = (2j+1)(2z-Mb) t_j - j (N^2 - j^2) t_{j-1}
        a = 2*j+1; b = j*(N*N - j*j)
        tn = [(a*(2*z-Mb)*tc[z] - b*tm[z]) // (j+1) for z in range(N)]
        # exactness check of division
        assert all((a*(2*z-Mb)*tc[z] - b*tm[z]) % (j+1) == 0 for z in (0, 1, Mb//2, Mb))
        yield j+1, tn
        tm, tc = tc, tn

def scales(Mb, n):
    nus = []; G = np.zeros((n+1, Mb+1))
    for j, t in gram_rows(Mb, n):
        nu = max(abs(x) for x in t)
        nus.append(nu); G[j] = np.array([x/nu for x in t], dtype=float)  # exact int / int -> float
    return nus, G

def certify(Mb, n, theta_a, eps=1e-6):
    nus, G = scales(Mb, n)
    obj = -(G[:, 2] - 2*G[:, 1])
    Aub = np.vstack([G[:, 1:].T, -G[:, 1:].T])
    th = float(theta_a)
    bub = np.concatenate([(1-eps)*np.ones(Mb), (th-eps)*np.ones(Mb)])
    res = linprog(obj, A_ub=Aub, b_ub=bub, A_eq=G[:, 0:1].T, b_eq=[1.0], bounds=[(None, None)]*(n+1), method='highs')
    assert res.status == 0, res.message
    a = [Fr(x).limit_denominator(2**62) for x in res.x]
    # q(z) = sum_j a_j t_j(z)/nu_j ; enforce q(0)=1 exactly by fixing a_0
    q = [Fr(0)]*(Mb+1)
    t0vals = []
    for j, t in gram_rows(Mb, n):
        cj = a[j]/nus[j]
        t0vals.append(t[0])
        for z in range(Mb+1):
            q[z] += cj*t[z]
    # adjust constant term: t_0 = 1 everywhere
    shift = 1 - q[0]
    q = [v + shift for v in q]
    assert q[0] == 1
    mx = max(q); mn = min(q)
    ok = (mx <= 1) and (mn >= -theta_a)
    S = q[2] - 2*q[1] + 1
    return dict(ok=ok, S=S, max=mx, min=mn, a=a, nus=nus, shift=shift, q=q)

T = lambda D: (D-1)*(D-2)
def WB(g): return Fr(T(g+4)-12, 3) + Fr(3, 2)
def theta(g): return (Fr(T(g+2), 3) + Fr(3, 2))/WB(g)

if __name__ == '__main__':
    Mb = int(sys.argv[1]); n = int(sys.argv[2]); ga = Fr(sys.argv[3])
    theta_a = theta(ga)
    r = certify(Mb, n, theta_a)
    print('Mb', Mb, 'n', n, 'theta_a', float(theta_a), 'ok', r['ok'], 'max', float(r['max']), 'min+theta', float(r['min']+theta_a), 'S', float(r['S']), flush=True)
    with open(f'cache/cert_M{Mb}_n{n}_ga{ga}.json', 'w') as f:
        json.dump(dict(Mb=Mb, n=n, ga=str(ga), theta_a=str(theta_a), S=str(r['S']), ok=bool(r['ok']),
                       a=[str(x) for x in r['a']], nus=[str(x) for x in r['nus']], shift=str(r['shift'])), f)
