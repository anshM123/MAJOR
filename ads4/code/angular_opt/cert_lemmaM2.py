"""Fast exact certificate for Lemma-M bounds (integer arithmetic).
q(z) = Q(z)/2^(60+E), Q(z) = sum_j A_j t_j(z) 2^(E-e_j) + shift, with t_j the integer Gram polynomials on {0..Mb},
nu_j = 2^e_j >= max|t_j|.  Exact check: Q(0) = 2^(60+E); Q(z) <= 2^(60+E); Q(z) >= -theta_a 2^(60+E) on z = 0..Mb.
Usage: python3 cert_lemmaM2.py Mb n ga   (theta_a = theta(ga))."""
import numpy as np, sys, json, math, time
from fractions import Fraction as Fr
from scipy.optimize import linprog
from cert_lemmaM import gram_rows, WB, theta

def run(Mb, n, ga, eps=1e-6):
    t0 = time.time()
    es = []; G = np.zeros((n+1, Mb+1))
    for j, t in gram_rows(Mb, n):
        m = max(abs(x) for x in t); e = m.bit_length()
        es.append(e)
        G[j] = np.array([math.ldexp(float(x >> max(0, e-60)) if e > 60 else float(x), -min(e, 60)) for x in t])
    th = theta(Fr(ga)); thf = float(th)
    obj = -(G[:, 2] - 2*G[:, 1])
    Aub = np.vstack([G[:, 1:].T, -G[:, 1:].T])
    bub = np.concatenate([(1-eps)*np.ones(Mb), (thf-eps)*np.ones(Mb)])
    res = linprog(obj, A_ub=Aub, b_ub=bub, A_eq=G[:, 0:1].T, b_eq=[1.0], bounds=[(None, None)]*(n+1), method='highs')
    assert res.status == 0, res.message
    A = [int(round(x*2**60)) for x in res.x]
    E = max(es)
    Q = [0]*(Mb+1)
    for j, t in gram_rows(Mb, n):
        f = A[j] << (E - es[j])
        for z in range(Mb+1): Q[z] += f*t[z]
    one = 1 << (60+E)
    shift = one - Q[0]
    Q = [v + shift for v in Q]
    mx = max(Q); mn = min(Q)
    ok = (Q[0] == one) and (mx <= one) and (Fr(mn, one) >= -th)
    S = Fr(Q[2] - 2*Q[1], one) + 1
    out = dict(Mb=Mb, n=n, ga=str(ga), theta_a=str(th), ok=bool(ok), S=str(S), Sfloat=float(S),
               maxq=float(Fr(mx, one)), minq=float(Fr(mn, one)), E=E, es=es, A=[str(a) for a in A], shift=str(shift),
               secs=time.time()-t0)
    return out

if __name__ == '__main__':
    Mb = int(sys.argv[1]); n = int(sys.argv[2]); ga = sys.argv[3]
    out = run(Mb, n, ga)
    print({k: out[k] for k in ['Mb', 'n', 'ga', 'ok', 'Sfloat', 'maxq', 'minq', 'secs']}, flush=True)
    json.dump(out, open(f'cache/certI_M{Mb}_n{n}_ga{ga}.json', 'w'))
