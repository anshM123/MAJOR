"""check_B.py: [DV] Lemma B (chain_bounds_multi) vs exact |S^(5/2 - j + i nu)| (n-recurrence for 5/2, 7/2 + exact beta-recurrence,
acb), random states; M0, M1 = exact values."""
import sys, math, random
sys.path.insert(0, '.')
from common2 import n_of_L
from fastS import S_vec
from flint import arb, acb, ctx
import bchain as BC
ctx.prec = 300
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
worst = -1e9
for it in range(int(sys.argv[2]) if len(sys.argv) > 2 else 8):
    a = random.choice([998.0, random.uniform(998, 1100), random.uniform(1100, 5000)]); L = random.uniform(3.75, 5.0)
    nu = random.uniform(300, 1100); b = random.choice([0, 4]); s = 2*b + 1
    n = n_of_L(L, a, b)
    N = arb(n) + arb(a); D = N*(N + s); d = 1 - arb(n*(n + s))/D
    S0, lg = S_vec(n, a, s, [complex(2.5, nu), complex(3.5, nu)])
    Sd = {0: acb(complex(S0[0]))*arb(lg[0]).exp(), 1: acb(complex(S0[1]))*arb(lg[1]).exp()}
    al = acb(2.5, nu)
    for j in range(0, 30):
        beta = al - j
        Sd[-(j + 1)] = ((d*D + beta*(s - 2*beta - 1))*Sd[-j] + beta*(beta + 1)*Sd[-j + 1])/(D - beta*(beta - s))
    ex = [float(abs(Sd[-j]*d**(al - j)).mid()) for j in range(-1, 30)]
    df = float(d.mid())
    bd = BC.chain_bounds_multi(2.5, nu - 0.3, nu + 0.3, ex[1]*1.000001, ex[0]*1.000001, 29, df*0.999, df*1.001, a - 0.5, a + 0.5)
    r = max(math.log(ex[k]/bd[k]) for k in range(2, 31))
    worst = max(worst, r)
    print(f'a={a:.1f} b={b} L={L:.3f} nu={nu:.1f}: max_j log(exact/bound) = {r:.2f};  j=6: exact {ex[7]:.3e} bound {bd[7]:.3e}', flush=True)
print('worst log(exact/bound) =', worst, '(must be <= 0)')
