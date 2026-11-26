"""
thmE.py [certified scalar evaluation]: closed form of Theorem E' (ANALYTIC_SLIVER §5.2), b = 0:
  |d^alpha S_n(alpha)| <= 2 c0 X^sigma exp((|alpha|/tau)[asin k + (|alpha|/tau)((1-k^2)^{-1/2} - 1)]) e^{sigma/((a+1)sqrt(1-k^2))} e^{4|alpha|/(a^2(1-k^2)^{3/2})},
  X = 2/(1 + sqrt(1-k^2)), k = tau/a; c0 <= 1 + 2|alpha|/N^2 (|w_0| = |alpha|/(N(N+1)) <= 1e-5 on the sliver).
thmE_bound(sigma, tau, a, N=None) -> arb upper bound (all n, all real N >= n + a).  Random adversarial test against exact multipliers.
"""
import sys, random
from flint import arb, acb, ctx
def thmE_bound(sigma, tau, a):
    ctx.prec = 128
    s = arb(sigma); t = arb(tau); a = arb(a)
    k = t/a
    assert k.upper() < 1
    al = (s*s + t*t).sqrt()
    q = al/t if tau > 0 else arb(1)
    sq = (1 - k*k).sqrt()
    X = 2/(1 + sq)
    e1 = q*(k.asin() + q*(-(1 - k*k).log()/2))
    e2 = s/((a + 1)*sq)
    c0 = 1 + 2*al/(a*a)        # N >= a
    return 2*c0*X**s*(e1 + e2).exp()
if __name__ == '__main__':
    sys.path.insert(0, '.')
    random.seed(7)
    worst = 0.0
    for trial in range(int(sys.argv[1]) if len(sys.argv) > 1 else 40):
        n = random.choice([1001, 1500, 2500, 5000, 12000])
        a = random.uniform(998, 3600)
        sig = random.choice([0.25, 1.25, 2.5, 3.5, 4.5])
        tau = random.uniform(0.5, 0.95*a if a < 1100 else min(1100, 0.95*a))
        ctx.prec = 8000
        N = arb(n) + arb(a); al = acb(sig, tau)
        Sm, S = acb(1), acb(1)
        for m in range(n):
            eps = 2*(2*m + 1)/((N - m)*(N + m + 1)); Sm, S = S, Sm + eps*al*S
        d = 1 - arb(n*(n + 1))/(N*(N + 1))
        ex = float(abs(d**al*S).mid())
        b = float(thmE_bound(sig, tau, a).upper())
        worst = max(worst, ex/b)
        print(f'n={n} a={a:.1f} alpha={sig}+{tau:.1f}i: exact {ex:.4g}  bound {b:.4g}  ratio {ex/b:.3f}', flush=True)
    print('max exact/bound =', worst)
