"""
energy_dv.py [DV]: the b = 0 spin recurrence at fixed N and the quadratic-energy bound (Theorem E, ANALYTIC_SLIVER §6).
  S_{n+1} = S_{n-1} + eps_n alpha S_n,  eps_n = 2(2n+1)/((N-n)(N+n+1)),  S_{-1} = S_0 = 1,
  S_n(alpha) = L_{n,N}((1-z)^{-alpha}),  d_n = 1 - n(n+1)/(N(N+1)).
Bound: |d^alpha S_n|^2 <= 2 (1-y_f)^{-2} exp( sum_k eps_k [sigma |y_k| + eps_k sigma^2] / (1-|y_k|) ),  y_k = eps_k tau/2, y_f = y_{n-1}.
"""
import sys, math
from flint import arb, acb, ctx
ctx.prec = 8000
def direct(n, N, al):
    s = 1; C0 = n*(n+s); D = N*(N+s)
    U = arb(1); t = acb(1); S = acb(1)
    for k in range(1, n+1):
        U = U*(C0 - (k-1)*(k-1+s))/(D - (k-1)*(k-1+s))
        t = t*(al + k - 1)/k
        S += U*t
    return S
def run(nf, a, sig, tau, check=False):
    N = arb(nf) + arb(a); al = acb(sig, tau)
    Sm, S = acb(1), acb(1)
    lsum = 0.0; ys = []
    for n in range(0, nf):
        eps = 2*(2*n + 1)/((N - n)*(N + n + 1))
        y = float((eps*tau/2).mid()); e = float(eps.mid())
        lsum += e*(sig*abs(y) + e*sig*sig)/(1 - abs(y))
        ys.append(y)
        Sm, S = S, Sm + eps*al*S
    d = 1 - arb(nf*(nf + 1))/(N*(N + 1))
    Sh = abs(d**al*S)
    yf = ys[-1]
    bound = math.sqrt(2*math.exp(lsum))/(1 - yf)
    out = f'n={nf} a={a} alpha={sig}+{tau}i: |S^| = {float(Sh.mid()):.6g}   bound = {bound:.4g}  (y_f = {yf:.4f}, sum = {lsum:.3f})'
    if check:
        Sd = direct(nf, N, al)
        out += f'   [direct check rel err {float((abs(Sd - S)/abs(S)).mid()):.1e}]'
    return out
if __name__ == '__main__':
    print(run(30, 20.5, 2.5, 7.3, check=True))
    print(run(200, 998.5, 1.25, 300, check=True))
    for nf, a in [(1001, 998), (5000, 998), (100000, 998), (2000, 3000), (20000, 1500)]:
        for sig in (0.25, 2.5):
            for tau in (10, 300, 700, 900, 990):
                print(run(nf, a, sig, tau))
