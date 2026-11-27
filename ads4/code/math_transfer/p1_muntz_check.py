"""P1: Theorem C's exact-data factor beta_N = prod n/(n+2 gamma) is the Muntz(-Szasz) product.

Mellin picture: F(s) = int_0^1 f(x) x^s x^{2gamma-1} dx is the Hardy space H^2 of {Re s > -gamma} (Paley-Wiener/Mellin
isometry, up to 2*pi), with ||F||^2 = int |f|^2 x^{2gamma-1} dx.  Data F(n) = <f, x^n>, target F(0) = <f, 1>.
The worst |F(0)| over ||f|| <= 1 and F(1..N) = 0 is dist(1, span{x^1..x^N}) in L^2(x^{2gamma-1}dx), which by the
Cauchy (Gram) determinant equals ||1|| * prod_n n/(n+2gamma) = (2gamma)^{-1/2} beta_N.
We compute the Gram-determinant distance in 60-digit arithmetic and compare.
Output: out_p1_muntz.txt
"""
import mpmath as mp

mp.mp.dps = 60


def dist(gamma, N):
    # <x^a, x^b> = int x^{a+b+2gamma-1} dx = 1/(a+b+2gamma)
    idx = list(range(0, N + 1))
    G = mp.matrix(N + 1, N + 1)
    for i, a in enumerate(idx):
        for j, b in enumerate(idx):
            G[i, j] = 1 / mp.mpf(a + b + 2 * gamma)
    Gs = G[1:, 1:]
    return mp.sqrt(mp.det(G) / mp.det(Gs))


if __name__ == "__main__":
    out = []
    for gamma in (mp.mpf("0.25"), mp.mpf("0.5"), mp.mpf(1), mp.mpf("1.5")):
        for N in (1, 2, 5, 10, 20):
            d = dist(gamma, N)
            beta = mp.fprod(mp.mpf(n) / (n + 2 * gamma) for n in range(1, N + 1))
            pred = beta / mp.sqrt(2 * gamma)
            out.append(f"gamma={mp.nstr(gamma, 3):>4} N={N:2d}  Gram distance={mp.nstr(d, 15):>20}  "
                       f"(2g)^-1/2*beta_N={mp.nstr(pred, 15):>20}  rel.diff={mp.nstr(abs(d / pred - 1), 3)}")
    txt = "\n".join(out)
    print(txt)
    open(__file__.replace("p1_muntz_check.py", "out_p1_muntz.txt"), "w").write(txt + "\n")
