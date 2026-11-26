"""tailcert/tc_author_dv.py [DV]: spot checks of the author's Lemma R' inputs (ANGULAR_ROUTE §82) at small parameters:
 (i) S_b(alpha) = Gamma(N+1+alpha)/(Gamma(N+1)Gamma(alpha)) int_0^1 q^{alpha-1}(1-q)^a K(q) dq, K(q) = B(a,n+s)^{-1} int t^{a-1}(1-t)^{n+s-1}(1-qt)^n dt
 (ii) |K(q)| <= 2|1+q|^{-a} on sample points of the closed first-quadrant unit disc (as stated in §37 for n >= 1024 a^2; here sampled)."""
import mpmath as mp
mp.mp.dps = 30
def S_dir(beta, n, a, s):
    N = n + a; U = mp.mpf(1); tot = 0; poch = mp.mpf(1)
    for k in range(n + 1):
        tot += U*poch; poch *= (beta + k)/(k + 1); U *= mp.mpf((n - k)*(n + s + k))/((N - k)*(N + s + k))
    return tot
def Kq(q, n, a, s):
    return mp.quad(lambda t: t**(a - 1)*(1 - t)**(n + s - 1)*(1 - q*t)**n, [0, 1])/mp.beta(a, n + s)
for (n, a, s) in ((12, 4, 1), (15, 6, 9)):
    N = n + a
    for beta in (mp.mpf(2.5), mp.mpc(0.25, 3)):
        rep = mp.gamma(N + 1 + beta)/(mp.gamma(N + 1)*mp.gamma(beta))*mp.quad(lambda q: q**(beta - 1)*(1 - q)**a*Kq(q, n, a, s), [0, 0.5, 1])
        print(f'(i) n={n} a={a} s={s} beta={beta}: direct {mp.nstr(S_dir(beta, n, a, s), 12)}  author-rep {mp.nstr(rep, 12)}')
n, a, s = 1024*16, 4, 1
worst = 0
for r in (0.3, 0.7, 1.0):
    for th in (0, 0.4, 0.8, 1.2, 1.5707):
        q = r*mp.expj(th)
        v = abs(Kq(q, n, a, s))*abs(1 + q)**a
        worst = max(worst, v)
print(f'(ii) n = 1024 a^2 = {n}, a = {a}: max |K(q)| |1+q|^a over samples = {mp.nstr(worst, 6)} (author: <= 2)')
