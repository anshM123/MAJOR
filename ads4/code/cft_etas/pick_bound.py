"""
Theorem C (Euclidean-to-Kubo stability) and Proposition C' (sharpness): numerics.

F analytic in H_g = {Re s > -g}, |F| <= M on H_g, |F(n)| <= eps for n = 1..N  ==>
   |F(0)| <= eps*Lam_N + (M + eps*Lam_N)*beta_N,
   beta_N = prod_{n<=N} n/(n+2g) = Gamma(N+1)Gamma(1+2g)/Gamma(N+1+2g),
   Lam_N  = sum_k prod_{j!=k} (k+j+2g)/|k-j|.
Sharpness: F_sh(s) = M * B_N(s) * ((a-g)/(a+s))^k, B_N = prod (s-n)/(s+n+2g), a = N, k = N, N = ceil(log2(M/eps)).
"""
import mpmath as mp
mp.mp.dps = 60
a0 = mp.log(3 + 2*mp.sqrt(2))

def beta(N, g):
    return mp.gamma(N + 1)*mp.gamma(1 + 2*g)/mp.gamma(N + 1 + 2*g)

def Lam(N, g):
    tot = 0
    for k in range(1, N + 1):
        p = mp.mpf(1)
        for j in range(1, N + 1):
            if j != k:
                p *= (k + j + 2*g)/mp.mpf(abs(k - j))
        tot += p
    return tot

def upper(eps, M, g, Nmax=200):
    best = (mp.inf, None)
    for N in range(1, Nmax):
        L = Lam(N, g)
        b = eps*L + (M + eps*L)*beta(N, g)
        if b < best[0]:
            best = (b, N)
        if eps*L > 10*best[0]:
            break
    return best

def sharp(eps, M, g):
    N = int(mp.ceil(mp.log(M/eps)/mp.log(2)))
    a = mp.mpf(N); k = N
    def F(s):
        B = mp.mpf(1)
        for n in range(1, N + 1):
            B *= (s - n)/(s + n + 2*g)
        return M*B*((a - g)/(a + s))**k
    # checks: Matsubara values beyond N, boundary sup
    mats = max(abs(F(n)) for n in range(N + 1, N + 400))
    bnd = max(abs(F(-g + 1j*t)) for t in [0.01*i for i in range(0, 3000, 7)])
    return F(0), mats, bnd, N

if __name__ == '__main__':
    print('a0 = log(3+2 sqrt2) =', mp.nstr(a0, 8), ';  (1/N) log Lam_N(g=1) at N=120:', mp.nstr(mp.log(Lam(120, 1))/120, 8))
    print('%-6s %-6s %-10s | %-12s %-4s | %-12s %-12s %-12s | %-8s' % ('gamma', 'M', 'eps', 'upper', 'N*', 'sharp F(0)', 'max|F(n>N)|', 'sup|F| bdry', 'ratio'))
    for g in (0.5, 1.0, 1.5):
        for (M, eps) in ((3, 1e-12), (3, 1e-16), (3, 1e-24), (3, 1e-40)):
            ub, Nst = upper(mp.mpf(eps), mp.mpf(M), mp.mpf(g))
            F0, mats, bnd, N = sharp(mp.mpf(eps), mp.mpf(M), mp.mpf(g))
            print('%-6.2f %-6g %-10.0e | %-12s %-4d | %-12s %-12s %-12s | %-8s' % (g, M, eps, mp.nstr(ub, 5), Nst,
                  mp.nstr(F0, 5), mp.nstr(mats, 3), mp.nstr(bnd, 5), mp.nstr(ub/F0, 4)))
    print()
    print('rate check: upper*(log(M/eps))^(2g)/M  and  sharp*(log(M/eps))^(2g)/M  (g=1, M=1):')
    for e in (1e-10, 1e-20, 1e-40, 1e-80, 1e-160):
        ub, Nst = upper(mp.mpf(e), 1, 1)
        F0, _, _, _ = sharp(mp.mpf(e), 1, 1)
        Lg = mp.log(1/mp.mpf(e))
        print('eps=%-8.0e  upper*L^2 = %-10s  sharp*L^2 = %-10s' % (e, mp.nstr(ub*Lg**2, 5), mp.nstr(F0*Lg**2, 5)))
