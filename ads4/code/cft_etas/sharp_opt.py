"""Optimised sharpness family for Theorem C: F(s) = M B_N(s) ((a-g)/(a+s))^k ; maximise |F(0)| s.t. |F(n)| <= eps for all n>N.
(|F| <= M on Re s >= -g holds for every N, a > g, k >= 0 by the maximum principle; checked on the boundary line.)"""
import mpmath as mp
mp.mp.dps = 50

def logB(N, g, s):
    return sum(mp.log(abs((s - n)/(s + n + 2*g))) for n in range(1, N + 1))

def best(eps, M, g):
    L = mp.log(M/eps)
    out = (mp.mpf(0), None)
    for N in range(1, int(2*L) + 2):
        lb0 = mp.log(M) + logB(N, g, 0)
        for a in [mp.mpf(x) for x in (0.5*N + g + 0.5, N + g, 2*N + g, 4*N + g, 8*N + g, 16*N + g)]:
            # smallest k with max_{n>N} log|F(n)| <= log eps ; check n up to large range (the tail decreases)
            def worst(k):
                m = -mp.inf
                for n in list(range(N + 1, N + 60)) + [int(N*1.5**j) + 1 for j in range(2, 40)]:
                    v = mp.log(M) + logB(N, g, n) + k*mp.log((a - g)/(a + n))
                    m = max(m, v)
                return m
            lo, hi = 0, 1
            while worst(hi) > mp.log(eps):
                hi *= 2
                if hi > 10**6: break
            while hi - lo > 1:
                mid = (lo + hi)//2
                if worst(mid) > mp.log(eps): lo = mid
                else: hi = mid
            k = hi
            val = mp.exp(lb0 + k*mp.log((a - g)/a))
            if val > out[0]:
                out = (val, (N, float(a), k))
    return out

if __name__ == '__main__':
    for g in (0.5, 1.0, 1.5):
        for e in (1e-12, 1e-24, 1e-48, 1e-96):
            v, par = best(mp.mpf(e), mp.mpf(1), mp.mpf(g))
            Lg = mp.log(1/mp.mpf(e))
            print('g=%.1f eps=%-7.0e  sharp |F(0)| = %-11s (N,a,k)=%-18s  |F(0)|*L^(2g) = %s' % (g, e, mp.nstr(v, 5), par, mp.nstr(v*Lg**(2*g), 5)))
