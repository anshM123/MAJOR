"""DV checks of the Loc''' weight bounds (K0, P_W tail, Chernoff Beta tail) against exact values (mpmath)."""
import sys, math
sys.path.insert(0, '../analytic2')
import mpmath as mp
import loc5
from common2 import n_of_L
mp.mp.dps = 50
for L in (3.6, 3.0, 2.0, 1.0):
    for a in (998, 1500, 4000):
        for b in (0, 4):
            s = 2*b + 1
            n = n_of_L(L, a, b); N = n + a
            d = 1 - mp.mpf(n*(n + s))/(N*(N + s)); Lx = float(-mp.log(d))
            t_lo, t_hi = loc5.t_range(Lx - 1e-9, Lx + 1e-9)
            t = mp.mpf(a)/N; rho = 1 - t
            assert t_lo <= t <= t_hi, (t, t_lo, t_hi)
            g = mp.mpf(a + 1)/2
            IW = mp.hyp2f1(g, g, 1, rho**2)
            lC = mp.loggamma(N + 1) - mp.loggamma(a + 1) - mp.loggamma(n + 1)
            K0 = mp.exp(-n*mp.log(rho) + mp.log(IW) - lC)
            out = [f'L={Lx:.3f} a={a} b={b}: K0={mp.nstr(K0, 6)} (bound {loc5.K0_bound():.4f})']
            for th0 in (2*t_hi, 5*t_hi, 0.3):
                W = lambda th: (1 - 2*rho*mp.cos(th) + rho**2)**(-(a + 1)/mp.mpf(2))
                num = 2*mp.quad(W, [th0, th0 + 4*t, th0 + 20*t, mp.pi]); den = 2*mp.pi*IW
                ex = mp.log(num/den)
                out.append(f'th0={th0:.4f}: log P_W ex {float(ex):.2f} bd {loc5.lfrac(th0, t_lo, t_hi, a):.2f}')
            for u0 in (2*t_hi, 5*t_hi, 0.5):
                M = N + s - 1; uu = mp.mpf(u0)
                lp = [mp.loggamma(M + 1) - mp.loggamma(k + 1) - mp.loggamma(M - k + 1) + k*mp.log(uu) + (M - k)*mp.log(1 - uu) for k in range(a)]
                m_ = max(lp); ex = m_ + mp.log(sum(mp.exp(x - m_) for x in lp))
                out.append(f'u0={u0:.4f}: log P_B ex {float(ex):.2f} bd {loc5.lchern(u0, t_hi, a):.2f}')
            print('  '.join(out), flush=True)
