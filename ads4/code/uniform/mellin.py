"""
Exact action of the CHMRS functional C_{k,nu} on a heavy scalar (J = 0) block, d = 3, identical external scalars Df,
from the Mellin representation of CHMRS section 3.5:
   C_{k,nu} = 2^{d-2-2k} pi^{3/2} Gamma((d-1)/2) a_Df / gamma_{2Df+k-1}(nu)^2 * int dT/(4 pi i) Gamma(T/2)^2 gamma_{2Df+k-1-T}(nu) Bhat_{k,T}
   Bhat_{k,T}[G_{Delta,0}] = (-1)^{k/2} sum_n Q_n / [(Df - tau/2 - n)_{k/2} ((tau - T)/2 + n)_{k/2}] [1/(2Df - tau - 2n) + 1/(T - tau - 2n)]
with the J = 0 Mack coefficients  Q_n = Q_0 (1 + Delta/2 - Df)_n^2 / (n! (Delta - d/2 + 1)_n),
Q_0 = -2 Gamma(Delta)/(Gamma(Df - Delta/2)^2 Gamma(Delta/2)^4)  (fixed by the collinear block; derivation in UNIFORM_REMAINDER.md),
and heavy density  C_{k,nu}[Delta,0] = C_{k,nu}[G]/(2 sin^2(pi(tau-2Df)/2) q_{Delta,0}).
The overall normalisation and sign are NOT assumed correct a priori: v2_mellin.py checks the result against the CHMRS
order-8 series (validated against exact blocks in v1) at small nu/m, where both are accurate.
Even k only (the Mellin B-hat sum rules are defined for even k).
"""
import mpmath as mp

def heavy_C(k, nu, Delta, Df, dps=30, nmax=None, tol=None):
    with mp.workdps(dps):
        d = 3
        Delta, Df, nu = mp.mpf(Delta), mp.mpf(Df), mp.mpf(nu)
        tau = Delta
        h = mp.mpf(d)/2
        kk = k//2
        # Mack coefficients (log-scale to avoid overflow), sign of Q0 from 1/Gamma(Df - tau/2)^2 (> 0 unless pole)
        logQ0 = mp.log(2) + mp.loggamma(Delta) - 2*mp.log(abs(mp.gamma(Df - Delta/2))) - 4*mp.loggamma(Delta/2)
        sQ0 = -1
        if nmax is None:
            nmax = int(Delta/2 + 40*mp.sqrt(Delta) + 400)
        # precompute n-dependent pieces
        A = 1 + Delta/2 - Df
        terms = []
        logr = mp.mpf(0)
        for n in range(nmax + 1):
            if n > 0:
                logr += 2*mp.log(abs(A + n - 1)) - mp.log(n) - mp.log(Delta - h + n)
                # sign of (A)_n^2 is +; (Delta-h+1)_n > 0
            c1 = mp.rf(Df - tau/2 - n, kk)
            terms.append((n, logQ0 + logr, c1))
        a = 2*Df + k - 1
        def gam(aa):
            return mp.gamma((aa - mp.mpf(1)/2 + 1j*nu)/2)*mp.gamma((aa - mp.mpf(1)/2 - 1j*nu)/2)
        # normalise the n-sum by its largest log-term to keep numbers O(1)
        Lmax = max(t[1] for t in terms)
        def Bhat(T):
            s = 0
            for (n, lq, c1) in terms:
                w = mp.exp(lq - Lmax)
                if w < mp.mpf(10)**(-dps):
                    continue
                s += w/(c1*mp.rf((tau - T)/2 + n, kk))*(1/(2*Df - tau - 2*n) + 1/(T - tau - 2*n))
            return (-1)**kk*sQ0*s
        c = mp.mpf(1)
        f = lambda t: mp.gamma((c + 1j*t)/2)**2*gam(a - c - 1j*t)*Bhat(c + 1j*t)
        # dT/(4 pi i) with T = c + i t: dT = i dt -> dt/(4 pi)
        I = mp.quad(f, [-mp.inf, -2*nu - 20, -nu, 0, nu, 2*nu + 20, mp.inf])/(4*mp.pi)
        pref = mp.mpf(2)**(d - 2 - 2*k)*mp.pi**(mp.mpf(3)/2)*mp.gamma(mp.mpf(d - 1)/2)
        aDf = 4*mp.gamma(Df)**2*mp.gamma(Df - mp.mpf(1)/2)**2
        CG = pref*aDf/gam(a)**2*I                      # times exp(Lmax)
        # heavy density: divide by 2 sin^2 and q
        import skel
        lq = skel.log_q(Delta, 0, Df)
        s2 = 2*mp.sin(mp.pi*(tau - 2*Df)/2)**2
        val = CG*mp.exp(Lmax - lq)/s2
        return val
