"""
tp/tcE_common.py: TAIL_CERT's multiplier bound K_osc (Theorem P) sharpened by Theorem E' (ANALYTIC_SLIVER §5.2) where valid.
Both are certified upper bounds of |S^_w(sigma + i tau)| (w-normalisation, TAIL), for all tail states with L >= L0, a >= 998, b in {0,4}:
  * Theorem E' (b = 0, d = w): closed form, uniform in spin and in L, for tau <= 0.95 a; increasing-in-kappa factors evaluated with
    tau as an Arb interval (enclosure over the cell), worst case a = 998 (every factor is nonincreasing in a).
  * b = 4: helicity identity S^(9)_{n,N}(beta) = sum_i C(4,i)(-1)^i S^(1)_{n+4,N+4}(beta-i) / U^(1)_4(n+4,N+4)  ([DV] 1e-80, ANALYTIC_SLIVER §5.4),
    |w^beta| <= |d^beta| (w <= d for b = 4), d/d0 = (D+20)/D, U_4 >= (1-d0)^4 (1-1e-6); i >= 1 terms carry d0^i <= e^{-i L0}, with
    |S^(1)(beta - i)| <= 2^i max(E'(sigma), E'(sigma+1)) d0^{-0} (beta-recurrence growth <= 2 per step for |beta| <= 1500, TAIL §2.3).
The tailcert modules are imported read-only; the Theorem-P table is cached in analytic/cache.
"""
import os, sys, math, json
from flint import arb, acb, ctx
ADS4 = '/home/user/MAJOR/ads4'
TC = os.path.join(ADS4, 'code', 'tailcert')
sys.path.insert(0, TC)
import tc_power as P
CACHE = os.path.join(ADS4, 'code', 'analytic', 'cache')
_tab = os.path.join(CACHE, f'ehat_d_table_L{P.ETAB_L0:g}.json')
def _Ehat_table():
    if P._ETAB is None:
        if os.path.exists(_tab):
            P._ETAB = json.load(open(_tab))
        else:
            P._ETAB = [P.Ehat_d(P.ETAB_U0 + P.ETAB_DU*i, 0.92, P.ETAB_DELTA) for i in range(P.ETAB_N)]
            json.dump(P._ETAB, open(_tab, 'w'))
    return P._ETAB
P.Ehat_table = _Ehat_table
A_MIN = 998
KAPPA_MAX = 0.95


def thmE(sigma, t0, t1, a=A_MIN):
    """Theorem E' (b = 0) upper bound valid for every tau in [t0, t1], sigma >= 0, every a >= a_given; None if t1 > 0.95 a.
    Sharpened form (ANALYTIC_SLIVER §5.2): 2 c0 e^{sigma/((a+1)sqrt(1-k^2))} X^sigma exp(q asin k) (1-k^2)^{-q^2/2}, q = |alpha|/tau.
    Every factor of the closed form is nondecreasing in kappa = tau/a, in q = |alpha|/tau and in |alpha|, so it is evaluated (in Arb, upward)
    at kappa_max = t1/a, q_max = sqrt(sigma^2 + t1^2)/max(t0, tiny), |alpha|_max = sqrt(sigma^2 + t1^2)."""
    if t1 > KAPPA_MAX*a:
        return None
    ctx.prec = 128
    s = arb(sigma); aa = arb(a)
    k = arb(t1)/aa
    al = (s*s + arb(t1)**2).sqrt()
    q = al/arb(max(t0, 1e-3)) if t0 > 0 else None
    if q is None:                       # t0 = 0: use |alpha|/tau * asin(kappa) <= |alpha| asin(k)/(k a) form: (|alpha|/tau) asin(tau/a) <= |alpha|/a * (asin k / k)
        q1 = al/aa*(k.asin()/k); q2 = None
    sq = (1 - k*k).sqrt()
    X = 2/(1 + sq)
    if q is not None:
        e1 = q*(k.asin() + q*(-(1 - k*k).log()/2))
    else:
        # (|alpha|/tau)^2 (1/2) log(1/(1-k^2)) <= (|alpha|/tau)^2 k^2/(2(1-k^2)) = (|alpha|/a)^2/(2 sq^2)
        e1 = q1 + (al/aa)**2/(2*sq*sq)
    e2 = s/((aa + 1)*sq)
    c0 = 1 + 2*al/(aa*aa)
    return float((2*c0*X**s*(e1 + e2).exp()).upper())


def thmE_b(sigma, t0, t1, L0):
    """bound for both channels (max over b = 0, 4) of |S^_w|, using the helicity identity for b = 4."""
    e0 = thmE(sigma, t0, t1)
    if e0 is None:
        return None
    e1 = thmE(sigma + 1, t0, t1)
    if e1 is None:
        return None
    d0 = math.exp(-L0)*1.0001
    m = max(e0, e1)
    corr = sum(math.comb(4, i)*d0**i*2**i*m for i in range(1, 5))
    b4 = (1 + 1e-5)*(e0 + corr)/((1 - d0)**4*(1 - 1e-6))
    return max(e0, b4)


def K_osc(sigma, t1, t2, L0):
    """min(Theorem P, Theorem E') -- same signature/return convention as tc_power.K_osc."""
    b, res = P.K_osc(sigma, t1, t2, L0)
    e = thmE_b(sigma, t1, t2, L0)
    if e is not None and e < b:
        return e, [('E', e)]
    return b, res


def logK_osc(sigma, t1, t2, L0):
    e = thmE_b(sigma, t1, t2, L0)
    l = P.logK_osc(sigma, t1, t2, L0)
    return min(l, math.log(e)) if e is not None else l


if __name__ == '__main__':
    L0 = float(sys.argv[1]) if len(sys.argv) > 1 else 16
    for sg in (0.25, 1.25, 3.5, 4.5):
        for (a_, b_) in ((0.25, 1), (10, 11), (45, 46), (300, 301), (600, 601), (900, 901), (940, 941)):
            bp, _ = P.K_osc(sg, a_, b_, L0)
            print(f'sigma={sg} tau in [{a_},{b_}]: Theorem P {bp:.4g}   Theorem E\' (both b) {thmE_b(sg, a_, b_, L0):.4g}', flush=True)
