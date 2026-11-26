"""
tailaudit/ta_power_check.py [DV]: independent evaluation of S^(beta) = w^beta S(beta) at tail states and comparison with TAIL_CERT's
certified bounds (tc_power.K_osc = Theorem P, tc_power.NB_bound = Lemma NB).

Method (independent of tc_power's estimates, but uses the exact representation verified in ta_identities.py):
  S^(beta) = w^beta kappa (2a)^{-beta} int_path x^{beta-1} E_V[F_V(x)] dx,   F_V(x) = ((1 - lam_V x/2a)/(1 + x/2a))^a,
  lam_V = 1/V - 1, V ~ Beta(N+1, n+s)  (E_V by 24-point Gauss-Hermite on the Beta's mean/variance; V-outside-window and the
  K > N binomial tail are < e^{-10^6} here and dropped);  path x: i[0, Y1] -> arc |x| = Y1 -> [Y1, Xend] (real), analytic in between.
  'ft' = fixed-twist limit (lam = 1, w(N+s) -> 2a).  Regime A (tau <= 0.52 a): Y1 = 2 tau, double precision is enough (no
  cancellation); regime B (tau > 0.52 a): Y1 = 1.84 a, mpmath at 70 digits (fixed twist only).
Validation: direct ladder sum at n = 3e5, a = 998 (L = 5), against the same code.
Output: cache/ta_power_check.txt
"""
import os, sys, math, cmath
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(ADS4, 'code', 'tailcert'))  # tc_power is imported read-only
sys.dont_write_bytecode = True
out = []
def pr(s):
    print(s, flush=True); out.append(s)

GLX, GLW = np.polynomial.legendre.leggauss(24)
GHX, GHW = np.polynomial.hermite.hermgauss(24)
GHW = GHW/np.sqrt(np.pi)


def state(a, b, L):
    """integer n with w = 1 - C/T closest to e^{-L}; exact integers where a is integer."""
    s = 2*b + 1
    n = int(round(2*a*math.exp(L)))
    for _ in range(60):
        N = mp.mpf(n) + a; T = N*(N + s) + b*(b + 1); w = mp.mpf(a)*(N + n + s)/T
        Lc = -mp.log(w)
        n = int(n*mp.e**(float(L - Lc)))
    N = mp.mpf(n) + a; T = N*(N + s) + b*(b + 1); w = mp.mpf(a)*(N + n + s)/T
    return n, w


def panels(lo, hi, width):
    k = max(1, int(math.ceil((hi - lo)/width)))
    e = np.linspace(lo, hi, k + 1)
    xs = ((e[1:] - e[:-1])[:, None]*(GLX[None, :] + 1)/2 + e[:-1, None]).ravel()
    ws = ((e[1:] - e[:-1])[:, None]*GLW[None, :]/2).ravel()
    return xs, ws


def J_numpy(beta, a, lams, lw, Y1=None):
    """e^{+pi tau/2} * int_path x^{beta-1} E[F] dx  (regime A path, Y1 = 2 tau), lams/lw: lambda nodes & weights."""
    sig, tau = beta.real, beta.imag
    Y1 = max(2*tau, 2.0) if Y1 is None else Y1
    def Fbar(x, off=0):
        z = x/(2*a)
        base = -np.log1p(z)
        acc = 0
        for lam, wt in zip(lams, lw):
            acc = acc + wt*np.exp(a*(np.log1p(-lam*z) + base) + off)
        return acc
    tot = 0j
    # leg x = i y, y in [eps, Y1] in u = ln y; [0, eps] analytically (F = 1 + O(eps))
    eps = 1e-14
    tot += cmath.exp((beta)*math.log(eps) + 1j*math.pi/2*sig)/beta          # int_0^eps y^{beta-1} i^{beta} dy * e^{pi tau/2}
    us, ws = panels(math.log(eps), math.log(Y1), min(0.25, 1.0/(tau + 1)))
    y = np.exp(us)
    x = 1j*y
    # x^{beta-1} dx = e^{(beta-1)(u + i pi/2)} i e^u du ; times e^{pi tau/2}
    lg = beta*us + 1j*math.pi/2*(sig - 1) + 1j*math.pi/2       # (beta-1)(i pi/2) + pi tau/2 = i pi/2 (sig-1) ; i = e^{i pi/2}
    tot += np.sum(ws*np.exp(lg)*Fbar(x))
    # arc x = Y1 e^{i phi}, phi: pi/2 -> 0
    ph, wp = panels(0.0, math.pi/2, min(0.05, 1.0/(tau + 1)))
    x = Y1*np.exp(1j*ph)
    lg = (beta - 1)*(math.log(Y1) + 1j*ph) + math.pi*tau/2 + np.log(1j*Y1) + 1j*ph
    tot -= np.sum(wp*Fbar(x, lg))
    # real axis [Y1, Xend]
    Xend = min(1.84*a, Y1 + 120 + math.pi*tau/2)
    if Xend > Y1:
        xr, wr = panels(Y1, Xend, min(0.5, 2.0/(tau + 1)*Y1))
        lg = (beta - 1)*np.log(xr) + math.pi*tau/2
        tot += np.sum(wr*Fbar(xr.astype(complex), lg))
    return tot


def Shat(beta, a, b, L=None, ft=False, Y1=None):
    beta = complex(beta)
    sig, tau = beta.real, beta.imag
    if ft:
        lams, lw = [1.0], [1.0]
        lpre = mp.mpc(0)
    else:
        s = 2*b + 1
        with mp.workdps(90):
            n, w = state(a, b, L)
            N = mp.mpf(n) + a; M = N + n + s
            mu = (N + 1)/(M + 1); var = mu*(1 - mu)/(M + 2)
            V = [mu + mp.sqrt(2*var)*float(z) for z in GHX]
            lams = [float(1/v - 1) for v in V]; lw = list(GHW)
            mb = mp.mpc(beta)
            lpre = mb*mp.log(w) + mp.loggamma(N + s) - mp.loggamma(N + s - mb) - mb*mp.log(2*a)
    J = J_numpy(beta, a, lams, lw, Y1)
    lpre = lpre - mp.loggamma(mp.mpc(beta)) - mp.pi*tau/2
    return complex(mp.exp(lpre))*J


def Shat_mp_ft(beta, a, dps=70):
    """fixed-twist S^ for regime B (tau > 0.52a): path Y1 = 1.84a, mpmath."""
    mp.mp.dps = dps
    beta = mp.mpc(beta); tau = beta.imag
    Y1 = mp.mpf(1.84)*a
    F = lambda x: mp.exp(a*(mp.log1p(-x/(2*a)) - mp.log1p(x/(2*a))))
    eps = mp.mpf(10)**-25
    tot = mp.exp(beta*mp.log(eps) + 1j*mp.pi/2*beta.real)/beta
    f_leg = lambda u: mp.exp(beta*u + 1j*mp.pi/2*beta.real)*F(1j*mp.exp(u))
    nseg = int(float((mp.log(Y1) - mp.log(eps))*(tau + 1)/4)) + 1
    ed = mp.linspace(mp.log(eps), mp.log(Y1), nseg + 1)
    tot += sum(mp.quad(f_leg, [ed[i], ed[i + 1]]) for i in range(nseg))
    f_arc = lambda ph: mp.exp((beta - 1)*(mp.log(Y1) + 1j*ph) + mp.pi*tau/2)*1j*Y1*mp.exp(1j*ph)*F(Y1*mp.exp(1j*ph))
    nseg = int(float(tau/2)) + 10
    ed = mp.linspace(0, mp.pi/2, nseg + 1)
    tot -= sum(mp.quad(f_arc, [ed[i], ed[i + 1]]) for i in range(nseg))
    f_re = lambda x: mp.exp((beta - 1)*mp.log(x) + mp.pi*tau/2)*F(x)
    tot += mp.quad(f_re, mp.linspace(Y1, 2*a, 40))
    val = mp.exp(-mp.loggamma(beta) - mp.pi*tau/2)*tot
    mp.mp.dps = 30
    return complex(val)


def direct(beta, n, a, b):
    mp.mp.dps = int(abs(beta.imag)*0.7) + 40
    s = 2*b + 1; N = n + a; T = N*(N + s) + b*(b + 1); C = (n + b)*(n + b + 1)
    w = mp.mpf(T - C)/T
    U = mp.mpf(1); poch = mp.mpf(1); tot = mp.mpc(0); bb = mp.mpc(beta)
    for k in range(n + 1):
        tot += U*poch
        poch *= (bb + k)/(k + 1)
        U *= mp.mpf((n - k)*(n + s + k))/((N - k)*(N + s + k))
    v = complex(w**bb*tot)
    mp.mp.dps = 30
    return v, float(-mp.log(w))


if __name__ == '__main__':
    mp.mp.dps = 30
    import tc_power as P
    # ---- validation against the direct ladder sum (n = 3e5, a = 998, b = 0: L = 5.0)
    nval = int(sys.argv[1]) if len(sys.argv) > 1 else 300000
    for beta in (complex(0.25, 10), complex(2.5, 0), complex(1.25, 3)):
        dv, Lv = direct(beta, nval, 998, 0)
        # same state through the representation: build lams for this exact n
        a = 998; s = 1; N = mp.mpf(nval) + a; M = N + nval + s
        mu = (N + 1)/(M + 1); var = mu*(1 - mu)/(M + 2)
        lams = [float(1/(mu + mp.sqrt(2*var)*float(z)) - 1) for z in GHX]
        T = N*(N + s); w = mp.mpf(a)*(N + nval + s)/T; mb = mp.mpc(beta)
        lpre = mb*mp.log(w) + mp.loggamma(N + s) - mp.loggamma(N + s - mb) - mb*mp.log(2*a) - mp.loggamma(mb) - mp.pi*beta.imag/2
        rv = complex(mp.exp(lpre))*J_numpy(beta, a, lams, list(GHW))
        pr(f'validation n={nval}, a=998, L={Lv:.3f}, beta={beta}: direct sum {dv:.10f}  representation {rv:.10f}  |diff| {abs(dv-rv):.2e}')
    # ---- Lemma NB: real and near-real beta
    pr('')
    pr('Lemma NB: |S^ - 1| computed vs certified bound (NB_bound uses a = 998, s = 9)')
    for a in (998, 10000, 1000000):
        for b in (0, 4):
            for L in (18, 20, 25, 40):
                for beta in (2.5, 1.75, complex(1.75, 0.25), complex(0.25, 0.25), 3.5):
                    v = Shat(beta, a, b, L)
                    bd = P.NB_bound(complex(beta).real, complex(beta).imag)
                    ft = (complex(beta).real*(complex(beta).real + 1)*(complex(beta).real + 2))/(12*a*a)
                    if (a == 998 and b in (0, 4)) or L == 20:
                        pr(f'  a={a:>7} b={b} L={L} beta={complex(beta)}: |S^-1| = {abs(v-1):.4e}   bound {bd:.4e}   ratio {abs(v-1)/bd:.4f}'
                           f'   (fixed-twist (s)_3/(12a^2) = {ft:.4e})')
                    assert abs(v - 1) <= bd + 2e-11, (a, b, L, beta, v, bd)
    pr('  -> all NB checks pass (assert)')
    # ---- Theorem P on the leg (sigma = 1/4, 5/4) and k = 3 (sigma = 7/2, 9/2)
    pr('')
    pr('Theorem P: |S^(sigma+it)| (finite n, regime A path) and fixed twist vs certified K_osc(sigma, t, t+0.01, L0=20)')
    worst = 0
    for sig in (0.25, 1.25, 3.5, 4.5):
        for t in (0.25, 1, 3, 10, 45, 100, 300, 500):
            bd, _ = P.K_osc(sig, t, t + 0.01, 20)
            vals = []
            for a in (998, 10000, 1000000):
                if t > 0.52*a:
                    continue
                for b in (0, 4):
                    for L in (20, 25, 40):
                        vals.append(abs(Shat(complex(sig, t), a, b, L)))
            ftv = abs(Shat(complex(sig, t), 998, 0, ft=True))
            worst = max(worst, max(vals)/bd)
            pr(f'  sigma={sig} t={t}: max_states |S^| = {max(vals):.4f}  (min {min(vals):.4f}), fixed twist a=998 {ftv:.4f};  bound {bd:.4g}  '
               f'ratio {max(vals)/bd:.3f}')
    pr(f'  worst |S^|/bound (regime A, finite n) = {worst:.3f}')
    pr('')
    pr('Theorem P, regime B (a = 998, tau > 0.52a): |S^| on the rho = 0.92 path (double precision; the arc amplification e^{a Ehat}')
    pr('  is <= e^{6} for tau <= 1000, so double precision suffices) vs certified K_osc(sigma, t, t+0.01, 20)')
    for sig in (0.25, 1.25, 3.5, 4.5):
        for t in (550, 600, 700, 800, 900, 998):
            vals = [abs(Shat(complex(sig, t), 998, b, L, Y1=1.84*998)) for b in (0, 4) for L in (20, 25)]
            ftv = abs(Shat(complex(sig, t), 998, 0, ft=True, Y1=1.84*998))
            # cross-check of the path independence in the overlap: regime-A path at the same point
            alt = abs(Shat(complex(sig, t), 998, 0, ft=True))
            bd, _ = P.K_osc(sig, t, t + 0.01, 20)
            pr(f'  sigma={sig} t={t}: max_states |S^| = {max(vals):.4f}, fixed twist {ftv:.4f} (regime-A path {alt:.4f});  bound {bd:.4g};  '
               f'ratio {max(vals)/bd:.2e}')
            assert max(vals) <= bd
    open(os.path.join(HERE, 'cache', 'ta_power_check.txt'), 'w').write('\n'.join(out) + '\n')
