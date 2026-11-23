"""
hk_toy (HK_CLOSE.md §4): exactly solvable de-shift toys.  Truncated (compact-support) weights versus analytically cut-off weights.

Toy k = 2 ("one shift", the structure of Result D's C'_EE):
    raw per-nu image  a(nu, beta) = (1/2)[P_{-1/2+i nu+1} + P_{-1/2+i nu-1}](cosh beta)        (grows like e^{beta/2})
    Mehler-Dirichlet: W_chi(t) = int_0^oo chi(nu) (1/2)[cos((nu-i)t) + cos((nu+i)t)] dnu = cosh(t) C[chi](t),
    image I(beta) = (sqrt2/pi) int_0^beta W(t) (cosh beta - cosh t)^(-1/2) dt.
    Exact de-shift to the pure image P_{-1/2}(cosh beta) (W = 1):  C[chi] = sech t  <=>  chi_2(nu) = sech(pi nu/2).
Toy k = 4 ("shifts 1 and 3", the structure of Result Q's C4^ds):
    a(nu, beta) = sum_{j = +-1, +-3} r_j P_{-1/2+i nu+j},  r_{+-1} = 3/8, r_{+-3} = 1/8  =>  W = C[chi] (3/4 cosh t + 1/4 cosh 3t) = C[chi] cosh^3 t.
    Raw growth e^{5 beta/2}.  Exact de-shift (W = 1):  C[chi] = sech^3 t  <=>  chi_4(nu) = (1/2)(1 + nu^2) sech(pi nu/2).
Both exact weights decay like e^{-(pi/2) nu}, the rate measured for the real weights (V4FINAL §2.4: kappa -> pi/2).

Compared, for each toy:
    (T) truncation   chi_T = chi 1[nu < nu_c]                                   (what any stored grid weight is)
    (C) cutoff       chi_C = chi c(nu),  c = (1/2)[erf((nu+nu_c)/a) - erf((nu-nu_c)/a)]   (even, entire, bounded on strips)
    residual dW(t) = W_X(t) - 1 for t in [0, 40], and image residual dI(beta) at beta = 10, 20, 30, 40.
Prediction (HK_CLOSE §3): (T) dW ~ eps e^{(j_max) t}  (raw rate: e^{t} toy 2, e^{3t} toy 4);
                          (C) dW = O(eps) bounded (toy 2) and O(eps e^{(1+0)t}) at worst (toy 4; here even bounded, since the toy chi is
                              analytic in the whole strip |Im nu| < 3).
Usage: OMP_NUM_THREADS=1 python3 hk_toy.py  -> cache/hk_toy.txt   (~2-4 min, one core, < 300 MB)
"""
import os, sys, time
sys.dont_write_bytecode = True
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
OUT = []
mp.mp.dps = 110


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


def chi2(v):
    return mp.sech(mp.pi*v/2)


def chi4(v):
    return (1 + v*v)*mp.sech(mp.pi*v/2)/2


def S2(t):
    return mp.cosh(t)


def S4(t):
    return mp.cosh(t)**3


# ---- closed-form tail  int_{nu_c}^oo chi(nu) cos(nu t) dnu  (sech = 2 sum (-1)^n e^{-(2n+1) x}) ----
def tail(kind, nuc, t, nterms=40):
    s = mp.mpc(0)
    for n in range(nterms):
        z = (2*n + 1)*mp.pi/2 - 1j*t          # int e^{-z nu} e^{i t nu}...: cos = Re e^{i nu t}; e^{-lam nu + i t nu} = e^{-z nu}
        e = mp.e**(-z*nuc)
        if kind == 2:
            I = e/z
        else:   # (1 + nu^2)/2
            I = (e*(1/z) + e*(nuc**2/z + 2*nuc/z**2 + 2/z**3))/2
        s += 2*(-1)**n*I
    return mp.re(s)


# ---- cutoff: 1 - c(nu) on nu > 0 ----
def one_minus_c(v, nuc, a):
    return (mp.erfc((v + nuc)/a) + mp.erfc((nuc - v)/a))/2


_GLC = {}


def legendre_nodes(n):
    """Gauss-Legendre nodes/weights on [-1, 1] at current precision (cached)."""
    if n in _GLC:
        return _GLC[n]
    # mpmath has no public GL generator in all versions; build by Newton on P_n.
    xs, ws = [], []
    for i in range(1, n + 1):
        x = mp.cos(mp.pi*(i - mp.mpf(1)/4)/(n + mp.mpf(1)/2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for k in range(2, n + 1):
                p0, p1 = p1, ((2*k - 1)*x*p1 - (k - 1)*p0)/k
            dp = n*(x*p1 - p0)/(x*x - 1)
            dx = p1/dp
            x -= dx
            if abs(dx) < mp.mpf(10)**(-mp.mp.dps + 5):
                break
        p0, p1 = mp.mpf(1), x
        for k in range(2, n + 1):
            p0, p1 = p1, ((2*k - 1)*x*p1 - (k - 1)*p0)/k
        dp = n*(x*p1 - p0)/(x*x - 1)
        xs.append(x); ws.append(2/((1 - x*x)*dp*dp))
    _GLC[n] = (xs, ws)
    return xs, ws


def cut_residual_nodes(kind, nuc, a, U, h=mp.mpf('0.125'), n=40):
    """nodes nu_k and weights w_k chi(nu)(1 - c(nu)) on [0, U] (composite GL), for int_0^U chi (1-c) cos(nu t)."""
    GLX, GLW = legendre_nodes(n)
    ch = chi2 if kind == 2 else chi4
    nodes, wts = [], []
    m = int(mp.ceil(U/h))
    for i in range(m):
        lo, hi = i*h, min((i + 1)*h, U)
        for x, w in zip(GLX, GLW):
            v = (hi + lo)/2 + (hi - lo)/2*x
            nodes.append(v); wts.append(w*(hi - lo)/2*ch(v)*one_minus_c(v, nuc, a))
    return nodes, wts


def dC_cut(nodes, wts, kind, U, t):
    """int_0^oo chi (c - 1) cos(nu t) = -[ int_0^U chi(1-c) cos + int_U^oo chi cos ]  (1 - c = 1 to < 1e-300 beyond U)."""
    s = mp.fsum(w*mp.cos(v*t) for v, w in zip(nodes, wts))
    return -(s + tail(kind, U, t))


def abel(Wfun, beta, n=120):
    """(sqrt2/pi) int_0^beta W(t)(cosh beta - cosh t)^(-1/2) dt, t = beta - v^2, GL in v."""
    xs, ws = legendre_nodes(n)
    vm = mp.sqrt(beta)
    s = mp.mpf(0)
    for x, w in zip(xs, ws):
        v = vm*(x + 1)/2; wv = vm*w/2; t = beta - v*v
        s += wv*2*v/mp.sqrt(mp.cosh(beta) - mp.cosh(t))*Wfun(t)
    return mp.sqrt(2)/mp.pi*s


def run(kind, nuc, a):
    S = S2 if kind == 2 else S4
    jmax = 1 if kind == 2 else 3
    U = nuc + 14*a
    t0 = time.time()
    nodes, wts = cut_residual_nodes(kind, mp.mpf(nuc), mp.mpf(a), mp.mpf(U))
    ch = chi2 if kind == 2 else chi4
    eps_tail = mp.quad(lambda v: abs(ch(v)), [nuc, mp.inf])
    pr(f'\n== toy k = {2 if kind == 2 else 4}: nu_c = {nuc}, a = {a}; eps_tail = int_{{nu_c}}^oo |chi| = {mp.nstr(eps_tail, 4)}'
       f'  ({len(nodes)} nodes, {time.time()-t0:.0f}s)')
    pr('     t    |dW_T(t)| (truncation)   |dW_T| e^{-jmax t}   |dW_C(t)| (analytic cutoff)   |dW_C|/eps_tail')
    rows = []
    for t in (0, 2, 5, 10, 15, 20, 25, 30, 35, 40):
        t = mp.mpf(t)
        dWT = -S(t)*tail(kind, mp.mpf(nuc), t)
        dWC = S(t)*dC_cut(nodes, wts, kind, U, t)
        rows.append((t, dWT, dWC))
        pr(f'  {float(t):5.1f}   {mp.nstr(abs(dWT), 4):>14}   {mp.nstr(abs(dWT)*mp.e**(-jmax*t), 4):>14}   '
           f'{mp.nstr(abs(dWC), 4):>14}   {mp.nstr(abs(dWC)/eps_tail, 4):>10}')
    # growth rates (log-slope between t = 20 and 40)
    rT = (mp.log(abs(rows[-1][1])) - mp.log(abs(rows[5][1])))/20
    rC = (mp.log(abs(rows[-1][2]) + mp.mpf(10)**-200) - mp.log(abs(rows[5][2]) + mp.mpf(10)**-200))/20
    pr(f'  growth rate of |dW| on t in [20, 40]: truncation {mp.nstr(rT, 4)} (raw rate {jmax}), cutoff {mp.nstr(rC, 4)}')
    # image residuals
    pr('   beta   P_{-1/2}(cosh b)   |dI_T|          |dI_C|          |dI_C|/P_{-1/2}')
    for b in (10, 20, 30, 40):
        b = mp.mpf(b)
        P = mp.legenp(-0.5, 0, mp.cosh(b), type=3)
        dIT = abel(lambda t: -S(t)*tail(kind, mp.mpf(nuc), t), b)
        dIC = abel(lambda t: S(t)*dC_cut(nodes, wts, kind, U, t), b, n=60)
        pr(f'  {float(b):5.1f}   {mp.nstr(mp.re(P), 5):>12}   {mp.nstr(abs(dIT), 4):>12}   {mp.nstr(abs(dIC), 4):>12}   '
           f'{mp.nstr(abs(dIC)/abs(P), 4):>10}')
    pr(f'  [{time.time()-t0:.0f}s]')


if __name__ == '__main__':
    T0 = time.time()
    pr('# hk_toy: truncation vs analytic cutoff of an exact de-shift weight (HK_CLOSE.md §4)')
    pr('# chi_2 = sech(pi nu/2) (W = cosh t C[chi] = 1);  chi_4 = (1+nu^2) sech(pi nu/2)/2 (W = cosh^3 t C[chi] = 1); exact image P_{-1/2}.')
    for kind in (2, 4):
        run(kind, 20, 2)
        run(kind, 30, 2.5)
    pr(f'\n[total {time.time()-T0:.0f}s]')
    open(os.path.join(CACHE, 'hk_toy.txt'), 'w').write('\n'.join(OUT) + '\n')
