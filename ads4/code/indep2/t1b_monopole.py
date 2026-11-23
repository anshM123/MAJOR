"""Target 1b: the monopole cost formula of PW_IMAGES (4.1),
     cost(A_g[F]) := Re sum_k dnu X_k A[F/c](nu_k)/Lam  =  sum_k dnu psi2 F/c /Lam  +  D F(i/2),
with A[f] = (R_- f(nu+i) + R_+ f(nu-i))/(R_- + R_+) (EE kernels of prl/data/kernels.json), X, psi2 from
prl/data/d_abel_X_{1000,300}.txt, c(mu) = 32(4mu^2+1)(4mu^2+9)/((4mu^2+25)(4mu^2+81)).
Independent route: partial fractions 1/c = 1/32 + (15/2)/(4mu^2+1) - (9/2)/(4mu^2+9), so
   defect[P_t/c] = (15/2) defect[P_t/(4mu^2+1)] + (analytic pieces, defect 0)
and P_t(+-i/2) = 1 makes defect[(P_t - 1)/(4mu^2+1)] = 0, hence D = (15/2) defect[1/(4mu^2+1)], t-independent.
Output: cache/t1b_monopole.txt
"""
import json, time
import numpy as np
import sympy as sp
import mpmath as mp

DATA = '../../prl/data/'
out = open('cache/t1b_monopole.txt', 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

kj = json.load(open(DATA + 'kernels.json'))['D_k2_EE']
nus = sp.Symbol('nu')
Rm = sp.lambdify(nus, sp.sympify(kj['R_minus'], locals={'nu': nus}), 'numpy')
Rp = sp.lambdify(nus, sp.sympify(kj['R_plus'], locals={'nu': nus}), 'numpy')

def cfun(mu):
    m2 = 4 * mu * mu
    return 32 * (m2 + 1) * (m2 + 9) / ((m2 + 25) * (m2 + 81))

def load(Dg):
    a = np.loadtxt(DATA + f'd_abel_X_{Dg}.txt')
    nu = a[:, 0]; X = a[:, 1] + 1j * a[:, 2]; psi2 = a[:, 5]
    dnu = nu[1] - nu[0]
    Lam = np.sqrt((Dg - 2.0) * (Dg - 1.0))
    return nu, X, psi2, dnu, Lam

# ---------------------------------------------------------------- Legendre P_{-1/2+i mu}(cosh t), complex mu (Mehler)
def mehlerP(mu, t, nper=4.0, nmin=200, chunk=512):
    """P_{-1/2+i mu}(cosh t) = (sqrt2/pi) int_0^t cos(mu s)/sqrt(cosh t - cosh s) ds, s = t - u^2."""
    mu = np.asarray(mu, complex)
    if t == 0:
        return np.ones_like(mu)
    n = int(max(nmin, nper * np.max(np.abs(mu.real)) * t / np.pi + 100))
    x, w = np.polynomial.legendre.leggauss(n)
    U = np.sqrt(t)
    u = 0.5 * U * (x + 1); wu = 0.5 * U * w
    s = t - u * u
    # cosh t - cosh s = 2 sinh((t+s)/2) sinh((t-s)/2), (t - s) = u^2
    den = np.sqrt(2 * np.sinh((t + s) / 2) * np.sinh(u * u / 2))
    g = np.where(u > 0, 2 * u / np.where(den > 0, den, 1), 2 / np.sqrt(np.sinh(t)))
    wt = wu * g * np.sqrt(2) / np.pi
    outv = np.empty_like(mu)
    for i0 in range(0, mu.size, chunk):
        m = mu[i0:i0 + chunk]
        outv[i0:i0 + chunk] = np.cos(np.outer(m, s)) @ wt
    return outv

mp.mp.dps = 20
for (muv, t) in [(3.7 + 1j, 2.0), (250.3 - 1j, 1.0), (0.4 + 0j, 6.0), (1000.2 + 1j, 0.3)]:
    ref = complex(mp.legenp(-0.5 + 1j * muv, 0, mp.cosh(t), type=3))
    val = mehlerP(np.array([muv]), t)[0]
    log(f'Mehler vs mpmath legenp: mu={muv}, t={t}: rel err {abs(val/ref-1):.1e}')

def make_functional(Dg):
    nu, X, psi2, dnu, Lam = load(Dg)
    rm, rp = Rm(nu + 0j), Rp(nu + 0j)
    def cost(fp, fm):
        """Re sum X (R_- f(nu+i) + R_+ f(nu-i))/(R_- + R_+) dnu / Lam,  fp = f(nu+i), fm = f(nu-i)."""
        return float(np.real(np.sum(X * (rm * fp + rp * fm) / (rm + rp))) * dnu / Lam)
    def target(f0):
        return float(np.sum(psi2 * f0) * dnu / Lam)
    return nu, Lam, cost, target

for Dg in (1000, 300):
    t0 = time.time()
    nu, Lam, cost, target = make_functional(Dg)
    log(f'\n===== Delta_g = {Dg}: Lam = {Lam:.10f}, grid {len(nu)} points, nu_max = {nu[-1]:.2f}')
    # rational test functions: exact values at nu +- i
    def defect_rat(h):
        return cost(h(nu + 1j), h(nu - 1j)) - target(h(nu + 0j))
    d1 = defect_rat(lambda z: 1 / (4 * z * z + 1))
    d9 = defect_rat(lambda z: 1 / (4 * z * z + 9))
    dconst = defect_rat(lambda z: 0 * z + 1)
    dc = defect_rat(lambda z: 1 / cfun(z))
    D_pf = 7.5 * d1
    log(f'defect[1] (constant, analytic)            = {dconst:+.3e}   (target value {target(np.ones_like(nu)):.6f})')
    log(f'defect[1/(4nu^2+9)] (analytic in strip)   = {d9:+.3e}')
    log(f'defect[1/(4nu^2+1)] (poles at +-i/2)      = {d1:+.9e}')
    log(f'D from partial fractions = 7.5*defect[1/(4nu^2+1)] = {D_pf:+.7e};  direct defect[1/c] = {dc:+.7e};  D*Lam^2 = {D_pf*Lam**2:.1f}')
    # Legendre test functions
    ts = [0.0, 0.01, 0.03, 0.1, 0.3, 1.0, 2.0, 3.0, 6.0] if Dg == 1000 else [0.0, 0.1, 1.0, 3.0]
    log(' t      defect[P_t] (de-shift id.)   Omega2(t)=int psi2 P_t   Omega_c(t)=int psi2 P_t/c   defect[P_t/c]     Omega~(t)=Omega_c+D   defect[(P_t-1)/(4nu^2+1)]')
    Ps = {}
    for t in ts:
        P0 = mehlerP(nu + 0j, t); Pp = mehlerP(nu + 1j, t); Pm = mehlerP(nu - 1j, t)
        Ps[t] = (P0, Pp, Pm)
        dP = cost(Pp, Pm) - target(P0)
        Om2 = target(P0)
        cp, cm, c0 = cfun(nu + 1j), cfun(nu - 1j), cfun(nu + 0j)
        Omc = target(P0 / c0)
        dPc = cost(Pp / cp, Pm / cm) - Omc
        f = lambda z, P: (P - 1) / (4 * z * z + 1)
        dloc = cost(f(nu + 1j, Pp), f(nu - 1j, Pm)) - target(f(nu + 0j, P0))
        log(f'{t:5.2f}  {dP:+.3e}               {Om2:+.4e}             {Omc:+.4e}               {dPc:+.7e}   {Omc + D_pf:+.4e}           {dloc:+.2e}')
    # mixtures F = sum w_k (1 + a_k nu^2) P_{t_k}
    rng = np.random.default_rng(11)
    tl = [t for t in ts if t > 0]
    for trial in range(4):
        idx = rng.choice(len(tl), size=3, replace=False)
        w = rng.normal(size=3); a = rng.uniform(-0.3, 0.3, size=3)
        Fp = sum(w[i] * (1 + a[i] * (nu + 1j)**2) * Ps[tl[k]][1] for i, k in enumerate(idx))
        Fm = sum(w[i] * (1 + a[i] * (nu - 1j)**2) * Ps[tl[k]][2] for i, k in enumerate(idx))
        F0 = sum(w[i] * (1 + a[i] * nu**2) * Ps[tl[k]][0] for i, k in enumerate(idx))
        lhs = cost(Fp / cfun(nu + 1j), Fm / cfun(nu - 1j))
        Fhalf = sum(w[i] * (1 - a[i] / 4) for i in range(3))
        rhs = target(F0 / cfun(nu + 0j)) + D_pf * Fhalf
        log(f'mixture {trial}: t = {[tl[k] for k in idx]}, w = {np.round(w,3)}, a = {np.round(a,3)}: cost = {lhs:+.9e}, formula = {rhs:+.9e}, rel diff {abs(lhs-rhs)/max(abs(lhs),abs(rhs)):.1e}, F(i/2) = {Fhalf:+.4f}')
    log(f'[{time.time()-t0:.0f}s]')
