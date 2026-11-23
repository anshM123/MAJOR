"""
pwimg/pw1_magic.py (PW_IMAGES.md §2): the magic-block identity for the b-shifted (0, b) continuous-spin blocks that make up the EE
kernel (EESTEP_CLOSE §1), and the resulting exact kernel-level de-shift of the raw EE functional.

Claim (M_b):  for b = 0..4 and every complex l = -1/2 + i nu,
     2 eta (1 + r^2) Gh^{Dp=1}_b(r, eta; nu) = (1 + i/(2nu)) Gh^{Dp=0}_b(r, eta; nu + i) + (1 - i/(2nu)) Gh^{Dp=0}_b(r, eta; nu - i)
as power series in r (Gh = sum_n r^n sum_m B[n,m] P_{l+m}(eta); 1 + i/(2nu) = 2l/(2l+1), 1 - i/(2nu) = 2(l+1)/(2l+1)).
Since f = v' - u' = eta (1+r^2)/(2r) in the H7 frame (checked here), (M_b) for all b gives, with E^{Dp}_mu = sum_b C_b v^{b/2} Gh^{Dp}_b:
     f K^EE_nu = N(nu) (r/4) [ (1 + i/2nu) E^0_{nu+i} + (1 - i/2nu) E^0_{nu-i} ]        (every polarisation set, all orders in r).
Checks:
 (a) EXACT: sympy, l symbolic, b = 0..4, all coefficients of r^n P_{l+m}, n <= NSYM (rational-function identity -> [P] to that order).
 (b) mpmath 50 digits, n <= 40, random complex l, b = 0..4 (and b = 7, 10 as a b-independence probe).
 (c) end to end: the directly integrated EE kernel K^EE_nu (eestep/ee_quad.py, 2-d quadrature of the light-ray integral) against
     N(nu)(r/4)/f [alpha E^0_{nu+i} + beta E^0_{nu-i}] with the Dp = 0 blocks evaluated at complex nu +- i; N(nu) from one point.
 (d) analyticity: Gh^{Dp=0}_b(mu) at the would-be poles mu = +-i/2, +-3i/2 (zeros of the Casimir denominators) and evenness in mu.
Output: cache/pw1_magic.txt
"""
import os, sys, io, time
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import mpmath as mp
import pw_common as PC

NSYM = int(os.environ.get('PW_NSYM', '7'))
OUT = io.StringIO()


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.write(s + '\n')


# ------------------------------------------------------------------ generic recursion (same algebra as ee_blocks.block_coeffs_b)
def leg1(j):
    f_up = (j + 1)/(2*j + 1); f_dn = j/(2*j + 1); c = j*(j + 1)/(2*j + 1)
    return f_up, f_dn, -c, c


def leg2(j):
    e_up = (j + 1)*(j + 2)/((2*j + 1)*(2*j + 3))
    e_mid = (j + 1)**2/((2*j + 1)*(2*j + 3)) + j**2/((2*j + 1)*(2*j - 1))
    e_dn = j*(j - 1)/((2*j + 1)*(2*j - 1))
    c = j*(j + 1)/(2*j + 1)
    return e_up, e_mid, e_dn, -c*(j + 2)/(2*j + 3), c*j/(2*j - 1) - c*(j + 1)/(2*j + 3), c*(j - 1)/(2*j - 1)


def coeffs(Jp, b, N, Dp, simp=lambda x: x):
    """dict B[(n, m)] for the (Dp, Jp) continuous-spin (0, b) block (ee_blocks recursion, generic arithmetic)."""
    C = Dp*(Dp - 3) + Jp*(Jp + 1)
    B = {(0, 0): 1}
    get = lambda n, m: B.get((n, m), 0)
    for n in range(1, N + 1):
        rhs = {}
        def add(m, v):
            rhs[m] = rhs.get(m, 0) + v
        for (dn, kind) in ((2, 'a'), (4, 'b'), (6, 'c')):
            if n < dn:
                continue
            E = Dp + n - dn
            for m in range(-(n - dn), n - dn + 1):
                v = get(n - dn, m)
                if v == 0:
                    continue
                j = Jp + m
                K = E*(E - 3) + j*(j + 1) - C
                if kind == 'c':
                    add(m, -(K/2 + 3*E)*v); continue
                a2, aD, a1 = ((-4*(K/2 + E), 4, (K/2 + E)) if kind == 'a' else (4*(K/2 + 2*E), -4, -(K/2 + 2*E)))
                eu, em, ed, du, dm, dd = leg2(j)
                add(m, (a2*em + aD*dm + a1)*v); add(m + 2, (a2*eu + aD*du)*v); add(m - 2, (a2*ed + aD*dd)*v)
        if b != 0:
            for (dn, ae, aD) in ((1, None, 4*b), (3, None, 0), (5, None, -4*b)):
                if n < dn:
                    continue
                E = Dp + n - dn
                aeta = {1: -4*b*E, 3: 8*b*E, 5: -4*b*E}[dn]
                for m in range(-(n - dn), n - dn + 1):
                    v = get(n - dn, m)
                    if v == 0:
                        continue
                    fu, fd, gu, gd = leg1(Jp + m)
                    add(m + 1, (aeta*fu + aD*gu)*v); add(m - 1, (aeta*fd + aD*gd)*v)
        for m in range(-n, n + 1):
            if (m - n) % 2 or m not in rhs:
                continue
            E = Dp + n; j = Jp + m
            K = E*(E - 3) + j*(j + 1) - C
            B[(n, m)] = simp(-rhs[m]/(K/2))
    return B


def eta_times(Bd, Jp, simp=lambda x: x):
    """coefficients of 2 eta (1 + r^2) sum B[n,m] r^n P_{Jp+m}."""
    out = {}
    for (n, m), v in Bd.items():
        fu, fd, _, _ = leg1(Jp + m)
        for nn in (n, n + 2):
            out[(nn, m + 1)] = out.get((nn, m + 1), 0) + 2*fu*v
            out[(nn, m - 1)] = out.get((nn, m - 1), 0) + 2*fd*v
    return {k: simp(v) for k, v in out.items()}


def identity_residual(Jp, b, N, simp=lambda x: x, one=1):
    lhs = eta_times(coeffs(Jp, b, N, 1*one, simp), Jp, simp)
    Bp = coeffs(Jp - 1, b, N, 0*one, simp)      # nu + i  <->  l - 1: P_{l-1+m'} -> m = m' - 1
    Bm = coeffs(Jp + 1, b, N, 0*one, simp)      # nu - i  <->  l + 1
    al = 2*Jp/(2*Jp + 1); be = 2*(Jp + 1)/(2*Jp + 1)
    res = {}
    keys = set(k for k in lhs if k[0] <= N) | set((n, m - 1) for (n, m) in Bp) | set((n, m + 1) for (n, m) in Bm)
    for (n, m) in keys:
        if n > N:
            continue
        r = lhs.get((n, m), 0) - al*Bp.get((n, m + 1), 0) - be*Bm.get((n, m - 1), 0)
        res[(n, m)] = simp(r)
    return res, lhs


def main():
    t0 = time.time()
    pr(__doc__.split('Output')[0])
    # f = v' - u' in the H7 frame
    worst = 0
    for (r, beta) in ((0.1, 0.5), (0.3, 1.2), (0.05, 2.0), (0.6, 0.1)):
        xh = np.array([np.cosh(beta), np.sinh(beta), 0.]); yh = np.array([1., 0, 0])
        md = lambda a, c: -(-(a - c)[0]**2 + (a - c)[1]**2 + (a - c)[2]**2)
        x1, x3, x2, x4 = -r*yh, r*yh, -xh, xh
        u = md(x1, x2)*md(x3, x4)/(md(x1, x3)*md(x2, x4)); v = md(x1, x4)*md(x2, x3)/(md(x1, x3)*md(x2, x4))
        worst = max(worst, abs((v - u)/(np.cosh(beta)*(1 + r*r)/(2*r)) - 1))
    pr(f'(0) f = v\'-u\' = eta(1+r^2)/(2r) in the H7 frame: max rel dev {worst:.1e}')
    # (a) symbolic
    l = sp.Symbol('l')
    simp = lambda x: sp.cancel(sp.together(x))
    for b in range(5):
        ts = time.time()
        res, lhs = identity_residual(l, b, NSYM, simp, one=sp.Integer(1))
        nz = [k for k, v in res.items() if v != 0]
        pr(f'(a) b = {b}: n <= {NSYM}: {len(res)} coefficients of r^n P_(l+m), {len(nz)} nonzero residuals (exact rational functions of l)'
           f'  [{time.time() - ts:.0f}s]')
    # (b) high precision numerics
    mp.mp.dps = 50
    rng = np.random.default_rng(5)
    for b in (0, 1, 2, 3, 4, 7, 10):
        worst = 0; scale = 0
        for _ in range(2):
            nu = mp.mpc(rng.uniform(0.2, 8), rng.uniform(-0.9, 0.9))
            Jp = mp.mpf(-0.5) + 1j*nu
            res, lhs = identity_residual(Jp, b, 40, one=mp.mpf(1))
            for (n, m), v in res.items():
                sc = abs(lhs.get((n, m), 0)) + 1e-300
                worst = max(worst, float(abs(v)/max(sc, mp.mpf(10)**-30)))
        pr(f'(b) b = {b:2d}: n <= 40, 50 digits, 2 random complex nu: max |residual|/|lhs coefficient| = {worst:.1e}')
    # (d) analyticity at the Casimir zeros, evenness
    pr('(d) Dp = 0 blocks near mu0 = -i/2, -3i/2 (rho = 0.4, beta = 0.8): |Gh(mu0 + eps)| for eps = 1e-2, 1e-4, 1e-6 (finite => removable)')
    r, beta = 0.4*np.exp(-0.8), 0.8
    for mu0 in (-0.5j, -1.5j, 0.5j, 1.5j):
        row = []
        for b in range(5):
            row.append('/'.join(f'{abs(PC.Gh_val(mu0 + e, b, r, beta, 50)):.6f}' for e in (1e-2, 1e-4, 1e-6)))
        pr(f'    mu0 = {mu0}: ' + '  '.join(f'b={b}: {x}' for b, x in enumerate(row)))
    ev = max(abs(PC.Gh_val(mu, b, r, beta, 50) - PC.Gh_val(-mu, b, r, beta, 50))/abs(PC.Gh_val(mu, b, r, beta, 50))
             for mu in (0.7, 1.3 + 0.4j, 0.2j, 2.9 - 0.7j) for b in range(5))
    pr(f'    evenness Gh^0_b(mu) = Gh^0_b(-mu): max rel dev {ev:.1e} (4 mu, b = 0..4)')
    # (c) end to end with the direct quadrature
    pr('(c) end to end: K^EE_nu (ee_quad direct 2-d quadrature)  vs  N(nu) (r/4)/f [a E0(nu+i) + b E0(nu-i)]  (N(nu) fixed at point 0)')
    rng = np.random.default_rng(23)
    pts = [(0.25, 0.6), (0.55, 1.4), (0.85, 2.2), (0.4, 0.3)]
    for ip in range(2):
        pa = tuple(rng.uniform(0, 2*np.pi, 4))
        for nu in (0.7, 2.5, 6.0):
            ratios = []
            for (rho, beta) in pts:
                r = rho*np.exp(-beta)
                pol = PC.Q.pols(beta, *pa)
                parts = PC.Q.kernel_parts(r, beta, nth=128, ds=0.02)
                cb, tail, odd, q = PC.Q.laurent_cb(r, beta, pol)
                ins = PC.Q.dressing(parts['th'], r, beta, pol)
                K = PC.Q.integrate(parts, [nu], ins)[0]
                f = np.cosh(beta)*(1 + r*r)/(2*r)
                Ep, _ = PC.E_val(nu + 1j, r, beta, pol, N=140)
                Em, _ = PC.E_val(nu - 1j, r, beta, pol, N=140)
                pred = (r/4)/f*((1 + 0.5j/nu)*Ep + (1 - 0.5j/nu)*Em)
                ratios.append(K/pred)
            ratios = np.array(ratios)
            pr(f'    pol set {ip}, nu = {nu}: K/pred at 4 points / first: max dev {np.max(np.abs(ratios/ratios[0] - 1)):.1e}')
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(PC.CACHE, 'pw1_magic.txt'), 'w').write(OUT.getvalue())


if __name__ == '__main__':
    main()
