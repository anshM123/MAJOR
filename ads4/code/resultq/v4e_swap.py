"""
v4e_swap (RESULT_Q_V4FINAL.md §2-§3): the order of the nu-integral, the heavy OPE sum and the light side.

(1) Per heavy state (leading Regge, heavy-normalised): the nu-outside image
        I(beta) = int dnu chi(nu) sum_j R_j(nu) P_{-1/2 + i(nu - i j)}(cosh beta)        (= P_{-1/2}(cosh beta) for the value weight)
    versus the absolute integral  A(beta) = int dnu |chi| sum_j |R_j| |P_{-1/2 + i(nu - i j)}(cosh beta)|,  and the raw per-nu image
    max_nu |sum_j R_j P(nu - i j)|.  A/|I| measures what Fubini over (nu, states) costs: the state sum of A converges iff the RAW corner
    does (J0 < 1/2 under K1'), the state sum of I iff the de-shifted one does (J0 < 7/2).
    P_{-1/2+i mu}(cosh b) = (sqrt2/pi) int_0^b cos(mu t)(cosh b - cosh t)^(-1/2) dt (Mehler-Dirichlet, complex mu), t = b - v^2.
(2) Position-normalised weight w = chi N_com(0)/N_com(nu): local growth rate of log|w| vs pi/2 (the celestial-sector argument).
(3) H_nuout for the J0 >= 4 EFT light sides (P): exact moments sum_j int chi R_j (nu - i j)^(2n) = (-1)^n Phi^(2n)(0) and the ABSOLUTE
    moments int |chi| sum_j |R_j| |nu - i j|^(2n) (finite => the nu-outside light side of any polynomial (P)-term converges absolutely).
(4) Single-family obstruction vs summed light side: the MFT [TT]_0 edge term I_t(nu) (EE closed form R_t, V4prime) per nu, heavy-
    normalised, and the truncated nu-outside de-shift sum; versus I_s + I_t, whose integrand R_s + R_t vanishes pointwise.
Usage: python3 v4e_swap.py -> cache/v4e_swap.txt (~1 min)
"""
import os, sys, time, json
sys.dont_write_bytecode = True
for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(v, '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import sympy as sp
import mpmath as mp
import v4d_core as D

CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
L = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def Rfuns():
    nus = sp.Symbol('nu')
    Rj = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']['R_j']
    return {int(j): sp.lambdify(nus, sp.sympify(v, locals={'nu': nus, 'I': sp.I}), 'numpy') for j, v in Rj.items()}


def Pleg(mu, b, nq=600):
    """P_{-1/2 + i mu}(cosh b) for an array of complex mu (Mehler-Dirichlet, t = b - v^2, Gauss-Legendre in v)."""
    x, w = np.polynomial.legendre.leggauss(nq)
    vmax = np.sqrt(b)
    v = 0.5*vmax*(x + 1); wv = 0.5*vmax*w
    t = b - v**2
    den = np.sqrt(np.cosh(b) - np.cosh(t))
    with np.errstate(all='ignore'):
        jac = np.where(v > 0, 2*v/den, 2/np.sqrt(np.sinh(b)))
    return (np.sqrt(2)/np.pi)*(np.cos(np.outer(mu, t)) @ (wv*jac))


def part1(R):
    pr('(1) per-state nu-outside image vs absolute integral (XP weights; heavy-normalised leading-Regge images)')
    W = D.weights()
    for nm in ('xp_value', 'xp_deriv'):
        nu, chi, dnu = W[nm]
        sel = np.abs(chi)*(1 + nu)**8 > 1e-40
        nu, chi = nu[sel], chi[sel]
        Rv = {j: R[j](nu) for j in R}
        pr(f'  {nm}: nu <= {nu[-1]:.1f} ({len(nu)} nodes)')
        pr('      beta   I(beta)          target            A(beta)       A/|target|   |raw image| at nu = 1, 4 (heavy-normalised, per nu)')
        for b in (0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            nq = int(200 + 60*b*np.sqrt(nu[-1]))
            img = 0j; A = 0.0; raw = np.zeros(len(nu), complex)
            for j in (-3, -1, 1, 3):
                P = Pleg(nu - 1j*j, b, nq=nq)
                img += dnu*np.sum(chi*Rv[j]*P)
                A += dnu*np.sum(np.abs(chi)*np.abs(Rv[j])*np.abs(P))
                raw += Rv[j]*P
            if nm == 'xp_value':
                tgt = complex(mp.legenp(-0.5, 0, mp.cosh(b), type=3))
            else:   # Phi = -t^2/2  ->  -D(beta) = -(sqrt2/2pi) int_0^b t^2 (cosh b - cosh t)^(-1/2) dt
                x, w = np.polynomial.legendre.leggauss(400); vm = np.sqrt(b); v = 0.5*vm*(x + 1); wv = 0.5*vm*w; t = b - v**2
                with np.errstate(all='ignore'):
                    jac = np.where(v > 0, 2*v/np.sqrt(np.cosh(b) - np.cosh(t)), 2/np.sqrt(np.sinh(b)))
                tgt = -(np.sqrt(2)/np.pi)*np.sum(wv*jac*t**2/2)
            k1 = np.argmin(np.abs(nu - 1.0)); k4 = np.argmin(np.abs(nu - 4.0))
            pr(f'     {b:4.1f}  {img.real:+.8e}  {tgt.real:+.8e}  {A:.4e}   {A/abs(tgt):.3e}      {abs(raw[k1]):.3e}  {abs(raw[k4]):.3e}')
    pr('    growth: A(beta) and |raw| ~ e^{5 beta/2} (|P_{5/2+i nu}|); I(beta) ~ e^{-beta/2} beta (value), beta^3 e^{-beta/2} (deriv).')


def part2():
    pr('(2) position-normalised weight w = chi N(0)/N(nu): local exponential rates (d log|.|/dnu over windows)')
    W = D.weights()
    lN0 = D.logNcom(0.0)
    for nm in ('xp_value', 'xp_deriv'):
        nu, chi, dnu = W[nm]
        lw = np.log(np.abs(chi) + 1e-300) + lN0 - D.logNcom(nu)
        lc = np.log(np.abs(chi) + 1e-300)
        out = []
        for a, b in ((10, 20), (20, 30), (30, 40), (50, 60), (70, 80), (90, 100), (100, 110)):
            m = (nu >= a) & (nu <= b) & (np.abs(chi) > 0)
            if m.sum() < 10:
                continue
            sc = np.polyfit(nu[m], lc[m], 1)[0]; sw = np.polyfit(nu[m], lw[m], 1)[0]
            out.append(f'[{a},{b}]: chi {-sc:.3f}, w {sw:+.3f}')
        pr(f'    {nm}: ' + ';  '.join(out) + f'   (pi/2 = {np.pi/2:.3f}, pi = {np.pi:.3f})')
    pr('    reading: |w| grows at rate kappa = pi - rate(chi) > pi/2 on the whole computed range; a nu-inside (kernel-level) position kernel'
       ' int w C_nu dnu exists only as an analytic functional needing a celestial-angle rotation theta > kappa > pi/2, i.e. beyond the'
       ' (13) light-cone cut arg x13^2 = pi.')


def part3(R):
    pr('(3) H_nuout for polynomial (P) light sides: exact and absolute moments of the de-shift weights')
    W = D.weights()
    for nm in ('xp_value', 'xp_deriv'):
        nu, chi, dnu = W[nm]
        Rv = {j: R[j](nu) for j in R}
        ex, ab = [], []
        for n in range(4):
            ex.append(sum(dnu*np.sum(chi*Rv[j]*(nu - 1j*j)**(2*n)) for j in R))
            ab.append(sum(dnu*np.sum(np.abs(chi)*np.abs(Rv[j])*np.abs(nu - 1j*j)**(2*n)) for j in R))
        pr(f'    {nm}: exact moments n=0..3: ' + ' '.join(f'{x.real:+.6f}' for x in ex) +
           ';  absolute moments: ' + ' '.join(f'{x:.3e}' for x in ab))
    pr('    targets: value (1, 0, 0, 0); deriv (0, 1, 0, 0) [(-1)^n Phi^(2n)(0)].  Finite absolute moments: |L_nu| <= sum_terms |c| sum_j |R_j|'
       ' |p(nu - ij)| is polynomially bounded in heavy normalisation and its nu-outside de-shift converges absolutely.')


def part4():
    pr('(4) single light family vs the summed light side (MFT, k = 4, EE): per-nu I_t(nu) = int C_nu R_t, heavy-normalised')
    lN0 = D.logNcom(0.0)
    Rt = lambda x: D.RtEE(np.cos(x))
    Rs = lambda x: -D.RtEE(np.cos(x))
    rows = []
    for nu in (0.5, 2, 5, 8, 12, 16, 20, 25, 30):
        It = D.pair_exact(float(nu), Rt, N=2**14)
        Is = D.pair_exact(float(nu), Rs, N=2**14)
        rows.append((nu, It, Is))
        pr(f'    nu = {nu:5.1f}: |I_t| = {abs(It):.3e}   |I_t| N(0)/N(nu) = {abs(It)*np.exp(lN0 - D.logNcom(nu)):.3e}   |I_s + I_t| = {abs(Is + It):.1e}')
    W = D.weights()
    for nm in ('xp_value',):
        nu, chi, dnu = W[nm]
        w = chi*np.exp(lN0 - D.logNcom(nu))
        sel = nu <= 40.0
        It = np.array([D.pair_exact(float(v), Rt, N=2**13) for v in nu[sel][::5]])
        nus = nu[sel][::5]; ws = w[sel][::5]; dd = 5*dnu
        out = []
        for nuc in (5, 10, 20, 30, 40):
            m = nus <= nuc
            out.append(f'nu_c={nuc}: {abs(dd*np.sum(ws[m]*It[m])):.3e}')
        pr(f'    nu-outside de-shift of I_t alone ({nm}), truncated: ' + '  '.join(out) + '   -> diverges with nu_c;'
           ' the summed light side R_s + R_t = 0 pointwise gives exactly 0 for every weight and every order of integration.')


def main():
    t0 = time.time()
    pr('# v4e_swap: nu-order, Fubini over heavy states, H_nuout')
    R = Rfuns()
    part1(R)
    part2()
    part3(R)
    part4()
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(CACHE, 'v4e_swap.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
