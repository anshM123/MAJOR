"""
hk_real (HK_CLOSE.md §3, §5): the de-shift identity at large t on the ACTUAL weights of Results D and Q, read-only.

Q (k = 4; XP weights of functionals_q7_xprec.json, 64 digits, nu_k = (k+1/2)/25, k < 3072; R_j of functionals_q4_deshift.json):
    F_j(t) := int_R chi(nu) R_j(nu) e^{-i nu t} dnu  (chi even, R_j(-nu) = R_{-j}(nu)),   W(t) = (1/2) sum_j e^{-j t} F_j(t)
    G(t)   := e^{3t} F_{-3}(t) + e^{t} F_{-1}(t),   W(t) = (G(t) + G(-t))/2.
    The identity W = Phi on t >= 0 forces |F_{-3}(t)| <= e^{-2t}(sup|G| + ||chi R_{-1}||_1)  (Lemma S4 of HK_CLOSE):
    printed: |F_{-3}| e^{2t}, |F_{-1}|, |G(t) - 2 Phi(t)|, |W - Phi| for t in [0, 45].
D (k = 2; Abel weights of spinfunc/cache/u3_chi_EE_{300,1000}.pkl, double precision):
    fhat(t) := C[aY](t) - i S[bY](t),  E(t) = cosh t C[aY] - i sinh t S[bY] = (e^{t} fhat(t) + e^{-t} fhat(-t))/2,  Phi = C[psi_2].
    Identity forces g(t) := e^{t} fhat(t) = 2 Phi(t) - e^{-t} fhat(-t)  (EE_FOUNDATIONS §3.1).  Printed: |g|, |g - 2Phi|, |E - Phi|.
Also: the flat-model check (HK_CLOSE §6) of per-state contributions beyond the light-cone edge x = nu/m > 1, for weights with an
exponential tail e^{-kappa nu} (uncut) and with the analytic cutoff.
Usage: OMP_NUM_THREADS=1 python3 hk_real.py -> cache/hk_real.txt   (~3-6 min, one core)
"""
import os, sys, json, time, pickle
from fractions import Fraction
sys.dont_write_bytecode = True
import numpy as np
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ADS4 = os.path.dirname(os.path.dirname(HERE))
RQ = os.path.join(ADS4, 'code', 'resultq')
SF = os.path.join(ADS4, 'code', 'spinfunc', 'cache')
CACHE = os.path.join(HERE, 'cache'); os.makedirs(CACHE, exist_ok=True)
OUT = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)


# ------------------------------------------------------------------ Q
def q_part():
    mp.mp.dps = 90
    nus = sp.Symbol('nu')
    Rj = json.load(open(os.path.join(RQ, 'functionals_q4_deshift.json')))['EE3']['R_j']
    Rf = {int(j): sp.lambdify(nus, sp.sympify(v, locals={'nu': nus, 'I': sp.I}), 'mpmath') for j, v in Rj.items()}
    x = json.load(open(os.path.join(RQ, 'functionals_q7_xprec.json')))
    for key, Phi in (('C4_value', lambda t: mp.mpf(1)), ('C4_deriv', lambda t: -t*t/2)):
        t0 = time.time()
        fr = Fraction(x[key]['dnu']); dnu = mp.mpf(fr.numerator)/fr.denominator
        chi = [mp.mpf(v) for v in x[key]['chi']]
        nu = [(k + mp.mpf(1)/2)*dnu for k in range(len(chi))]
        R = {j: [Rf[j](v) for v in nu] for j in Rf}
        l1 = {j: mp.fsum(dnu*abs(c*r) for c, r in zip(chi, R[j]))*2 for j in R}
        pr(f'\n== Q {key}: {len(chi)} nodes, nu_max = {mp.nstr(nu[-1], 6)}; ||chi R_j||_1 (full line): ' +
           ', '.join(f'j={j}: {mp.nstr(l1[j], 4)}' for j in sorted(l1)))
        pr(f'   chi(nu_k)/nu_k^2 at the first nodes: ' + ', '.join(mp.nstr(chi[k]/nu[k]**2, 6) for k in range(4)) +
           '   (chi = O(nu^2): chi R_j regular at 0)')
        pr('      t     |F_-3| e^{2t}     |F_-3| e^{3t}     |F_-1|           |G - 2 Phi|      |W - Phi|')

        def F(j, t):
            s = mp.mpc(0)
            for c, v, r, rm in zip(chi, nu, R[j], R[-j]):
                e = mp.expjpi(-v*t/mp.pi)        # e^{-i v t}
                s += c*(r*e + rm/e)
            return dnu*s
        for t in (0, 1, 2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 35, 40, 45):
            t = mp.mpf(t)
            Fm3, Fm1, F1, F3 = F(-3, t), F(-1, t), F(1, t), F(3, t)
            G = mp.e**(3*t)*Fm3 + mp.e**t*Fm1
            W = (G + mp.e**(-3*t)*F3 + mp.e**(-t)*F1)/2
            pr(f'  {float(t):5.1f}   {mp.nstr(abs(Fm3)*mp.e**(2*t), 5):>12}   {mp.nstr(abs(Fm3)*mp.e**(3*t), 5):>12}   '
               f'{mp.nstr(abs(Fm1), 5):>12}   {mp.nstr(abs(G - 2*Phi(t)), 5):>12}   {mp.nstr(abs(W - Phi(t)), 5):>12}')
        pr(f'   [{time.time()-t0:.0f}s]')



# ------------------------------------------------------------------ Q: cutoff vs truncation applied to the real weights
def qcut_part(nuc=50, a=3):
    mp.mp.dps = 90
    nus = sp.Symbol('nu')
    Rj = json.load(open(os.path.join(RQ, 'functionals_q4_deshift.json')))['EE3']['R_j']
    Rf = {int(j): sp.lambdify(nus, sp.sympify(v, locals={'nu': nus, 'I': sp.I}), 'mpmath') for j, v in Rj.items()}
    x = json.load(open(os.path.join(RQ, 'functionals_q7_xprec.json')))
    nuc, a = mp.mpf(nuc), mp.mpf(a)
    pr(f'\n== Q, real XP weights cut INSIDE their grid at nu_c = {nuc}: truncation chi 1[nu < nu_c] vs analytic cutoff chi c(nu),')
    pr(f'   c = (erf((nu+nu_c)/a) - erf((nu-nu_c)/a))/2, a = {a}.  Prediction: truncation residual ~ chi(nu_c) e^{{3t}}; cutoff residual')
    pr('   bounded (~ eps (1+|Phi|)) until the stored weights\' own floor 1.7e-54 e^{3t} takes over (t ~ 25).')
    for key, Phi in (('C4_value', lambda t: mp.mpf(1)), ('C4_deriv', lambda t: -t*t/2)):
        fr = Fraction(x[key]['dnu']); dnu = mp.mpf(fr.numerator)/fr.denominator
        chi = [mp.mpf(v) for v in x[key]['chi']]
        nu = [(k + mp.mpf(1)/2)*dnu for k in range(len(chi))]
        R = {j: [Rf[j](v) for v in nu] for j in Rf}
        cc = [(mp.erf((v + nuc)/a) - mp.erf((v - nuc)/a))/2 for v in nu]
        tr = [mp.mpf(1) if v < nuc else mp.mpf(0) for v in nu]
        k0 = int(nuc/dnu)
        pr(f'  {key}: |chi(nu_c)| = {mp.nstr(abs(chi[k0]), 3)},  int_{{nu>nu_c}} |chi| = {mp.nstr(mp.fsum(abs(c) for c in chi[k0:])*dnu, 3)}')
        pr('      t     |W_trunc - Phi|   |W_cut - Phi|    |W_full - Phi|')

        def W(mult, t):
            s = mp.mpc(0)
            for c, m, v, r in zip(chi, mult, nu, zip(*(R[j] for j in (-3, -1, 1, 3)))):
                if m == 0:
                    continue
                cv, sv = mp.cos(v*t), mp.sin(v*t)
                for j, rj in zip((-3, -1, 1, 3), r):
                    s += c*m*rj*(cv*mp.cosh(j*t) + 1j*sv*mp.sinh(j*t))
            return dnu*s
        one = [mp.mpf(1)]*len(chi)
        for t in (0, 2, 5, 8, 11, 14, 17, 20, 23, 26, 30):
            t = mp.mpf(t)
            pr(f'  {float(t):5.1f}   {mp.nstr(abs(W(tr, t) - Phi(t)), 4):>12}   {mp.nstr(abs(W(cc, t) - Phi(t)), 4):>12}   '
               f'{mp.nstr(abs(W(one, t) - Phi(t)), 4):>12}')


# ------------------------------------------------------------------ D
def d_part():
    for Dg in (300, 1000):
        d = pickle.load(open(os.path.join(SF, f'u3_chi_EE_{Dg}.pkl'), 'rb'))
        nu, Y, a, b, phi, Lam = d['nu'], d['Y'], d['a'], d['b'], d['phi'], float(d['Lam'])
        dnu = nu[1] - nu[0]
        aY, bY = (a*Y).astype(np.clongdouble), (b*Y).astype(np.clongdouble)
        nl = nu.astype(np.longdouble); ph = phi.astype(np.longdouble)
        pr(f'\n== D, Delta_g = {Dg}: Lambda = {Lam:.3f}, {len(nu)} nodes, nu_max = {nu[-1]:.2f} (double-precision data)')
        pr(f'   |X| at nu = Lam, Lam+10, Lam+20: ' + ', '.join(f'{abs((a*Y)[np.argmin(abs(nu - v))]):.2e}' for v in (Lam, Lam + 10, Lam + 20)))
        pr('      t      |g(t)| = e^t|fhat|   |g - 2 Phi|      |Phi(t)|        |E - Phi|        floor 1e-16 e^t ||X||_1')
        nX = float(np.sum(np.abs(a*Y))*dnu)
        for t in (0, 1, 2, 4, 6, 8, 10, 14, 18, 22, 26, 30, 34, 38):
            tl = np.longdouble(t)
            c, s = np.cos(nl*tl), np.sin(nl*tl)
            C = np.sum(aY*c)*dnu; S = np.sum(bY*s)*dnu; P = np.sum(ph*c)*dnu
            fh = C - 1j*S
            g = np.exp(tl)*fh
            E = np.cosh(tl)*C - 1j*np.sinh(tl)*S
            pr(f'  {t:5.1f}   {abs(g):14.5e}   {abs(g - 2*P):12.4e}   {abs(P):12.4e}   {abs(E - P):12.4e}   {1e-16*np.exp(t)*nX:10.2e}')


# ------------------------------------------------------------------ flat-model x > 1 check
def flat_part():
    from scipy.special import log_ndtr

    def log_half_erfc(z):          # log( erfc(z)/2 ) = log_ndtr(-sqrt2 z), no underflow
        return log_ndtr(-np.sqrt(2.0)*z)
    pr('\n== flat model beyond the light-cone edge (HK_CLOSE §6) [M]: per-state image at x = nu/m > 1 grows like e^{2J arccosh x};')
    pr('   twist gap => J <= (m^2/(Dg-2) - Dg + 1)/2.  Worst exponent  max_{m >= Lam, nu > m} [ log w(nu) + 2 J_max(m) arccosh(nu/m) ];')
    pr('   w(nu) = weight envelope relative to its O(1) bulk.  m_crit = smallest m with a positive exponent (per-state tail > O(1)).')
    for Dg in (300, 1000):
        Lam = np.sqrt((Dg - 2)*(Dg - 1))
        ms = Lam*np.linspace(1.0, 8.0, 1401)
        Jm = np.maximum(0, (ms**2/(Dg - 2) - Dg + 1)/2)
        cases = (('uncut tail, kappa = pi/2', lambda v: -np.pi/2*v),
                 ('uncut tail, kappa = 1.46 (measured at nu ~ 100)', lambda v: -1.46*v),
                 ('Q: cutoff nu_c = 150, a = 6', lambda v: -np.pi/2*v + log_half_erfc((v - 150)/6)),
                 ('D: tail 3e-9 e^{-pi/2 (nu-Lam)}, cutoff Lam + 100, a = 10',
                  lambda v: np.log(3e-9) - np.pi/2*np.maximum(0, v - Lam) + log_half_erfc((v - Lam - 100)/10)),
                 ('D: same tail, NO cutoff', lambda v: np.log(3e-9) - np.pi/2*np.maximum(0, v - Lam)))
        for label, logw in cases:
            best = (-np.inf, None, None); mcrit = None
            x = 1.0 + np.geomspace(1e-6, 11.0, 6000)
            for m, J in zip(ms, Jm):
                val = logw(x*m) + 2*J*np.arccosh(x)
                k = np.argmax(val)
                if val[k] > 0 and mcrit is None:
                    mcrit = m/Lam
                if val[k] > best[0]:
                    best = (val[k], m/Lam, x[k])
            mc = f'{mcrit:.2f} Lam' if mcrit is not None else 'none (m <= 8 Lam)'
            pr(f'   Dg = {Dg}: {label:55s}: worst exponent {best[0]:+10.1f} (m = {best[1]:.2f} Lam, x = {best[2]:.3f}); m_crit = {mc}')


if __name__ == '__main__':
    T0 = time.time()
    pr('# hk_real: de-shift identity at large t on the stored weights; flat-model x > 1 check (HK_CLOSE.md)')
    parts = sys.argv[1:] or ['q', 'qcut', 'd', 'flat']
    if 'q' in parts:
        q_part()
    if 'd' in parts:
        d_part()
    if 'flat' in parts:
        flat_part()
    if 'qcut' in parts:
        qcut_part()
    pr(f'\n[total {time.time()-T0:.0f}s]')
    open(os.path.join(CACHE, 'hk_real' + ('' if len(sys.argv) < 2 else '_' + '_'.join(sys.argv[1:])) + '.txt'), 'w').write('\n'.join(OUT) + '\n')
