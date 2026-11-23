"""Target A (iii): f_pos(+-3i/2) = 0 for Q's published k = 4 weights (ads4/prl/data), independent evaluation.
 (1) exact values R_{-3}^{EE3}(-3i/2), R_{-1}^{EE3}(-3i/2) from kernels.json (sympy, exact).
 (2) continuum weights chi = Q B (c4v, c4d): analytic structure near +-3i/2 -- direct values, contour integrals
     (1/2 pi i) oint chi dnu around +-3i/2 (claim 0) and around +-i/2 (control: Q has a pole there).
 (3) the nu-normalisation 1/N_com, N_com ~ G(11/4 - i nu/2)^2 G(11/4 + i nu/2)^2: regular and non-zero at +-3i/2.
 (4) decay of F_{-1}(t) = int_R chi R_{-1} e^{-i nu t} dnu: a pole of chi at -3i/2 would give rate 3/2;
     for the continuum weights (quadrature) and for the Stage-7 XP grid weights (midpoint sum, nu_k = (k+1/2)/25);
     control: chi_ctrl = chi (nu^2+4)/(nu^2+9/4) (pole at +-3i/2 inserted).
 (5) bonus: the de-shift identity W[chi](t) = sum_j int_0^oo chi R_j cos((nu - i j)t) dnu = Phi(t) at a few t.
"""
import json, os, sys, time, math
import mpmath as mp
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'prl', 'data')
OUT = os.path.join(HERE, 'cache', 'A3_fpos.txt')
lines = []
def say(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); lines.append(s)

mp.mp.dps = 100

def load_R():
    kj = json.load(open(os.path.join(DATA, 'kernels.json')))
    nu = sp.Symbol('nu')
    R = {int(j): sp.sympify(e, locals={'nu': nu, 'I': sp.I}) for j, e in kj['Q_k4_EE3']['R_j'].items()}
    return nu, R

def load_weight(piece):
    d = json.load(open(os.path.join(DATA, f'q_deshift_continuum_{piece}.json')))
    h, sig = mp.mpf(d['h']), mp.mpf(d['sigma'])
    labels = d['labels'] if isinstance(d['labels'], list) else eval(d['labels'])
    coef = d['coefficients'] if isinstance(d['coefficients'], list) else eval(d['coefficients'])
    x = [(lab, mp.mpf(c)) for lab, c in zip(labels, coef)]
    p = 1 if piece == 'c4v' else 3
    s2 = 3 + 2 * mp.sqrt(2)
    def Q(v):
        return v**2 * (v**2 + 1) * (v**2 + 4) / ((v**2 + 9)**p * (4 * v**2 + 1) * (4 * v**2 + 25) * (4 * v**2 + s2**2))
    def g(z): return mp.exp(-z**2 / (2 * sig**2))
    def B(v):
        s = 0
        for (kind, k), c in x:
            if kind == 'g': s += c * g(v)
            elif kind == 're': s += c * (g(v - k * h) + g(v + k * h))
            else: s += c * 1j * (g(v - k * h) - g(v + k * h))
        return s
    return (lambda v: Q(v) * B(v)), Q, B

def main():
    t0 = time.time()
    nu, R = load_R()
    say('(1) exact values at nu = -3i/2 (sympy, from prl/data/kernels.json):')
    for j in (-3, -1, 1, 3):
        val = sp.nsimplify(sp.simplify(R[j].subs(nu, -sp.Rational(3, 2) * sp.I))) if True else None
        try:
            val = sp.simplify(R[j].subs(nu, -sp.Rational(3, 2) * sp.I))
        except Exception as e:
            val = f'pole ({e})'
        say(f'    R_{j:+d}(-3i/2) = {val}')
    say('    (R_{+3} has the factor (2 nu + 3i) in its denominator: limit -> ' +
        str(sp.limit(R[3] * (2 * nu + 3 * sp.I), nu, -sp.Rational(3, 2) * sp.I)) + ' x 1/(2nu+3i))')
    Rf = {j: sp.lambdify(nu, R[j], 'mpmath') for j in R}

    for piece in ('c4v', 'c4d'):
        chi, Q, B = load_weight(piece)
        say(f'\n=== weight {piece} ===')
        # (2) values and contour integrals
        for c0 in (mp.mpc(0, 1.5), mp.mpc(0, -1.5)):
            v = chi(c0)
            say(f'  chi({mp.nstr(c0, 3)}) = {mp.nstr(v, 15)}   (Q = {mp.nstr(Q(c0), 10)}, B = {mp.nstr(B(c0), 10)})')
        def oint(c0, rad=mp.mpf('0.1'), n=64):
            s = 0
            for k in range(n):
                th = 2 * mp.pi * (k + mp.mpf(1) / 2) / n
                z = c0 + rad * mp.exp(1j * th)
                s += chi(z) * 1j * rad * mp.exp(1j * th)
            return s * 2 * mp.pi / n / (2j * mp.pi)
        for c0, lab in ((mp.mpc(0, 1.5), '+3i/2'), (mp.mpc(0, -1.5), '-3i/2'), (mp.mpc(0, 0.5), '+i/2 (control)'),
                        (mp.mpc(0, 2.5), '+5i/2 (control)'), (mp.mpc(0, (3 + 2 * mp.sqrt(2)) / 2), '+2.914i (control)')):
            say(f'  (1/2 pi i) oint chi dnu around {lab:18s}: {mp.nstr(oint(c0), 6)}')
        # (3) normalisation
        Ncom = lambda v: mp.gamma(mp.mpf(11) / 4 - 1j * v / 2)**2 * mp.gamma(mp.mpf(11) / 4 + 1j * v / 2)**2
        say(f'  1/N_com(+-3i/2) = {mp.nstr(1 / Ncom(mp.mpc(0, 1.5)), 10)}, {mp.nstr(1 / Ncom(mp.mpc(0, -1.5)), 10)};'
            f'  1024/N_com = 2^21 pi^2 / N_com-check: n(3i/2) = {mp.nstr(1024 / Ncom(mp.mpc(0, 1.5)), 8)}')
        # position-space residue: Res_{3i/2} [w c_l]  with w = chi/N_com: w regular -> 0; show |chi| on circle bounded
        # (4) F_{-1}(t) decay -- Gauss-Legendre panels on [0, 70] (chi < 1e-25 of its max beyond 60)
        say('  F_{-1}(t) = int_R chi R_{-1} e^{-i nu t} dnu (chi even); local decay rate -dlog|F|/dt')
        xg, wg = mp.gauss_legendre(24) if hasattr(mp, 'gauss_legendre') else (None, None)
        if xg is None:
            import numpy as np
            xs, ws = np.polynomial.legendre.leggauss(24)
            # refine nodes to full precision by Newton on P_24
            xg, wg = [], []
            for x0 in xs:
                x = mp.mpf(x0)
                for _ in range(6):
                    p, dp = mp.legendre(24, x), mp.diff(lambda u: mp.legendre(24, u), x)
                    x = x - p / dp
                dp = mp.diff(lambda u: mp.legendre(24, u), x)
                xg.append(x); wg.append(2 / ((1 - x**2) * dp**2))
        nodes, wts = [], []
        for k in range(140):
            a, b_ = mp.mpf(k) / 2, mp.mpf(k + 1) / 2
            for x, w in zip(xg, wg):
                nodes.append((a + b_) / 2 + (b_ - a) / 2 * x); wts.append((b_ - a) / 2 * w)
        chis = [chi(v) for v in nodes]
        Rv = {j: [Rf[j](v) for v in nodes] for j in (-3, -1, 1, 3)}
        Rm1n = [Rf[-1](-v) for v in nodes]
        ts = [4, 6, 8, 10, 12]
        for lab, fac in (('weight', lambda v: 1), ('control (pole at +-3i/2 inserted)', lambda v: (v**2 + 4) / (v**2 + mp.mpf(9) / 4))):
            Fs = []
            facs = [fac(v) for v in nodes]
            for t in ts:
                s_ = 0
                for v, w, c, f, a, b2 in zip(nodes, wts, chis, facs, Rv[-1], Rm1n):
                    s_ += w * c * f * (a * mp.exp(-1j * v * t) + b2 * mp.exp(1j * v * t))
                Fs.append(s_)
            rates = [-(mp.log(abs(Fs[k + 1])) - mp.log(abs(Fs[k]))) / (ts[k + 1] - ts[k]) for k in range(len(ts) - 1)]
            say(f'    {lab:36s} |F(t)| at t={ts}: ' + ' '.join(mp.nstr(abs(F), 4) for F in Fs) +
                '; rates ' + ' '.join(mp.nstr(r_, 4) for r_ in rates))
        # (5) de-shift identity
        Phi = (lambda t: 1) if piece == 'c4v' else (lambda t: -mp.mpf(t)**2 / 2)
        for t in (0.5, 3, 8, 15):
            Wt = 0
            for k, (v, w, c) in enumerate(zip(nodes, wts, chis)):
                Wt += w * c * sum(Rv[j][k] * mp.cos((v - 1j * j) * t) for j in (-3, -1, 1, 3))
            say(f'  W(t={t}) = {mp.nstr(mp.re(Wt), 16)} (Im {mp.nstr(mp.im(Wt), 2)})   target {mp.nstr(Phi(t), 6)}   |W - Phi| = {mp.nstr(abs(Wt - Phi(t)), 3)}')
        say(f'  ({time.time()-t0:.0f} s)')

    # (4') Stage-7 XP grid weights: midpoint sums
    mp.mp.dps = 80
    say('\n=== Stage-7 XP grid weights (q_deshift_xp_C4_*.txt), F_{-1}(t) as the midpoint sum (1/25) sum_k [chi_k R_{-1}(nu_k) e^{-i nu_k t} + (nu -> -nu)]')
    for name in ('value', 'deriv'):
        rows = [l.split() for l in open(os.path.join(DATA, f'q_deshift_xp_C4_{name}.txt')) if not l.startswith('#')]
        nus = [mp.mpf(r_[1]) for r_ in rows]; chis = [mp.mpf(r_[2]) for r_ in rows]
        Rm = [Rf[-1](v) for v in nus]; Rp = [Rf[-1](-v) for v in nus]
        ts = [10, 15, 20, 25, 30, 36]
        Fs = []
        for t in ts:
            s = 0
            for v, c, a, b2 in zip(nus, chis, Rm, Rp):
                s += c * (a * mp.exp(-1j * v * t) + b2 * mp.exp(1j * v * t))
            Fs.append(s / 25)
        rates = [-(mp.log(abs(Fs[k + 1])) - mp.log(abs(Fs[k]))) / (ts[k + 1] - ts[k]) for k in range(len(ts) - 1)]
        R1 = complex(Rf[-1](mp.mpc(0, -1.5)))
        amp = [abs(F) * mp.exp(1.5 * t) / (2 * mp.pi * abs(R1)) for F, t in zip(Fs, ts)]
        say(f'  C4_{name}: |F_-1| at t={ts}: ' + ' '.join(mp.nstr(abs(F), 3) for F in Fs))
        say(f'      local rates ' + ' '.join(mp.nstr(r_, 4) for r_ in rates) +
            ';  implied residue |A(t)| = |F e^{3t/2}/(2 pi R_-1(-3i/2))|: ' + ' '.join(mp.nstr(a, 3) for a in amp))
    say(f'\n(elapsed {time.time()-t0:.0f} s)')
    open(OUT, 'w').write('\n'.join(lines) + '\n')

if __name__ == '__main__':
    main()
