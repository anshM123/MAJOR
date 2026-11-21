"""
v4d_strip (RESULT_Q_V4DPRIME.md §3): is the C4 de-shift weight analytic in the strip |Im nu| <= 3 (hypothesis 5)?

Test by Paley-Wiener: for the even extension chi(|nu|), and for g_j(nu) := chi(|nu|) R_j(nu) (full-line form of the Abel operator,
g_j(-nu) = chi R_{-j}(nu)), the Fourier transforms
    chih(t) = int chi(|nu|) e^{i nu t} dnu,       gh_j(t) = int g_j(nu) e^{i nu t} dnu = C[chi A_j](t) - S[chi B_j](t),
decay like e^{-w t} iff the function is analytic in 0 < Im nu < w (poles on Im nu = w give polynomial prefactors).  The Abel equation reads
    W(t) = (1/2) sum_{j=1,3} [e^{j t} gh_j(t) + e^{-j t} gh_j(-t)] = Phi(t).
Uses the Stage-7 XP weights at their full 64-digit precision (the only weights realised to t = 30); midpoint sums are exact
Fourier transforms of the grid function up to aliasing at t ~ 2 pi/dnu = 157.
Usage: python3 v4d_strip.py -> cache/v4d_strip.txt  (~2 min)
"""
import os, sys, time, json
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '2')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import mpmath as mp
import sympy as sp
import numpy as np
import v4d_core as D

CACHE = os.path.join(HERE, 'cache')
L = []


def pr(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); L.append(s)


def main():
    t0 = time.time()
    mp.mp.dps = 80
    nus = sp.Symbol('nu')
    Rj = json.load(open(os.path.join(HERE, 'functionals_q4_deshift.json')))['EE3']['R_j']
    Rf = {int(j): sp.lambdify(nus, sp.sympify(v, locals={'nu': nus, 'I': sp.I}), 'mpmath') for j, v in Rj.items()}
    pr('# v4d_strip: Paley-Wiener test of the strip analyticity of the C4 de-shift weights (XP, 64 digits)')
    for name, tgt in (('xp_value', lambda t: mp.mpf(1)), ('xp_deriv', lambda t: -t**2/2)):
        nu, chi, dnu = D.weights_mp(name, dps=80)
        R = {j: [Rf[j](v) for v in nu] for j in Rf}
        TS = [mp.mpf(t) for t in (0, 1, 2, 3, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30)]
        pr(f'## {name}: N = {len(nu)}, dnu = {mp.nstr(dnu, 6)}')
        pr('    t   |chih(t)|      rate_chi   |gh_3(t)|     rate_3   |gh_1(t)|     rate_1   e^{3t}|gh_3|   e^{t}|gh_1|   |W/Phi - 1|')
        prev = None
        rows = []
        for t in TS:
            c = [mp.cos(v*t) for v in nu]; s = [mp.sin(v*t) for v in nu]
            chih = 2*dnu*mp.fsum(chi[k]*c[k] for k in range(len(nu)))
            gh = {}
            for j in (1, 3):
                gp = dnu*mp.fsum(chi[k]*(R[j][k]*mp.mpc(c[k], s[k]) + R[-j][k]*mp.mpc(c[k], -s[k])) for k in range(len(nu)))
                gm = dnu*mp.fsum(chi[k]*(R[j][k]*mp.mpc(c[k], -s[k]) + R[-j][k]*mp.mpc(c[k], s[k])) for k in range(len(nu)))
                gh[j] = (gp, gm)
            W = (mp.exp(3*t)*gh[3][0] + mp.exp(-3*t)*gh[3][1] + mp.exp(t)*gh[1][0] + mp.exp(-t)*gh[1][1])/2
            ph = tgt(t)
            err = abs(W/ph - 1) if ph != 0 else abs(W)
            rows.append((t, chih, gh[3][0], gh[1][0], err))
            if prev is not None:
                dt = t - prev[0]
                rc = -mp.log(abs(chih)/abs(prev[1]))/dt; r3 = -mp.log(abs(gh[3][0])/abs(prev[2]))/dt; r1 = -mp.log(abs(gh[1][0])/abs(prev[3]))/dt
            else:
                rc = r3 = r1 = mp.nan
            pr(f'   {int(t):2d}  {mp.nstr(abs(chih), 4):>10s}  {mp.nstr(rc, 4):>8s}  {mp.nstr(abs(gh[3][0]), 4):>10s}  {mp.nstr(r3, 4):>8s}  '
               f'{mp.nstr(abs(gh[1][0]), 4):>10s}  {mp.nstr(r1, 4):>8s}  {mp.nstr(mp.exp(3*t)*abs(gh[3][0]), 4):>11s}  '
               f'{mp.nstr(mp.exp(t)*abs(gh[1][0]), 4):>11s}  {mp.nstr(err, 3):>9s}')
            prev = (t, chih, gh[3][0], gh[1][0])
        # cancellation structure: e^{3t} gh_3 + e^{t} gh_1 (the growing parts) vs Phi
        pr('    growing parts: t, e^{3t} gh_3(t), e^{t} gh_1(t), their sum/2 (must equal Phi up to decaying terms):')
        for (t, chih, g3, g1, err) in rows[::3]:
            a3 = mp.exp(3*t)*g3; a1 = mp.exp(t)*g1
            pr(f'     t = {int(t):2d}: {mp.nstr(a3, 6)}   {mp.nstr(a1, 6)}   {mp.nstr((a3 + a1)/2, 8)}')
    pr(f'[{time.time() - t0:.0f}s]')
    open(os.path.join(CACHE, 'v4d_strip.txt'), 'w').write('\n'.join(L) + '\n')


if __name__ == '__main__':
    main()
