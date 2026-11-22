"""repro_engine/summary.py -- final tables for REPRO_ENGINE.md: direct rho vs engine (K = 4, 5), s_d, FF/N, grid and precision checks,
order-0 check.   Usage: python3 summary.py -> cache/summary.txt
"""
import os
import numpy as np
import compare as CM
import sdchain as C
import direct_rho as DR

HERE = os.path.dirname(os.path.abspath(__file__))
A = '10-10-16-16-8-16-16'; B = '8-8-18-18-6-20-20'; Cg = '8-8-16-16-6-16-16'
NODES = (0.3, 0.5, 0.9, 1.3, 1.5, 1.7, 2.02, 2.5, 3.02, 3.5, 4.02, 5.02, 6.02, 8.02, 10.02)
out = []


def pr(s=''):
    print(s, flush=True); out.append(s)


def fexists(D, ns, prec, tag=''):
    return os.path.exists(os.path.join(HERE, 'cache', 'direct_D%g_n%s_%s%s.npz' % (D, ns, prec, tag)))


def gk_direct(D, ns, prec='f64', nsr=None, precr=None):
    nus, rho = CM.direct_ratio(D, ns, prec, nsr, precr)
    return nus, rho, C.gammaK(nus, rho, D)


def main():
    for D in (300, 1000):
        if not fexists(D, A, 'f64'):
            continue
        nus, rho, g = gk_direct(D, A)
        n = len(nus); sabs = g['sabs']
        e4n, e4r = C.engine_ratio(D, 4); g4 = C.gammaK(e4n, e4r, D)
        e5n, e5r = C.engine_ratio(D, 5); g5 = C.gammaK(e5n, e5r, D)
        hasB = fexists(D, B, 'f64')
        if hasB:
            _, rhoB, gB = gk_direct(D, B)
        pr(f'## Dg = {D}: direct grid A = {A} (B = {B} as convergence check), nu < {nus[-1] + 0.02:.2f}')
        pr('| nu | rho_dir/rho_dir(0.02) - rho_K5/rho_K5(0.02) (rel.) | Gamma_K/sum|R_j|: direct A | engine K=4 | engine K=5 | A - K5 | A - K4 | A - B |')
        pr('|---|---|---|---|---|---|---|---|')
        for v in NODES:
            i = int(round((v - 0.02)/0.04))
            if i >= n:
                continue
            dr = (rho[i]/rho[0])/(e5r[i]/e5r[0]) - 1
            a = g['GK'][i].real/sabs[i]; b4 = g4['GK'][i].real/sabs[i]; b5 = g5['GK'][i].real/sabs[i]
            ab = (g['GK'][i].real - gB['GK'][i].real)/sabs[i] if hasB else np.nan
            pr(f'| {nus[i]:.2f} | {abs(dr):.1e} | {a:+.4f} | {b4:+.4f} | {b5:+.4f} | {a - b5:+.4f} | {a - b4:+.4f} | {ab:+.4f} |')
        wins = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 6), (6, 8), (8, 10), (10, 12)]
        pr('')
        pr('| nu window | s_d direct (A) | engine K=4 | engine K=5 | direct - K5 | direct(A) - direct(B) | max abs(A - K5) per harmonic |')
        pr('|---|---|---|---|---|---|---|')
        for lo, hi in wins:
            sel = (nus >= lo) & (nus < hi)
            sa = C.sd(nus, g['GK'], lo=lo, hi=hi)[0]; s4 = C.sd(e4n, g4['GK'], lo=lo, hi=hi)[0]; s5 = C.sd(e5n, g5['GK'], lo=lo, hi=hi)[0]
            sb = C.sd(nus, gB['GK'], lo=lo, hi=hi)[0] if hasB else np.nan
            mx = np.max(np.abs((g['GK'][sel] - g5['GK'][:n][sel]).real/sabs[sel]))
            pr(f'| [{lo},{hi}) | {sa:+.5f} | {s4:+.5f} | {s5:+.5f} | {sa - s5:+.5f} | {sa - sb:+.1e} | {mx:.4f} |')
        numax = nus[-1] + 0.02
        sdA, svA = C.sd(nus, g['GK'], hi=numax)
        t5d, t5v = C.sd(e5n, g5['GK'], lo=numax); t4d, t4v = C.sd(e4n, g4['GK'], lo=numax)
        S4 = C.sd(e4n, g4['GK']); S5 = C.sd(e5n, g5['GK'])
        sd_tot = sdA + t5d; sv_tot = svA + t5v
        pr('')
        pr(f'* s_d (H0): direct(nu < {numax:.0f}) = {sdA:+.5f}; engine tail nu >= {numax:.0f}: K5 {t5d:+.5f} (K4 {t4d:+.5f}); '
           f'**total {sd_tot:+.5f}**;  engine K=4 {S4[0]:+.5f}, K=5 {S5[0]:+.5f}')
        if hasB:
            sdB = C.sd(nus, gB['GK'], hi=numax)[0] + t5d
            pr(f'* grid B total: {sdB:+.5f}  (A - B = {sd_tot - sdB:+.1e})')
        pr(f'* s_v: direct + tail = {sv_tot:+.4f}; engine K=4 {S4[1]:+.4f}, K=5 {S5[1]:+.4f}')
        for gauge in ('H2', 'H4'):
            sdg = C.sd(nus, g['GK'], gauge=gauge, hi=numax)[0] + C.sd(e5n, g5['GK'], gauge=gauge, lo=numax)[0]
            pr(f'* gauge {gauge}: direct + tail {sdg:+.5f}; engine K=5 {C.sd(e5n, g5["GK"], gauge=gauge)[0]:+.5f}')
        f_dir, H = C.ffn(D, sd_tot, sv_tot)
        f5, _ = C.ffn(D, S5[0], S5[1], H); f4, _ = C.ffn(D, S4[0], S4[1], H)
        sfail = H['FF']/(2*H['a'])
        pr(f'* FF/N (delta = 0.2, H0): direct {f_dir:+.5f}; engine K=4 {f4:+.5f}, K=5 {f5:+.5f};  s_fail = {sfail:.4f}, '
           f'margin s_fail - s_d = {sfail - sd_tot:.4f}')
        pr('')
    # precision check: f80 vs f64 on grid C (Dg = 300, nu < 4.1)
    if fexists(300, Cg, 'f80') and fexists(20000, Cg, 'f80'):
        nus, r64, g64 = gk_direct(300, Cg, 'f64')
        _, r80, g80 = gk_direct(300, Cg, 'f80')
        k = min(len(nus), len(g80['GK']))
        d = np.max(np.abs((g64['GK'][:k] - g80['GK'][:k]).real/g64['sabs'][:k]))
        s64 = C.sd(nus[:k], g64['GK'][:k], hi=4.1)[0]; s80 = C.sd(nus[:k], g80['GK'][:k], hi=4.1)[0]
        pr(f'* precision (grid C, Dg = 300, nu < 4.1): max |Gamma_K(f64) - Gamma_K(f80)| = {d:.1e} c2 per harmonic; '
           f's_d(nu<4.1): {s64:+.6f} (f64) vs {s80:+.6f} (f80)')
    # order-0 check: state <TTO> coefficients at 300, geometry at Delta_g = 1e6, 2e6 (Richardson in a), reference = direct 2e4 (grid C, f80)
    o = {}
    for Dgeo in (1e6, 2e6):
        if fexists(Dgeo, Cg, 'f80', '_coef300'):
            nu_, I, _ = DR.load(os.path.join(HERE, 'cache', 'direct_D%g_n%s_f80_coef300.npz' % (Dgeo, Cg)))
            o[Dgeo] = (nu_, I)
    if len(o) == 2 and fexists(20000, Cg, 'f80'):
        nu_, Iref, _ = DR.load(os.path.join(HERE, 'cache', 'direct_D20000_n%s_f80.npz' % Cg))
        I1, I2 = o[1e6][1], o[2e6][1]
        I0 = 2*I2/I2[0]*1.0 - I1/I1[0]                      # Richardson a -> 0 (normalised at nu = 0.02)
        g0 = C.gammaK(nu_, I0/Iref, 300.0)
        e0n, e0r = C.engine_ratio(300, 0); ge0 = C.gammaK(e0n, e0r, 300.0)
        k = len(nu_)
        d = np.max(np.abs((g0['GK'] - ge0['GK'][:k]).real/g0['sabs']))
        s0d = C.sd(nu_, g0['GK'], hi=4.1)[0]; s0e = C.sd(e0n, ge0['GK'], hi=4.1)[0]
        pr(f'* order 0 (state a -> 0 at fixed <TTO>(300) coefficients; Richardson 1e6/2e6): max |Gamma_K - engine K=0| on nu < 4.1 = '
           f'{d:.3f} c2/harmonic; s_d(nu<4.1): direct {s0d:+.4f} vs engine K=0 {s0e:+.4f} (engine K=0 total {C.sd(e0n, ge0["GK"])[0]:+.4f})')
    open(os.path.join(HERE, 'cache', 'summary.txt'), 'w').write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
