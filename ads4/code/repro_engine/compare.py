"""repro_engine/compare.py -- direct rho(nu) (direct_rho.py) vs the engine's rho (K = 0..5), per harmonic and in s_d.
Usage: python3 compare.py Delta ns_state[;ns_ref] prec [prec_ref]  -> stdout (+ appended to cache/compare.txt)
"""
import os, sys
import numpy as np
import direct_rho as DR
import sdchain as C

HERE = os.path.dirname(os.path.abspath(__file__))
SHOW = (0.02, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 2.02, 2.5, 3.02, 3.5, 4.02, 5.02, 6.02, 8.02, 10.02, 11.5)


def direct_ratio(D, ns, prec, nsr=None, precr=None):
    f = lambda d, n, p: os.path.join(HERE, 'cache', 'direct_D%g_n%s_%s.npz' % (d, n, p))
    nus, Ist, _ = DR.load(f(D, ns, prec))
    _, Iref, _ = DR.load(f(2e4, nsr or ns, precr or prec))
    return nus, Ist/Iref


def main(D, ns, prec, nsr=None, precr=None, quiet=False):
    nus, rho = direct_ratio(D, ns, prec, nsr, precr)
    nmax = len(nus)
    out = []
    pr = lambda s='': (out.append(s), None if quiet else print(s, flush=True))
    gd = C.gammaK(nus, rho, D)
    eng = {}
    for K in (0, 3, 4, 5):
        en, er = C.engine_ratio(D, K)
        ge = C.gammaK(en, er, D)
        eng[K] = (en, ge)
    sabs = gd['sabs']
    pr(f'# compare D={D:g}: direct n={ns} ({prec}; ref n={nsr or ns} {precr or prec}), nu < {nus[-1] + 0.02:.2f}')
    pr('  nu    | Gamma_K/sum|R_j|: direct   engine K=4   K=5   | direct - K5 (c2 per harmonic) | rho_dir/rho_K5 normalised - 1')
    r5 = eng[5][1]
    en5, er5 = C.engine_ratio(D, 5)
    for v in SHOW:
        i = int(round((v - 0.02)/0.04))
        if i >= nmax:
            continue
        a = gd['GK'][i].real/sabs[i]; b = eng[4][1]['GK'][i].real/sabs[i]; c = r5['GK'][i].real/sabs[i]
        dr = (rho[i]/rho[0])/(er5[i]/er5[0]) - 1
        pr(f'  {nus[i]:5.2f} | {a:+9.3f}  {b:+9.3f}  {c:+9.3f} | {a - c:+8.3f} | {dr.real:+.2e} {dr.imag:+.2e}')
    wins = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 6), (6, 8), (8, 10), (10, 12)]
    pr('  s_d (H0) by window: direct | engine K=4 | K=5 | direct-K5;  max |direct-K5| per harmonic in window')
    for lo, hi in wins:
        if lo >= nus[-1]:
            continue
        sel = (nus >= lo) & (nus < hi)
        sdd = C.sd(nus, gd['GK'], lo=lo, hi=hi)[0]
        s4 = C.sd(eng[4][0], eng[4][1]['GK'], lo=lo, hi=hi)[0]; s5 = C.sd(eng[5][0], r5['GK'], lo=lo, hi=hi)[0]
        mx = np.max(np.abs((gd['GK'][sel] - r5['GK'][:nmax][sel]).real/sabs[sel]))
        pr(f'     [{lo:2d},{hi:2d}): {sdd:+.5f} | {s4:+.5f} | {s5:+.5f} | {sdd - s5:+.5f};  {mx:.3f}')
    # spliced totals: direct below numax, engine (K = 5) above
    numax = nus[-1] + 0.02
    sdd = C.sd(nus, gd['GK'], hi=numax); svd = C.sd(nus, gd['GK'], hi=numax)[1]
    tail5 = C.sd(eng[5][0], r5['GK'], lo=numax); tail4 = C.sd(eng[4][0], eng[4][1]['GK'], lo=numax)
    tot = {}
    for K in (0, 3, 4, 5):
        tot[K] = C.sd(eng[K][0], eng[K][1]['GK'])
    pr(f'  s_d total: direct(nu<{numax:.2f}) + K5 tail = {sdd[0] + tail5[0]:+.5f}  (tail {tail5[0]:+.5f});  engine K=0/3/4/5: '
       + ' / '.join(f'{tot[K][0]:+.5f}' for K in (0, 3, 4, 5)))
    pr(f'  s_v total: direct + K5 tail = {sdd[1] + tail5[1]:+.4f};  engine K=4/5: {tot[4][1]:+.4f} / {tot[5][1]:+.4f}')
    open(os.path.join(HERE, 'cache', 'compare.txt'), 'a').write('\n'.join(out) + '\n')
    return dict(nus=nus, gd=gd, eng=eng, sd_dir=sdd[0] + tail5[0], sv_dir=sdd[1] + tail5[1], tot=tot, sd_low=sdd[0], tail=tail5[0])


if __name__ == '__main__':
    D = float(sys.argv[1]); ns = sys.argv[2].split(';')
    prec = sys.argv[3] if len(sys.argv) > 3 else 'f64'
    precr = sys.argv[4] if len(sys.argv) > 4 else None
    main(D, ns[0].replace(',', '-'), prec, ns[1].replace(',', '-') if len(ns) > 1 else None, precr)
