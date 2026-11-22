"""
qf_tables (RESULT_Q_FWD_CLOSE.md): collects every qf_qan pickle into the document's tables:
  per (channel, state): s_d per gauge (H0 +- truncation, H2, H3, H4, SB, XP), FF/N (K = 0 -> measured) for delta = 0.05, 0.2 (H0, SB),
  covariant split, per-harmonic Gamma_K/sum|R_j| at nu = 0.5, 1, 1.5, 2, 2.5, 3, 4, 6, 10.
Usage: python3 qf_tables.py -> cache/qf_tables.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np

GS = ('H0', 'H2', 'H3', 'H4', 'SB', 'XP')
NUP = (0.5, 1.02, 1.5, 2.02, 2.5, 3.02, 4.02, 6.02, 10.02)

if __name__ == '__main__':
    import bf_qan as BQ
    rows = []
    for fn in sorted(glob.glob(os.path.join(HERE, 'cache', 'qf_qan_*_D*_J*_Dg*.pkl'))):
        if '_xs' in fn:
            continue
        d = pickle.load(open(fn, 'rb'))
        Lam = np.sqrt((d['Dg'] - 2.)*(d['Dg'] - 1.))
        for h in ('OO', '++', '+-'):
            if h not in d['res']:
                continue
            rows.append((h, d['J'], d['Dg'], d['m']/Lam, d['D'], d['K'], d['res'][h], d['rs'][d['K']][h]['GK'], d['mode']))
    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    nus = (np.arange(1025) + 0.5)*0.04
    R = BQ.Rj(nus[:1024]); AR = sum(np.abs(v) for v in R.values())
    idx = [int(np.argmin(np.abs(nus - v))) for v in NUP]
    L = ['# qf_tables', '| channel | J | Dg | m/Lam | Delta | K | s_d H0 | H2 | H3 | H4 | SB | XP | FF/N d=0.05 (K=0 -> meas, H0; SB) | d=0.2 | dgamma | src |',
         '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
    for h, J, Dg, mL, D, K, r, GK, mode in rows:
        f0, fS = r['H0']['ff'], r['SB']['ff']
        L.append(f'| {h} | {J} | {Dg} | {mL:.4f} | {D:.4f} | {K} | {r["H0"]["sd"]:+.3f} ± {r["H0"]["ed"]:.3f} | '
                 + ' | '.join(f'{r[g]["sd"]:+.3f}' for g in GS[1:]) + ' | '
                 + f'{f0[0.05][0]:+.4f} → {f0[0.05][1]:+.4f} ± {f0[0.05][2]*r["H0"]["ed"]:.4f} ({fS[0.05][1]:+.4f}) | '
                 + f'{f0[0.2][0]:+.4f} → {f0[0.2][1]:+.4f} ± {f0[0.2][2]*r["H0"]["ed"]:.4f} ({fS[0.2][1]:+.4f}) | {r["cov"]["dgamma"]:+.3f} | {mode} |')
    L.append('')
    L.append('# per-harmonic Gamma_K/sum|R_j| at nu = ' + ' '.join(f'{v:.2f}' for v in NUP))
    for h, J, Dg, mL, D, K, r, GK, mode in rows:
        L.append(f'  {h:2s} J={J:2d} Dg={Dg:4d} m={mL:.3f}Lam: ' + ' '.join(f'{(GK[i]/AR[i]).real:+6.2f}' for i in idx))
    print('\n'.join(L))
    open(os.path.join(HERE, 'cache', 'qf_tables.txt'), 'w').write('\n'.join(L) + '\n')
