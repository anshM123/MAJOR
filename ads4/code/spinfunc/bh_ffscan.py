"""
bh_ffscan (BELOW_FOLD3.md §1.5): bg_ffscan (copied logic) with the +- channel added: min FF/N over Stage 6's state set
(stage6_stress.states, HYB/EXA) when the MEASURED derivative projections s_d (per gauge) are applied to
  '++' J = 0, 2  (bg_qan pickles, BELOW_FOLD2),  '++' J = 4 and '+-' J = 4 (bh_qan4 pickles, this note; measured at the gap state),
each taken m-independent within its (h, J) [M] (as bg_ffscan); every other state at K = 0.  FF -> FF - a 2u^3 s_d (value piece dropped,
Lam^-2-suppressed).  The s_d of the Delta_g at hand is used where measured, otherwise the other Delta_g's (flagged).
Error = max(truncation, spread over the measurements of that (h, J)).
Usage: python3 bh_ffscan.py  -> cache/bh_ffscan.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np

GS = ('H0', 'H2', 'H3', 'H4', 'SB', 'XP')


def measured():
    sd = {}
    for fn in glob.glob(os.path.join(HERE, 'cache', 'bg_qan_D*_J*.pkl')):
        J = int(fn.split('_J')[1].split('.')[0])
        P = pickle.load(open(fn, 'rb'))['P']; K = max(P)
        for g, (s, _) in P[K].items():
            sd.setdefault(('++', J), {}).setdefault(None, []).append((s, abs(s - P[K - 1][g][0]), g))
    for fn in glob.glob(os.path.join(HERE, 'cache', 'bh_qan4_D*_J*_Dg*.pkl')):
        d = pickle.load(open(fn, 'rb'))
        for h in ('++', '+-'):
            if h not in d['res']:
                continue
            for g in GS:
                r = d['res'][h][g]
                sd.setdefault((h, d['J']), {}).setdefault(d['Dg'], []).append((r['sd'], r['ed'], g))
    return sd


if __name__ == '__main__':
    import stage6_model as M, stage6_stress as ST
    from stage6_scan import XS, JCAP, BMAX
    SD = measured()
    lines = ['# bh_ffscan: min FF/N with measured s_d on ' + ', '.join(f'{h} J={J} (Dg {sorted(str(k) for k in v)})' for (h, J), v in sorted(SD.items()))]
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        flat = M.flat_table(name, XS, JCAP)
        for Dg in (300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = ST.states(name, Dg, flat, Lam)
            FF = np.zeros(len(rows)); NN = np.zeros(len(rows)); U = np.zeros(len(rows)); lab = []
            Js = np.zeros(len(rows), int); hs = []
            for i, (h, u, J, b2, b3, fw) in enumerate(rows):
                f, n, _, _ = M.hyb(F, h, u, float(J), b2, b3, tab, Lam, fw=fw)
                FF[i], NN[i], U[i], Js[i] = f, n, u, J; hs.append(h); lab.append((h, J, 1/np.sqrt(u)))
            hs = np.array(hs)
            base = np.min(FF/NN)
            for g in GS:
                dF = np.zeros(len(rows)); eF = np.zeros(len(rows)); used = []
                for (h, J), byDg in SD.items():
                    vals = byDg.get(Dg) or byDg.get(None) or [v for k in byDg for v in byDg[k]]
                    flag = '' if (Dg in byDg or None in byDg) else '*'
                    vals = np.array([(s, e) for (s, e, gg) in vals if gg == g])
                    if not len(vals):
                        continue
                    s, e = vals[:, 0].mean(), max(vals[:, 1].max(), np.ptp(vals[:, 0])/2)
                    sel = (hs == h) & (Js == J)
                    dF[sel] = -a*2*U[sel]**3*s; eF[sel] = a*2*U[sel]**3*e
                    used.append(f'{h}{J}{flag}:{s:+.3f}')
                r = (FF + dF)/NN; i = int(np.argmin(r))
                # minimum restricted to each measured channel
                per = []
                for (h, J) in sorted(SD):
                    sel = (hs == h) & (Js == J)
                    if sel.any():
                        k = np.where(sel)[0][np.argmin(r[sel])]
                        per.append(f'{h}J{J} {r[k]:+.4f}+-{eF[k]/NN[k]:.4f}@{lab[k][2]:.3f}')
                # universal-sigma_h extrapolation [M]: s_d = sigma_h s_0 on EVERY state (s_0 = mean ++ measurement of this gauge);
                # odd-J +- (parity-odd structure, not measured) bracketed with both signs (worst kept)
                s0 = np.mean([v[0] for (h_, J_), byDg in SD.items() if h_ == '++' for k_ in byDg for v in byDg[k_] if v[2] == g])
                odd = (hs == '+-') & (Js % 2 == 1)
                sig = np.where(hs == '++', 1.0, -1.0)
                ru = (FF - a*2*U**3*s0*sig)/NN
                ru = np.where(odd, np.minimum(ru, (FF - a*2*U**3*s0)/NN), ru)
                iu = int(np.argmin(ru)); ieven = int(np.argmin(np.where(odd, np.inf, ru)))
                lines.append(f'  delta={delta} Dg={Dg} {g}: universal sigma_h s_0 = {s0:+.3f} on all states: min {ru[iu]:+.4f} at {lab[iu][0]} J={lab[iu][1]} '
                             f'm={lab[iu][2]:.3f}; excluding odd-J +-: {ru[ieven]:+.4f} at {lab[ieven][0]} J={lab[ieven][1]} m={lab[ieven][2]:.3f}; '
                             f'odd-J +- worst sign: {np.min(ru[odd]) if odd.any() else np.nan:+.4f}')
                lines.append(f'  delta={delta} Dg={Dg} {g}: K=0 min {base:+.4f}; measured min FF/N = {r[i]:+.4f} +- {eF[i]/NN[i]:.4f} '
                             f'at {lab[i][0]} J={lab[i][1]} m={lab[i][2]:.3f}Lam | per channel: ' + '; '.join(per) + ' | s_d used ' + ' '.join(used))
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bh_ffscan.txt'), 'w').write('\n'.join(lines) + '\n')
