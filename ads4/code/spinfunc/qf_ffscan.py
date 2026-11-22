"""
qf_ffscan (RESULT_Q_FWD_CLOSE.md sect. 4): copy of po_ffscan (PARITY_ODD_CLOSE §2.4; itself bh_ffscan + J = 5) with
  * the qf_qan measurements added: OO J = 0, 2 (parity-odd even J), +- J = 6, 7, 8, 9 (and 11 if present), ++ J = 6, 8,
    ++ J = 0, 2 at 1.3 Lam;
  * an explicit OO channel: Stage 6's state set has no OO rows (OO was identified with ++); here every even-J ++ row is duplicated as
    an OO row (same HYB K = 0 value, as the model prescribes), carrying the MEASURED OO s_d where measured (J = 0, 2), else the ++ s_d
    of that J where measured (lemma extension, flagged 'L'), else K = 0.
Measured s_d applied m-independently within each (h, J) [M]; mean over the masses/Delta_g measured; error = max(truncation, half-spread).
Also: the universal-sigma_h extrapolation (s_d = sigma_h s_0 on every unmeasured state, sigma_OO = +1, odd-J +- with sigma = -1 now
measured) and per-channel minima with the failure threshold.
Usage: python3 qf_ffscan.py -> cache/qf_ffscan.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

GS = ('H0', 'H2', 'H3', 'H4', 'SB', 'XP')


def measured():
    sd = {}

    def add(h, J, Dg, s, e, g):
        sd.setdefault((h, J), {}).setdefault(Dg, []).append((s, e, g))
    for fn in glob.glob(os.path.join(HERE, 'cache', 'bg_qan_D*_J*.pkl')):
        J = int(fn.split('_J')[1].split('.')[0])
        P = pickle.load(open(fn, 'rb'))['P']; K = max(P)
        for g, (s, _) in P[K].items():
            add('++', J, None, s, abs(s - P[K - 1][g][0]), g)
    for fn in glob.glob(os.path.join(HERE, 'cache', 'bh_qan4_D*_J*_Dg*.pkl')) + glob.glob(os.path.join(HERE, 'cache', 'po_qan5_D*_J*_Dg*.pkl')):
        d = pickle.load(open(fn, 'rb'))
        for h in ('++', '+-'):
            if h in d['res']:
                for g in GS:
                    r = d['res'][h][g]; add(h, d['J'], d['Dg'], r['sd'], r['ed'], g)
    for fn in glob.glob(os.path.join(HERE, 'cache', 'qf_qan_*_D*_J*_Dg*.pkl')):
        if '_xs' in fn:
            continue
        d = pickle.load(open(fn, 'rb'))
        for h in ('++', '+-', 'OO'):
            if h in d['res']:
                if h == '++' and d['mode'] == 'oo' and d['J'] in (0, 2) and abs(d['m']/np.sqrt((d['Dg'] - 2.)*(d['Dg'] - 1.)) - 1.3) > 0.05:
                    continue            # ++ at the bg_eeq states is already in via bg_qan pickles (same images, REFCORR differs by 0.005 x gamma)
                for g in GS:
                    r = d['res'][h][g]; add(h, d['J'], d['Dg'], r['sd'], r['ed'], g)
    return sd


if __name__ == '__main__':
    import stage6_model as M, stage6_stress as ST
    from stage6_scan import XS, JCAP, BMAX
    SD = measured()
    lines = ['# qf_ffscan: min FF/N with measured s_d on ' + ', '.join(f'{h} J={J} (Dg {sorted(str(k) for k in v)}, n={sum(len(x) for x in v.values())//6})'
                                                                    for (h, J), v in sorted(SD.items()))]
    summary = {}
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a, c = F['a'], F['c']
        flat = M.flat_table(name, XS, JCAP)
        for Dg in (300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = ST.states(name, Dg, flat, Lam)
            rows = rows + [('OO',) + tuple(r[1:]) for r in rows if r[0] == '++']
            FF = np.zeros(len(rows)); NN = np.zeros(len(rows)); U = np.zeros(len(rows)); lab = []
            Js = np.zeros(len(rows), int); hs = []
            for i, (h, u, J, b2, b3, fw) in enumerate(rows):
                f, n, _, _ = M.hyb(F, '++' if h == 'OO' else h, u, float(J), b2, b3, tab, Lam, fw=fw)
                FF[i], NN[i], U[i], Js[i] = f, n, u, J; hs.append(h); lab.append((h, J, 1/np.sqrt(u)))
            hs = np.array(hs)
            base = np.min(FF/NN)
            for g in GS:
                dF = np.zeros(len(rows)); eF = np.zeros(len(rows)); used = []; meas = np.zeros(len(rows), bool)
                chans = sorted(set(SD) | {('OO', J) for (h, J) in SD if h == '++' and J % 2 == 0})
                for (h, J) in chans:
                    src, flag = (h, J), ''
                    if src not in SD:
                        src, flag = ('++', J), 'L'
                    byDg = SD[src]
                    vals = byDg.get(Dg, []) + byDg.get(None, [])
                    if not vals:
                        vals = [v for k in byDg for v in byDg[k]]; flag += '*'
                    vals = np.array([(s, e) for (s, e, gg) in vals if gg == g])
                    if not len(vals):
                        continue
                    s, e = vals[:, 0].mean(), max(vals[:, 1].max(), np.ptp(vals[:, 0])/2)
                    sel = (hs == h) & (Js == J)
                    dF[sel] = -a*2*U[sel]**3*s; eF[sel] = a*2*U[sel]**3*e; meas[sel] = True
                    used.append(f'{h}{J}{flag}:{s:+.3f}')
                r = (FF + dF)/NN; i = int(np.argmin(r))
                per = []
                for (h, J) in chans:
                    sel = (hs == h) & (Js == J)
                    if sel.any():
                        k = np.where(sel)[0][np.argmin(r[sel])]
                        per.append(f'{h}J{J} {r[k]:+.4f}+-{eF[k]/NN[k]:.4f}@{lab[k][2]:.3f}')
                # unmeasured channels: universal sigma_h s_0 (s_0 = mean ++ measurement; sigma_++ = sigma_OO = +1, sigma_+- = -1)
                s0 = np.mean([v[0] for (h_, J_), byDg in SD.items() if h_ == '++' for k_ in byDg for v in byDg[k_] if v[2] == g])
                sig = np.where(hs == '+-', -1.0, 1.0)
                ru = np.where(meas, r, (FF - a*2*U**3*s0*sig)/NN)
                iu = int(np.argmin(ru))
                un = ~meas
                iun = np.where(un)[0][np.argmin(ru[un])] if un.any() else -1
                # worst sign on all unmeasured states (|s_d| = 3 s_0, either sign): how far from failure
                rw = np.where(meas, np.inf, (FF - a*2*U**3*3*s0)/NN)
                iw = int(np.argmin(rw))
                lines.append(f'  delta={delta} Dg={Dg} {g}: K=0 min {base:+.4f}; measured min FF/N = {r[i]:+.4f} +- {eF[i]/NN[i]:.4f} '
                             f'at {lab[i][0]} J={lab[i][1]} m={lab[i][2]:.3f}Lam')
                lines.append(f'      per channel: ' + '; '.join(per))
                lines.append(f'      s_d used: ' + ' '.join(used))
                lines.append(f'      + unmeasured states at sigma_h s_0 (s_0 = {s0:+.3f}): min {ru[iu]:+.4f} at {lab[iu][0]} J={lab[iu][1]} m={lab[iu][2]:.3f}; '
                             f'lowest unmeasured state {ru[iun]:+.4f} at {lab[iun][0]} J={lab[iun][1]} m={lab[iun][2]:.3f}; '
                             f'unmeasured states alone at s_d = +3 s_0 (adverse sign, 3x size): min {rw[iw]:+.4f} at {lab[iw][0]} J={lab[iw][1]} m={lab[iw][2]:.3f}')
                summary[(delta, Dg, g)] = (r[i], eF[i]/NN[i], lab[i])
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'qf_ffscan.txt'), 'w').write('\n'.join(lines) + '\n')
    pickle.dump(summary, open(os.path.join(HERE, 'cache', 'qf_ffscan.pkl'), 'wb'))
