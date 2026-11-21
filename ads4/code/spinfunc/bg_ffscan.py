"""
bg_ffscan (BELOW_FOLD2.md §4.3): min FF/N over Stage 6's state set (stage6_stress.states, HYB/EXA as Stage 6/7) when the MEASURED
derivative projection s_d of the analytic-radial EE k = 4 image (bg_qan, per gauge) is applied to the '++' states of the measured
spins (J = 0, 2; s_d taken m-independent: measured at m = Lam and 1.05 Lam, Delta_g = 300 and 1000, spread quoted), all other states
at K = 0.  FF -> FF - a 2u^3 s_d  (the value piece, Lam^-2-suppressed, is dropped).
Usage: python3 bg_ffscan.py   (reads cache/bg_qan_D*_J*.pkl)  -> cache/bg_ffscan.txt
"""
import os, sys, pickle, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import numpy as np

if __name__ == '__main__':
    import stage6_model as M, stage6_stress as ST
    from stage6_scan import XS, JCAP, BMAX
    sd = {}
    for fn in glob.glob(os.path.join(HERE, 'cache', 'bg_qan_D*_J*.pkl')):
        J = int(fn.split('_J')[1].split('.')[0])
        P = pickle.load(open(fn, 'rb'))['P']; K = max(P)
        for g, (s, _) in P[K].items():
            sd.setdefault((J, g), []).append((s, abs(s - P[K - 1][g][0])))
    lines = ['# bg_ffscan: min FF/N with measured s_d on ++ J in ' + str(sorted(set(k[0] for k in sd)))]
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
            for g in ('H0', 'H2', 'H3', 'H4', 'SB', 'XP'):
                dF = np.zeros(len(rows)); eF = np.zeros(len(rows))
                for J in sorted(set(k[0] for k in sd)):
                    if (J, g) not in sd:
                        continue
                    vals = np.array(sd[(J, g)])
                    s, e = vals[:, 0].mean(), max(vals[:, 1].max(), np.ptp(vals[:, 0])/2)
                    sel = (hs == '++') & (Js == J)
                    dF[sel] = -a*2*U[sel]**3*s; eF[sel] = a*2*U[sel]**3*e
                r = (FF + dF)/NN; i = int(np.argmin(r))
                lines.append(f'  delta={delta} Dg={Dg} {g}: K=0 min {base:+.4f};  measured: min FF/N = {r[i]:+.4f} +- {eF[i]/NN[i]:.4f} '
                             f'at {lab[i][0]} J={lab[i][1]} m={lab[i][2]:.3f}Lam')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'bg_ffscan.txt'), 'w').write('\n'.join(lines) + '\n')
