"""
qf_tol (RESULT_Q_FWD_CLOSE.md sect. 1): per-channel tolerances of the Result Q heavy set (Stage-6 HYB state set, stage6_stress.states;
logic copied from po_ffscan) for the forward-derivative correction s_d: for every (h, J), J <= 13, at K = 0:
  min FF/N and its m/Lam;  the failure threshold s_fail = min over the channel's states of FF/(a 2u^3) (s_d needed for FF = 0);
  the masses where FF/N - a 2u^3 s/N is smallest for s = +-0.51 (the measured universal size).
Usage: python3 qf_tol.py -> cache/qf_tol.txt
"""
import os, sys
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
RQ = os.path.join(HERE, '..', 'resultq'); sys.path.insert(0, RQ)
import warnings
warnings.filterwarnings('ignore')
import numpy as np

if __name__ == '__main__':
    import stage6_model as M, stage6_stress as ST
    from stage6_scan import XS, JCAP, BMAX
    lines = ['# qf_tol: per-channel K = 0 minima, s_d failure thresholds, minima at s_d = +-0.51 (HYB state set)']
    for delta in (0.05, 0.2):
        name = M.FUNCS[delta]; F, c2, c3 = M.load(name); a = F['a']
        flat = M.flat_table(name, XS, JCAP)
        for Dg in (300, 1000):
            Lam = M.lam(Dg); tab = M.Tables(name, Lam, BMAX[Dg])
            rows = ST.states(name, Dg, flat, Lam)
            FF = np.zeros(len(rows)); NN = np.zeros(len(rows)); U = np.zeros(len(rows)); Js = np.zeros(len(rows), int); hs = []
            for i, (h, u, J, b2, b3, fw) in enumerate(rows):
                f, n, _, _ = M.hyb(F, h, u, float(J), b2, b3, tab, Lam, fw=fw)
                FF[i], NN[i], U[i], Js[i] = f, n, u, J; hs.append(h)
            hs = np.array(hs); mm = 1/np.sqrt(U); co = a*2*U**3/NN
            for h in ('++', '+-'):
                for J in range(0, 14):
                    sel = (hs == h) & (Js == J)
                    if not sel.any():
                        continue
                    ix = np.where(sel)[0]
                    r0 = FF[ix]/NN[ix]; k0 = ix[np.argmin(r0)]
                    sfail = FF[ix]/(a*2*U[ix]**3); kf = ix[np.argmin(sfail)]
                    rp = r0 - co[ix]*0.51; kp = ix[np.argmin(rp)]
                    rm = r0 + co[ix]*0.51; km = ix[np.argmin(rm)]
                    near = ix[mm[ix] < 1.6]
                    lines.append(f'  d={delta} Dg={Dg} {h} J={J:2d}: n={sel.sum():4d}  K0 min {FF[k0]/NN[k0]:+.4f}@{mm[k0]:.3f}  coef a2u^3/N there {co[k0]:.4f} '
                                 f'(max {co[ix].max():.4f})  s_fail(min) {sfail.min():+8.2f}@{mm[kf]:.3f}  | s=+0.51: {rp.min():+.4f}@{mm[kp]:.3f}  '
                                 f's=-0.51: {rm.min():+.4f}@{mm[km]:.3f} | gap m {mm[ix].min():.4f}, states m<1.6: {len(near)}')
    print('\n'.join(lines))
    open(os.path.join(HERE, 'cache', 'qf_tol.txt'), 'w').write('\n'.join(lines) + '\n')
