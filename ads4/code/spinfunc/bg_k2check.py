"""
bg_k2check (BELOW_FOLD2.md §3.3 (b)): EE k = 2 state/ref of the analytic-radial expansion vs the K2 Laplace rows above the fold
(cache/w1_pts_EE_J{J}_k2_b.pkl: Delta = 400, 1600, ref 2e4; nu = 19.92, 39.85, ...; Laplace orders M0, M2, M4, M6).
Compared: m^2 [D(nu) - D(nu_1)], D = log(I_st/I_ref) (Laplace: log sum_k M_k + logC - logK), nu_1 = 19.92.
Usage: python3 bg_k2check.py J  -> cache/bg_k2check.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import bg_eeq as EQ

if __name__ == '__main__':
    J = int(sys.argv[1])
    d = pickle.load(open(os.path.join(HERE, 'cache', f'w1_pts_EE_J{J}_k2_b.pkl'), 'rb'))
    nus = sorted(set(k[1] for k in d if k[0] == 400.0 and k[1] < 45))
    lines = []
    ref = EQ.images(EQ.DREF, J, 2, 3, nus=np.array(nus), kinds=('EE2',), verbose=False)
    for D in [D_ for D_ in (400.0, 1600.0) if all((D_, v) in d for v in nus)]:
        st = EQ.images(D, J, 4, 5, nus=np.array(nus), kinds=('EE2',), verbose=False)
        m = d[(D, nus[0])]['m']
        lap = []
        for nu in nus:
            r, rr = d[(D, nu)], d[(EQ.DREF, nu)]
            Ls = np.log(sum(r[k] for k in ('M0', 'M2', 'M4', 'M6')).ravel()[0]) + r['logC'] - r['logK']
            Lr = np.log(sum(rr[k] for k in ('M0', 'M2', 'M4', 'M6')).ravel()[0]) + rr['logC'] - rr['logK']
            lap.append(Ls - Lr)
        lap = np.array(lap)
        for Kt in (2, 3, 4):
            Ist, _ = EQ.value(st, 'EE2', Kt=Kt); Iref, _ = EQ.value(ref, 'EE2')
            bg = np.log(Ist/Iref)
            dd = m*m*((bg - bg[0]) - (lap - lap[0]))
            s = f'J={J} Delta={D} K={Kt}: m^2[(D_bg - D_bg(nu1)) - (D_K2 - D_K2(nu1))] at nu = ' + ' '.join(f'{v:.2f}:{x.real:+.3f}{x.imag:+.3f}i' for v, x in zip(nus, dd))
            s += '   | m^2 (D_bg - D_bg(nu1)) = ' + ' '.join(f'{x.real:+.1f}' for x in m*m*(bg - bg[0]))
            print(s, flush=True); lines.append(s)
    open(os.path.join(HERE, 'cache', 'bg_k2check.txt'), 'a').write('\n'.join(lines) + '\n')
