"""
a3: the a2 action shift restricted to the RESOLVED (k, nu) nodes, ++ only (the +- moments are not resolved: DIRECT_GRAVITON_NEARGAP.md sect. 3).
Resolved set (from c1 and the smoothness of R_grav in nu, a1 tables): k = 2: nu <= 1.92; k >= 3: nu <= NUCUT (default 6.6).
Kernels: NF02 with K* (even) / A (odd) [proxy], NF02 with S_ff (as defined; k = 2-4 unresolved), M02 with D1^2 (k = 2-4 unresolved).
Reports dF/Ns, the overall factor c* (median R_grav/R_skel on the resolved set), the shape shift dF_shape = dF - (c*-1) F_part, and the
pointwise residual |delta rho| after removing c*, against d_crit = F_ex / sum_all |w env| (uniform pointwise remainder that would cancel the
skeleton margin).
Usage: python3 a3_restricted.py cache/a1_*.pkl [NUCUT]
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.special import j0
from a2_action import weights, pn, KS, EPS, SIGMA, ek, Dh_one, rel
args = [a for a in sys.argv[1:] if a.endswith('.pkl')]
NUCUT = float(next((a for a in sys.argv[1:] if not a.endswith('.pkl')), 6.6))
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
for fn in args:
    a1 = pickle.load(open(fn, 'rb')); Dg, mfl, J, m = a1['Dg'], a1['mfl'], a1['J'], a1['m']; u = 1/mfl**2
    pr(f"## {os.path.basename(fn)}: m/Lam={mfl} J={J}; ++; resolved nodes k=2: nu<=1.92, k>=3: nu<={NUCUT}")
    for Df in (3.0, 2.0):
        e1 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e1_%d_Df%g.pkl' % (Dg, Df)), 'rb'))
        r0 = e1[(round(mfl, 6), J)]
        e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_%d_recs.pkl' % Dg), 'rb'))
        e3 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e3_%d_recs.pkl' % Dg), 'rb'))
        for fun, keven in (('NF02', 'KS'), ('NF02', 'FF'), ('M02', 'E')):
            if fun == 'NF02':
                ref = next(r for r in e7 if r['h'] == '++' and abs(r['mfl'] - mfl) < 1e-6 and r['J'] == J and float(r['Df']) == Df)
            else:
                ref = next(r for r in e3[Df][('M02', '++')] if abs(r['mfl'] - mfl) < 1e-6 and r['J'] == J)
            dF = Fpart = 0.0; rows = []; ratios = []; cov = 0.0; full = 0.0
            for k in KS:
                wk, env = weights(fun, k, u); full += np.sum(np.abs(wk*env))
                kern = 'O' if k % 2 else keven
                for j, i in enumerate(a1['inodes']):
                    nu = a1['NUS'][j]
                    if nu > (1.92 if k == 2 else NUCUT) + 1e-9:
                        continue
                    x = nu/m; y = x*x; P = a1['P'][j]
                    Rg = rel(*a1['grav'][(k, nu)][(('same', 'same'), kern)]).real
                    Rs = a1[('skel', Df)][(k, nu)]['R_ex']
                    M = ek(k, y)*Dh_one('++', J, y) + 2*EPS[k]*(P - j0((2*J + 1)*x)) + env[i]*r0['rho'][k][i]
                    G = EPS[k]*2*P*Rg
                    dF += wk[i]*(G - M); Fpart += wk[i]*M; cov += abs(wk[i]*env[i])
                    rows.append((k, nu, Rg, Rs, (G - M)/env[i], M/env[i])); ratios.append(Rg/Rs)
            c = float(np.median(ratios))
            # clean H4 metric (no odd-k envelope convention): shift eps_k 2P (R_grav/c* - R_skel) of the exact moments themselves
            dH4 = 0.0; devk = {}
            for k in KS:
                wk, env = weights(fun, k, u); kern = 'O' if k % 2 else keven
                for j, i in enumerate(a1['inodes']):
                    nu = a1['NUS'][j]
                    if nu > (1.92 if k == 2 else NUCUT) + 1e-9:
                        continue
                    Rg = rel(*a1['grav'][(k, nu)][(('same', 'same'), kern)]).real; Rs = a1[('skel', Df)][(k, nu)]['R_ex']
                    dH4 += wk[i]*EPS[k]*2*a1['P'][j]*(Rg/c - Rs)
                    devk[k] = max(devk.get(k, 0.0), abs(Rg/(c*Rs) - 1))
            res = [d - (c - 1)*Mm for (k, nu, Rg, Rs, d, Mm) in rows]
            Ns = ref['Ns']
            pr(f"  Df={Df:g} {fun}[even:{keven}]: F/Ns(skel) {ref['Fex']/Ns:+.4f}; resolved coverage {cov/full:.4f}; dF = {dF/Ns:+.2e} Ns;"
               f" c* = {c:.4f}; dF_shape = {(dF - (c - 1)*Fpart)/Ns:+.2e} Ns; max|delta rho| raw {max(abs(r[4]) for r in rows):.1e},"
               f" after c* {max(abs(x) for x in res):.1e}; d_crit = {ref['Fex']/full:.1e}")
            byk = {}
            for (k, nu, Rg, Rs, d, Mm), rr in zip(rows, res):
                byk.setdefault(k, []).append((nu, Rg, Rs, rr))
 
            pr(f"     clean H4 shift (exact moments, c* removed): {dH4/Ns:+.2e} Ns;  per k max |R_grav/(c* R_skel) - 1|: " + " ".join(f"k{k} {v:.1e}" for k, v in sorted(devk.items())))
            pr("     per k: max |R_grav - R_skel| , max residual after c*: " + "; ".join(
                f"k{k}: {max(abs(a[1] - a[2]) for a in v):.1e}, {max(abs(a[3]) for a in v):.1e}" for k, v in sorted(byk.items())))
open(os.path.join(HERE, 'cache', 'a3_summary.txt'), 'w').write("\n".join(out) + "\n")
