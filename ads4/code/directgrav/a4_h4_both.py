"""
a4: H4 / H1 metrics for BOTH helicities from the cleaned (t19) extraction r1 (K* even k, A odd k; D-safe Omega^{NF,SC} weights with K* as
the even kernel).  Targets (normalised to the leading Regge form 2P):
   ++ :  R_skel = ex_skel/(2P)                                              (H4 as used by exactskel e7)
   +- :  reading B:  R_skel + E_k (d^J_44 - P_J)/(2P)   (exactskel e3/e7 model);   reading A: R_skel   (for comparison)
Resolved node: cleaned K*/A model spread <= SPMAX and envelope width <= 10 SPMAX.  c* = median over resolved nodes of R_grav/R_target
(overall per-state normalisation; cannot change a sign).  Clean H4 shift dF = sum sigma^k u^{k-2} phi_k w eps_k 2P (R_grav/c* - R_target)
over resolved nodes, in units of N_s; also |dF| against the skeleton margin.  H1 check: R_grav itself (leading-formula ratio).
Usage: python3 a4_h4_both.py [SPMAX]  -> cache/a4_summary.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from a2_action import weights, KS, EPS, SIGMA, Dh_one, pn
from common import Ek_scalar
SPMAX = float(sys.argv[1]) if len(sys.argv) > 1 else 3e-3
out = []
def pr(s=""):
    print(s, flush=True); out.append(s)
e7 = pickle.load(open(os.path.join(HERE, '..', 'exactskel', 'cache', 'e7_100_recs.pkl'), 'rb'))
for tag in ('gw_100_1.040_4', 'gw_100_1.300_6'):
    a1 = pickle.load(open(os.path.join(HERE, 'cache', 'a1_%s.pkl' % tag), 'rb'))
    r1 = pickle.load(open(os.path.join(HERE, 'cache', 'r1_%s.pkl' % tag), 'rb'))
    mfl, J, m = a1['mfl'], a1['J'], a1['m']; u = 1/mfl**2
    pr(f"## {tag}: m/Lam={mfl} J={J} m={m:.2f}; resolved = spread <= {SPMAX:g}, envelope width <= {10*SPMAX:g}")
    for h, wave in (('++', ('same', 'same')), ('+-', ('opp', 'opp'))):
        for Df in (3.0, 2.0):
            ref = next(r for r in e7 if r['h'] == h and abs(r['mfl'] - mfl) < 1e-6 and r['J'] == J and float(r['Df']) == Df)
            for reading in (('B', 'A') if h == '+-' else ('-',)):
                rows = []
                for k in KS:
                    kn = 'A' if k % 2 else 'K*'
                    wk, env = weights('NF02', k, u)
                    for j, i in enumerate(a1['inodes']):
                        nu = a1['NUS'][j]
                        if (k, nu, wave, kn) not in r1:
                            continue
                        d = r1[(k, nu, wave, kn)]
                        res = d['spread'] <= SPMAX and (d['env'][1] - d['env'][0]) <= 10*SPMAX
                        x = nu/m; y = x*x; P = a1['P'][j]
                        Rs = a1[('skel', Df)][(k, nu)]['R_ex']
                        Rt = Rs + (Ek_scalar(k, y)*(Dh_one('+-', J, y) - Dh_one('++', J, y))/(2*P) if (h == '+-' and reading == 'B') else 0.0)
                        rows.append(dict(k=k, nu=nu, i=i, R=d['R'], Rt=Rt, res=res, w=wk[i], P=P, sp=d['spread']))
                rr = [r for r in rows if r['res']]
                if not rr:
                    pr(f"  {h} Df={Df:g} reading {reading}: no resolved nodes"); continue
                c = float(np.median([r['R']/r['Rt'] for r in rr]))
                dF = sum(SIGMA[h]**r['k']*r['w']*EPS[r['k']]*2*r['P']*(r['R']/c - r['Rt']) for r in rr)
                byk = {}
                for r in rows:
                    byk.setdefault(r['k'], []).append(r)
                cells = []
                for k, v in sorted(byk.items()):
                    vr = [r for r in v if r['res']]
                    numax = max([r['nu'] for r in vr], default=0)
                    dev = max([abs(r['R']/(c*r['Rt']) - 1) for r in vr], default=np.nan)
                    cells.append(f"k{k}: {len(vr)}/{len(v)} res (nu<={numax:.1f}) dev {dev:.1e}")
                pr(f"  {h} Df={Df:g} reading {reading}: c* = {c:.4f} (m^2(c*-1) = {m*m*(c-1):+.0f}); clean H4 shift dF = {dF/ref['Ns']:+.2e} Ns"
                   f" vs skeleton margin F/Ns = {ref['Fex']/ref['Ns']:+.4f}")
                pr("     " + "; ".join(cells))
            # H1 table (leading-formula ratio) at a few nodes, once per helicity
        pr(f"  H1 {h}: R_grav (cleaned, K* even / A odd) [spread] at nu ~ 0.03 / 1.3 / 3.5 / 6.6: " + "; ".join(
            f"k{k}: " + " ".join(f"{r1[(k, nu, wave, 'A' if k % 2 else 'K*')]['R']:.4f}[{r1[(k, nu, wave, 'A' if k % 2 else 'K*')]['spread']:.0e}]"
                                 for nu in a1['NUS'] if any(abs(nu - t) < 0.1 for t in (0.034, 1.314, 3.469, 6.599)) and (k, nu, wave, 'K*') in r1) for k in KS))
open(os.path.join(HERE, 'cache', 'a4_summary.txt'), 'w').write("\n".join(out) + "\n")
