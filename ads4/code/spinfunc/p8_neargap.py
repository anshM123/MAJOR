"""
p8_neargap (PARITY_ODD_ROUTEA.md sect. 11): odd-sector near-gap test at Delta_g = 100 with the p8 grids, and its RQMC error.
  (1) p6_odd_neargap.main with P6_DIRECT = 0 (sect. 6 rows), 'k2x' (+ extra rows, Laplace only), '1' (+ direct values on x >= 0.2);
  (2) RQMC error propagation: NMC draws of independent Gaussian perturbations of every direct Gamma node (sigma = replica error), p6 rerun
      on each (GRID_OVERRIDE) -> spread of min F/N_s per J and of t*;
  (3) the a-independent crystal sign (p7) on the direct Delta = 100 row and on the extra rows (row's own m, x = p sqrt(u), u = (Lam/m)^2).
Usage: L1_DEG=6 python3 p8_neargap.py [NMC [PARTS]]  (PARTS: digits of the parts to run, default 123) -> cache/p8_neargap_<PARTS>.txt
"""
import os, sys, pickle, copy
sys.dont_write_bytecode = True
os.environ.setdefault('L1_DEG', '6')
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.interpolate import CubicSpline
import p6_odd_neargap as P6
import p4_odd_analysis as A

NMC = int(sys.argv[1]) if len(sys.argv) > 1 else 40
PARTS = sys.argv[2] if len(sys.argv) > 2 else '123'      # '13': skip the RQMC propagation; '2': only the propagation
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def fmt(res):
    r = res[100]
    return (f"n = 1: {r['ok']}; min F/N_s {r['min']:+.4f} at (m/Lam, J, Df) = {r['wmin']}; t* = ({r['tstar'][0]:+.2f}, {r['tstar'][1]:+.2f}); "
            "per J: " + ", ".join(f"J={J}: {r['perJ'][J]:+.4f} (m/Lam {r['perJ_arg'][J]})" for J in (5, 7, 9)))


pr("# p8_neargap: odd-J (+-) near-gap test at Delta_g = 100, C2 + C3H, CHMRS de-shift, K2 chain; variants of the Delta_g = 100 rows")
for tag, lab in ([] if '1' not in PARTS else (('0', 'sect. 6 rows (Laplace K2; rows 100, 200, ...)'),
                 ('direct_100', 'MAIN: sect. 6 rows, direct RQMC values on x >= 0.2 of the Delta = 100 row, K2 below'),
                 ('direct_100_k3low_afix', 'as MAIN, x < 0.2 of the Delta = 100 row from K3 (J = 5, 7, 9), a_{+-} kept from K2'),
                 ('k2x', 'sect. 6 rows + rows m = 1.05/1.1/1.2/1.3 Lam (Laplace K2; generic a-fit, unstable there)'),
                 ('1', '+ direct RQMC values on x >= 0.2 of rows 100 and 1.05-1.3 Lam, K2 below (same a caveat)'))):
    P6.DIRECT = tag; P6.GRID_OVERRIDE.clear()
    try:
        res = P6.main([100], write=(tag != '0'), quiet=True)
        pr(f"  [{lab}] {fmt(res)}")
    except FileNotFoundError as ex:
        pr(f"  [{lab}] grid not built ({ex})")

# (2) noise propagation on the direct grid
G0 = pickle.load(open(os.path.join(HERE, 'cache', 'p4_grid_direct_100%s.pkl' % ('' if A.DESHIFT else '_noshift')), 'rb'))
rng = np.random.default_rng(11)
mins = {J: [] for J in (5, 7, 9)}; tl, th, gmin = [], [], []
P6.DIRECT = 'direct_100'
for it in range(NMC if '2' in PARTS else 0):
    P6.GRID_OVERRIDE.clear()
    for J in (5, 7, 9):
        g = copy.deepcopy(G0[J])
        for D, rec in g.items():
            for kd in ('C2', 'C3H'):
                e = np.nan_to_num(rec[(kd, 'err')], nan=0.0)
                rec[kd] = rec[kd] + e*rng.standard_normal(len(e))
        P6.GRID_OVERRIDE[J] = g
    r = P6.main([100], write=False, quiet=True)[100]
    for J in mins:
        mins[J].append(r['perJ'][J])
    tl.append(r['tstar'][0]); th.append(r['tstar'][1]); gmin.append(r['min'])
P6.GRID_OVERRIDE.clear()
if NMC and '2' in PARTS:
  pr(f"  RQMC propagation ({NMC} draws, independent Gaussian per direct node, sigma = replica error): min F/N_s {np.mean(gmin):+.4f} +- {np.std(gmin):.4f}; "
   + ", ".join(f"J={J}: {np.mean(v):+.4f} +- {np.std(v):.4f}" for J, v in mins.items())
   + f"; t* lower {np.mean(tl):+.2f} +- {np.std(tl):.2f}, upper {np.mean(th):+.2f} +- {np.std(th):.2f}")

# (3) a-independent crystal action per row: S = (sum psi2 G2 + sum psi3 G3)/(|sum psi2 G2| + |sum psi3 G3|), G = the reference-normalised
#     crystal value (state/ref ratio x env x D_ref; no a, no de-shift: P4_DESHIFT=0 grids), at the row's own m (x = p sqrt(u), u = (Lam/m)^2);
#     and the normalisation-free direct-vs-K2 gap on x >= 0.2, m^2 max|dG/env| / max|G_K2/env| (~ c2 units x LO shape).
import p8_analysis as P8
from spf_common import NearGap
NG = NearGap(100); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
A.DESHIFT = False
pr("  a-free crystal action per row (literal crystal integral, no de-shift): S = signed/abs; variants: K2 everywhere | direct x>=0.2 + K2 below"
   " | direct x>=0.2 + K3 below (where computed); and the a-free direct-K2 gap on x >= 0.2 (C2 / C3H)")
for J in ((5, 7, 9) if '3' in PARTS else ()):
    gd = P8.grid(J, lowk=2); gk = P8.grid(J, direct=False, lowk=2); g3 = P8.grid(J, lowk=3)
    for D in sorted(gd):
        if D > 140:
            continue
        Ss = []
        for g in (gk[D], gd[D], g3[D]):
            m = g['m']; u = (Lm/m)**2; x = pn*np.sqrt(u); a = g['a']
            Dh = A.Dc('+-', J, x)
            G2 = CubicSpline(g['xs'], g['C2'])(np.minimum(x, g['xs'][-1])); G3 = CubicSpline(g['xs'], g['C3H'])(np.minimum(x, g['xs'][-1]))
            f2 = (a*p2*A.env('+-', x, 2)*(Dh + G2/m**2)).sum(); f3 = (a*p3*u*A.env('+-', x, 3)*(-Dh + G3/m**2)).sum()
            Ss.append((f2 + f3)/(abs(f2) + abs(f3)))
        gaps = []
        for kd, k in (('C2', 2), ('C3H', 3)):
            g = gd[D]; m = g['m']; dm = g[(kd, 'direct')]
            sg = -1.0 if k == 3 else 1.0
            Gd = g['a']*(sg*A.Dc('+-', J, g['xs']) + g[kd]/m**2); Gk = g['a']*(sg*A.Dc('+-', J, g['xs']) + g[(kd, 'K2same_a')]/m**2)
            gaps.append(m*m*np.abs(Gd - Gk)[dm].max()/np.abs(Gk).max() if dm.any() else np.nan)
        pr(f"    J={J} Delta={D:9.4f} (m/Lam = {gd[D]['m']/Lm:.3f}): S = {Ss[0]:+.3f} | {Ss[1]:+.3f} | {Ss[2]:+.3f} (K3 nodes {g3[D]['nk3']}); "
           f"direct-K2 gap {gaps[0]:.0f} / {gaps[1]:.0f}; direct nodes {gd[D]['ndirect']}")
A.DESHIFT = True
open(os.path.join(HERE, 'cache', 'p8_neargap_%s.txt' % PARTS), 'w').write("\n".join(lines) + "\n")
