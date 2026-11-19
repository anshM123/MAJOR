"""
p8_neargap (PARITY_ODD_ROUTEA.md sect. 11): odd-sector near-gap test at Delta_g = 100 with the p8 grids, and its RQMC error.
  (1) p6_odd_neargap.main with P6_DIRECT = 0 (sect. 6 rows), 'k2x' (+ extra rows, Laplace only), '1' (+ direct values on x >= 0.2);
  (2) RQMC error propagation: NMC draws of independent Gaussian perturbations of every direct Gamma node (sigma = replica error), p6 rerun
      on each (GRID_OVERRIDE) -> spread of min F/N_s per J and of t*;
  (3) the a-independent crystal sign (p7) on the direct Delta = 100 row and on the extra rows (row's own m, x = p sqrt(u), u = (Lam/m)^2).
Usage: L1_DEG=6 python3 p8_neargap.py [NMC]  -> cache/p8_neargap.txt
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
lines = []


def pr(s=""):
    print(s, flush=True); lines.append(s)


def fmt(res):
    r = res[100]
    return (f"n = 1: {r['ok']}; min F/N_s {r['min']:+.4f} at (m/Lam, J, Df) = {r['wmin']}; t* = ({r['tstar'][0]:+.2f}, {r['tstar'][1]:+.2f}); "
            "per J: " + ", ".join(f"J={J}: {r['perJ'][J]:+.4f} (m/Lam {r['perJ_arg'][J]})" for J in (5, 7, 9)))


pr("# p8_neargap: odd-J (+-) near-gap test at Delta_g = 100, C2 + C3H, CHMRS de-shift, K2 chain; variants of the Delta_g = 100 rows")
for tag, lab in (('0', 'sect. 6 rows (Laplace K2; rows 100, 200, ...)'), ('k2x', '+ rows m = 1.05/1.1/1.2/1.3 Lam (Laplace K2 only)'),
                 ('1', '+ direct RQMC values on x >= 0.2 of rows 100 and 1.05-1.3 Lam, K2 below')):
    P6.DIRECT = tag; P6.GRID_OVERRIDE.clear()
    res = P6.main([100], write=(tag != '0'), quiet=True)
    pr(f"  [{lab}] {fmt(res)}")

# (2) noise propagation on the direct grid
G0 = pickle.load(open(os.path.join(HERE, 'cache', 'p4_grid_direct%s.pkl' % ('' if A.DESHIFT else '_noshift')), 'rb'))
rng = np.random.default_rng(11)
mins = {J: [] for J in (5, 7, 9)}; tl, th, gmin = [], [], []
P6.DIRECT = '1'
for it in range(NMC):
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
pr(f"  RQMC propagation ({NMC} draws, independent Gaussian per direct node, sigma = replica error): min F/N_s {np.mean(gmin):+.4f} +- {np.std(gmin):.4f}; "
   + ", ".join(f"J={J}: {np.mean(v):+.4f} +- {np.std(v):.4f}" for J, v in mins.items())
   + f"; t* lower {np.mean(tl):+.2f} +- {np.std(tl):.2f}, upper {np.mean(th):+.2f} +- {np.std(th):.2f}")

# (3) a-independent crystal sign per row (p7 formula), direct grid vs Laplace on the same rows
from spf_common import NearGap
NG = NearGap(100); F = NG.F; pn, pw, Lm = NG.pn, NG.pw, NG.Lam
GK = pickle.load(open(os.path.join(HERE, 'cache', 'p4_grid_k2x%s.pkl' % ('' if A.DESHIFT else '_noshift')), 'rb'))
pr("  a-independent crystal action F_cr/(a N) at each Delta_g = 100 row (row's own m; p7 formula): direct grid | Laplace K2 same rows")
for J in (5, 7, 9):
    for D in sorted(G0[J]):
        if D > 140:
            continue
        out = []
        for g in (G0[J][D], GK[J].get(D)):
            if g is None:
                out.append(np.nan); continue
            m = g['m']; u = (Lm/m)**2; x = pn*np.sqrt(u)
            p2 = F.psi2(pn)*pw; p3 = F.psi3(pn)*pw
            Dh = A.Dc('+-', J, x)
            G2 = CubicSpline(g['xs'], g['C2'])(np.minimum(x, g['xs'][-1])); G3 = CubicSpline(g['xs'], g['C3H'])(np.minimum(x, g['xs'][-1]))
            f2 = (p2*A.env('+-', x, 2)*(Dh + G2/m**2)).sum(); f3 = (p3*u*A.env('+-', x, 3)*(-Dh + G3/m**2)).sum()
            n2 = abs((p2*A.env('+-', x, 2)*Dh).sum()); n3 = abs((p3*u*A.env('+-', x, 3)*Dh).sum())
            out.append((f2 + f3)/(n2 + n3))
        pr(f"    J={J} Delta={D:.4f} (m/Lam = {G0[J][D]['m']/Lm:.3f}): {out[0]:+.4f} | {out[1]:+.4f}")
open(os.path.join(HERE, 'cache', 'p8_neargap.txt'), 'w').write("\n".join(lines) + "\n")
