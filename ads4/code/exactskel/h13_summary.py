"""
h13_summary: tables for NF02_LARGE_B.md from the caches of h13_hyb.py and h13_exact.py.
 (1) hybrid action (flat NF02 + leading AdS Regge curvature) on every integer J, per (Dg, helicity, m-band, bhat-band): min F_hyb/N_s;
 (2) asymptotic envelope check: F_hyb / [A1 M1/Lam^2 + A2 M2/(2 Lam^3)] on bhat >= 50, 100 (h13_asym);
 (3) exact scalar-skeleton remainder S on the sampled states: F/N_s, S/N_s, |S|/F_hyb, error estimates.
Output: cache/h13_summary.txt
"""
import os, sys, pickle, glob
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from h13_common import CACHE, Lam_of
import h13_asym as HA

out = []
def pr(s=''):
    print(s, flush=True); out.append(s)

# M1, M2 tables
pk = os.path.join(CACHE, 'h13_M12.pkl')
if os.path.exists(pk):
    BT, M1T = pickle.load(open(pk, 'rb'))
else:
    BT = np.unique(np.concatenate([np.logspace(-3, 0, 60), np.linspace(1.0, 6.0, 51)]))
    M1T = np.array([HA.Mn(b, 1) for b in BT])
    pickle.dump((BT, M1T), open(pk, 'wb'))
from scipy.interpolate import CubicSpline
_cs = CubicSpline(np.log(BT), M1T)
M1 = lambda b: _cs(np.log(b))
M2 = lambda b: -(1/(2*np.sqrt(2)))*(np.cosh(b) - 1)**-1.5

MB = ((1.5, 2.0), (2.0, 3.0), (3.0, 5.0), (5.0, 10.01))
BB = ((20, 40), (40, 70), (70, 120), (120, 200.01), (200, 400), (400, 1000.1))
pr('# (1) hybrid D-safe action F_hyb = flat NF02 + leading AdS Regge curvature, every integer J (h13_hyb.py)')
for fn in sorted(glob.glob(os.path.join(CACHE, 'h13_hyb_*_*_*.pkl'))):
    tag = os.path.basename(fn)[8:-4]
    Dg = int(tag.split('_')[0])
    db = pickle.load(open(fn, 'rb'))
    recs = [r for v in db.values() for r in v]
    L = Lam_of(Dg)
    pr(f"## Dg = {Dg}, bhat range {tag.split('_')[1]}-{tag.split('_')[2]}: {len(recs)} (state, helicity) records, "
       f"{len(set(r['mfl'] for r in recs))} masses in [{min(r['mfl'] for r in recs):.3f}, {max(r['mfl'] for r in recs):.2f}] Lam; "
       f"non-positive: {sum(r['Fh'] <= 0 for r in recs)}")
    for h in ('++', '+-'):
        R = [r for r in recs if r['h'] == h]
        w = min(R, key=lambda r: r['Fh']/r['Ns'])
        pr(f"   {h}: min F_hyb/N_s = {w['Fh']/w['Ns']:+.4f} at m/Lam = {w['mfl']:.3f}, J = {w['J']}, bhat = {w['bh']:.1f}")
        for mb in MB:
            row = []
            for bb in BB:
                v = [r for r in R if mb[0] <= r['mfl'] < mb[1] and bb[0] <= r['bh'] < bb[1]]
                row.append(f"{min(r['Fh']/r['Ns'] for r in v):.3f}" if v else '  -  ')
            pr(f"      m/Lam [{mb[0]:4.1f},{mb[1]:5.2f}): " + ' '.join(f"b[{bb[0]},{int(bb[1])}) {c}" for bb, c in zip(BB, row)))
    # asymptotic envelope
    for b0 in (50, 100, 200):
        v = [r for r in recs if r['bh'] >= b0]
        if not v:
            continue
        rat = []
        for r in v:
            u = 1/r['mfl']**2
            A1, A2 = HA.Acoef(r['h'], u, 1), HA.Acoef(r['h'], u, 2)
            asym = A1*M1(r['beta'])/L**2 + A2*M2(r['beta'])/(2*L**3)
            rat.append(r['Fh']/asym)
        rat = np.array(rat)
        pr(f"   (2) bhat >= {b0}: F_hyb / [A1 M1/L^2 + A2 M2/(2L^3)] in [{rat.min():.5f}, {rat.max():.5f}]  ({len(v)} records)")
    pr('')

pr('# (3) exact scalar-skeleton remainder (h13_exact.py): F = F_hyb + S')
for fn in sorted(glob.glob(os.path.join(CACHE, 'h13_ex_*.pkl'))):
    if 'OLD' in fn:
        continue
    db = pickle.load(open(fn, 'rb'))
    if not db:
        continue
    Dg = list(db.values())[0]['Dg']
    pr(f"## Dg = {Dg}: {len(db)} (state, skeleton) records")
    pr("   Df  m/Lam     J   bhat   beta   NP |   ++: F/Ns     S/Ns    |S|/F_hyb  |   +-: F/Ns     S/Ns    |S|/F_hyb  | n-sum err/Ns  odd-k/Ns  sum|.|/Ns")
    allr = []
    for key in sorted(db, key=lambda k: (k[2], k[1], k[0]), reverse=False):
        r = db[key]
        cells = []
        errs, odds, abss = [], [], []
        for h in ('++', '+-'):
            v = r['hel'].get(h)
            if v is None:
                cells.append(' ' * 36); continue
            cells.append(f"{v['Fex']/v['Ns']:+.5f}  {v['S']/v['Ns']:+.2e}  {abs(v['S']/v['Fh']):.1e}")
            errs.append(abs(v['S'] - v.get('Salt', v['S']))/v['Ns']); odds.append(abs(v['S'] - v['S3'])/v['Ns']); abss.append(v['Sabs']/v['Ns'])
            allr.append((v['Fex']/v['Ns'], h, r, v, max(abs(w['S'] - w.get('Salt', w['S']))/w['Ns'] for w in r['hel'].values())))
        flag = '  UNRELIABLE (n-sum)' if max(errs) > 1e-3 else ''
        pr(f"   {r['Df']:3g} {r['mfl']:6.3f} {r['J']:5d} {r['bh']:6.1f} {r['beta']:6.3f} {r.get('NP', 64):4d} | {cells[0]} | {cells[1]} | "
           f"{max(errs):.0e}  {max(odds):.0e}  {max(abss):.1e}{flag}")
    bad = [a for a in allr if a[4] > 1e-3]
    allr = [a for a in allr if a[4] <= 1e-3]
    pr(f"   reliable (n-sum err/N_s <= 1e-3 in both channels): {len(allr)} helicity records; unreliable: {len(bad)}")
    w = min(allr, key=lambda a: a[0])
    pr(f"   min F/N_s = {w[0]:+.5f} ({w[1]}, Df = {w[2]['Df']:g}, m/Lam = {w[2]['mfl']:.3f}, J = {w[2]['J']}, bhat = {w[2]['bh']:.1f});"
       f"  max |S|/F_hyb = {max(abs(a[3]['S']/a[3]['Fh']) for a in allr):.2e};  max |S|/N_s = {max(abs(a[3]['S']/a[3]['Ns']) for a in allr):.2e};"
       f"  non-positive: {sum(a[3]['Fex'] <= 0 for a in allr)}")
    # empirical scaling of the remainder: C_eff = |S|/F_hyb * m^2 (m in AdS units)
    for bb in ((20, 40), (40, 80), (80, 201)):
        v = [a for a in allr if bb[0] <= a[2]['bh'] < bb[1]]
        if v:
            w = max(v, key=lambda a: abs(a[3]['S']/a[3]['Fh'])*a[2]['m']**2)
            pr(f"   bhat in [{bb[0]},{bb[1]}): max |S|/F_hyb * m^2 = {abs(w[3]['S']/w[3]['Fh'])*w[2]['m']**2:.1f}  "
               f"({w[1]}, Df={w[2]['Df']:g}, m/Lam={w[2]['mfl']:.2f}, J={w[2]['J']});  max |S|/F_hyb * Lam^2 = "
               f"{max(abs(a[3]['S']/a[3]['Fh']) for a in v)*Lam_of(Dg)**2:.1f}")
    pr('')
open(os.path.join(CACHE, 'h13_summary.txt'), 'w').write('\n'.join(out) + '\n')

# (4) condensed markdown table for NF02_LARGE_B.md: per state, worst over both skeletons
md = []
for fn in sorted(glob.glob(os.path.join(CACHE, 'h13_ex_*.pkl'))):
    if 'OLD' in fn:
        continue
    db = pickle.load(open(fn, 'rb'))
    if not db:
        continue
    Dg = list(db.values())[0]['Dg']
    md.append(f"\n**Δ_g = {Dg}**\n")
    md.append("| `m/Λ` | `J` | `b̂` | `β` | `++` `F/N_s` | `+−` `F/N_s` | max `|S|/N_s` | max `|S|/F_hyb` | `n`-sum err/`N_s` | status |")
    md.append("|---|---|---|---|---|---|---|---|---|---|")
    states = sorted(set((k[0], k[1]) for k in db), key=lambda a: (a[1], a[0]))
    for mflt, J in states:
        recs = [db[k] for k in db if (k[0], k[1]) == (mflt, J)]
        cell = {}
        smax = sfmax = emax = 0.0
        good = []
        for r in recs:
            e = max(abs(v['S'] - v.get('Salt', v['S']))/v['Ns'] for v in r['hel'].values())
            emax = max(emax, e)
            if e > 1e-3:
                continue
            good.append(r['Df'])
            for h, v in r['hel'].items():
                cell[h] = min(cell.get(h, 9e9), v['Fex']/v['Ns'])
                smax = max(smax, abs(v['S']/v['Ns'])); sfmax = max(sfmax, abs(v['S']/v['Fh']))
        r = recs[0]
        f = lambda h: (f"{cell[h]:+.4f}" if h in cell else "—")
        st = ('ok (Δφ=' + ','.join(f"{d:g}" for d in sorted(good, reverse=True)) + ')') if good else '**unreliable** (n-sum)'
        md.append(f"| {r['mfl']:.2f} | {J} | {r['bh']:.1f} | {r['beta']:.3f} | {f('++')} | {f('+-')} | "
                  + (f"{smax:.1e} | {sfmax:.1e}" if good else "— | —") + f" | {emax:.0e} | {st} |")
open(os.path.join(CACHE, 'h13_exact_table.md'), 'w').write('\n'.join(md) + '\n')
