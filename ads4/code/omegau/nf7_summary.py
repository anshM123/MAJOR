"""nf7: tables from the nf4 scans (Omega^{NF,SC}, NF02) -> cache/nf7_summary.txt.
Same definitions as u7_summary.py: kappa = (F - E)/U (abs, Cs = 0) per band; controlled margin; n(kappa_req, delta) worst over h and Df = 3, 2;
binding states; negative-kappa states; conditional bound sqrt(t4^2 + aT^2) <= 240 sqrt(n^8 16 L_AdS(Lam/n)/k)/Lam^4 with NF02's psi2 (factored)."""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from scipy.integrate import quad
from nf_common import NF

F = NF()
NG = np.round(np.concatenate([np.arange(1.0, 1.3, 0.01), np.arange(1.3, 3.0, 0.05), [3.0]]), 4)
KREQ = (0.05, 0.1, 0.3, 0.5, 1.0, 2.0)
out = []
def pr(s):
    print(s, flush=True); out.append(s)
def L_AdS(Lam):
    return quad(lambda t: float(F.psi2(t))/(t*t + 2.25/Lam**2), 0, 1, limit=400, points=[1.5/Lam])[0]
def bound(Dg, n):
    if not np.isfinite(n):
        return np.nan, np.nan
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); X = n**8*16*L_AdS(Lam/n)       # k(0) = 1 (normalised)
    return X, 240*np.sqrt(X)/Lam**4
def marg(r, kreq, delta, mode):
    core = r['Fh'] - r['E'] if mode == 'abs' else r['Fh'] + r['S'] - r['Es']
    return core - kreq*r['U'] - delta*r['Ns']
def n_of(L, kreq, delta, mode='abs'):
    for n0 in NG:
        if all(marg(r, kreq, delta, mode) > 0 for v in L for r in v if r['mfl'] >= n0 - 1e-9):
            return n0
    return np.inf
R = {}
for Dg in (20, 50, 100, 300, 1000):
    for Df in ('Df3', 'Df2'):
        R[(Dg, Df)] = pickle.load(open(os.path.join(HERE, 'cache', f'recs_nf4_{Dg}_{Df}.pkl'), 'rb'))
pr("== Omega^{NF,SC} (NF02): kappa (abs) per band [1,1.1) [1.1,1.2) [1.2,1.3) [1.3,1.5) | m>=1.5: min (F-E)/Ns, min kappa | min F_hyb/Ns | max NU/Ns | max U/Ns")
BANDS = ((1.0, 1.1), (1.1, 1.2), (1.2, 1.3), (1.3, 1.5))
for (Dg, Df), RR in R.items():
    for h, v in RR.items():
        ks = [min([(r['Fh'] - r['E'])/r['U'] for r in v if b[0] <= r['mfl'] < b[1] and r['U'] > 0], default=np.inf) for b in BANDS]
        hi = [r for r in v if r['mfl'] >= 1.5]
        lo = min([(r['Fh'] - r['E'])/r['Ns'] for r in hi if r['U'] == 0], default=np.inf)
        ka = min([(r['Fh'] - r['E'])/r['U'] for r in hi if r['U'] > 0], default=np.inf)
        pr(f"  {Df} Dg={Dg:4d} {h}: " + " ".join(f"{k:7.3g}" for k in ks) + f" | {lo:+.3f} {ka:6.3g} | {min(r['Fh']/r['Ns'] for r in v):+.4f} | "
           f"{max(r['NU']/r['Ns'] for r in v):.3g} | {max(r['U']/r['Ns'] for r in v):.3g}")
for mode in ('abs', 'coh'):
    pr(f"== [{mode}] n(kappa_req, delta), worst over h and Df = 3, 2; kappa_req = " + ", ".join(map(str, KREQ)) + "  ; bound sqrt(t4^2+aT^2) at that n")
    for Dg in (20, 50, 100, 300, 1000):
        L = [R[(Dg, Df)][h] for Df in ('Df3', 'Df2') for h in ('++', '+-')]
        for delta in (0.0, 0.03):
            ns = [n_of(L, kq, delta, mode) for kq in KREQ]
            pr(f"  Dg={Dg:4d} delta={delta:.2f}: n = " + " ".join(f"{x:5.2f}" for x in ns) + "   bound: " + " ".join(f"{bound(Dg, x)[1]:.2g}" for x in ns))
        for kq, delta in ((1.0, 0.0), (0.3, 0.03)):
            n1 = n_of(L, kq, delta, mode)
            cand = [r for v in L for r in v if r['mfl'] >= min(n1, 3.0) - 1e-9]
            b = min(cand, key=lambda r: marg(r, kq, delta, mode)/r['Ns'])
            h = '++' if any(b is r for Df in ('Df3', 'Df2') for r in R[(Dg, Df)]['++']) else '+-'
            pr(f"     binding at n(kappa_req={kq}, delta={delta}) = {n1} (or worst state above 3 Lam): m/Lam={b['mfl']:.3f}, J={b['J']}, {h}: F/Ns={b['Fh']/b['Ns']:.3f} "
               f"E/Ns={b['E']/b['Ns']:.3f} S/Ns={b['S']/b['Ns']:+.4f} Es/Ns={b['Es']/b['Ns']:.4f} U/Ns={b['U']/b['Ns']:.3g} spare={marg(b, kq, delta, mode)/b['Ns']:+.3f}")
pr("== states with F_hyb - E <= 0 (abs criterion fails even before the uncontrolled weight); per (Df, Dg, h): count, worst, and how many have m >= 3 Lam:")
for (Dg, Df), RR in R.items():
    for h, v in RR.items():
        bad = [r for r in v if r['Fh'] - r['E'] <= 0]
        if bad:
            w = min(bad, key=lambda r: (r['Fh'] - r['E'])/r['Ns'])
            pr(f"  {Df} Dg={Dg} {h}: {len(bad)} states ({sum(r['mfl'] >= 3 for r in bad)} at m >= 3 Lam), worst {(w['Fh'] - w['E'])/w['Ns']:+.4f} at m/Lam={w['mfl']:.3f}, J={w['J']} (bhat={(2*w['J']+1)/w['mfl']:.1f}, E/Ns={w['E']/w['Ns']:.3g}, |S|/Ns={abs(w['S'])/w['Ns']:.3g})")
pr("== bound vs n (NF02): " + "; ".join(f"Dg={Dg}: " + ", ".join(f"n={n}: {bound(Dg, n)[0]:.3g} ({bound(Dg, n)[1]:.2g})" for n in (1.0, 1.1, 1.2, 1.3, 1.4)) for Dg in (100, 1000)))
open(os.path.join(HERE, 'cache', 'nf7_summary.txt'), 'w').write("\n".join(out) + "\n")

# ---- where the failures sit: thresholds restricted to impact parameter bhat <= B (the large-bhat tail of NF02's action is a cancellation
#      ~ 1/bhat^3, far below the per-node error estimates of the truncated fixed-J / fixed-eta series there)
pr("== thresholds restricted to bhat = (2J+1)/(m/Lam) <= B  [mode: n(kappa_req=1, delta=0) / n(0.3, 0.03) / n(0.1, 0)]")
for mode in ('abs', 'coh'):
    for B in (5.0, 10.0, 20.0, 40.0):
        row = []
        for Dg in (20, 50, 100, 300, 1000):
            L = [[r for r in R[(Dg, Df)][h] if (2*r['J'] + 1)/r['mfl'] <= B] for Df in ('Df3', 'Df2') for h in ('++', '+-')]
            row.append("/".join(f"{n_of(L, kq, dl, mode):.2f}" for kq, dl in ((1.0, 0.0), (0.3, 0.03), (0.1, 0.0))))
        pr(f"  [{mode}] bhat <= {B:4.0f}: " + "  ".join(f"Dg={Dg}: {x}" for Dg, x in zip((20, 50, 100, 300, 1000), row)))
pr("== failing states (abs, F - E <= 0) by bhat bin, all Dg, Df:")
for lo, hi in ((0, 5), (5, 10), (10, 20), (20, 30), (30, 41)):
    tot = bad = 0
    for (Dg, Df), RR in R.items():
        for h, v in RR.items():
            for r in v:
                b = (2*r['J'] + 1)/r['mfl']
                if lo <= b < hi:
                    tot += 1; bad += r['Fh'] - r['E'] <= 0
    pr(f"  bhat in [{lo},{hi}): {bad}/{tot}")
open(os.path.join(HERE, 'cache', 'nf7_summary.txt'), 'w').write("\n".join(out) + "\n")
pr("== kappa (abs) per band restricted to bhat <= 10, worst over h and Df; and min (F-E)/Ns on controlled states with bhat <= 10")
for Dg in (20, 50, 100, 300, 1000):
    v = [r for Df in ('Df3', 'Df2') for h in ('++', '+-') for r in R[(Dg, Df)][h] if (2*r['J'] + 1)/r['mfl'] <= 10]
    ks = [min([(r['Fh'] - r['E'])/r['U'] for r in v if b[0] <= r['mfl'] < b[1] and r['U'] > 0], default=np.inf) for b in BANDS + ((1.5, 3.0), (3.0, 10.01))]
    lo = min([(r['Fh'] - r['E'])/r['Ns'] for r in v if r['U'] == 0], default=np.inf)
    pr(f"  Dg={Dg:4d}: " + " ".join(f"{k:7.3g}" for k in ks) + f"  (bands [1,1.1) [1.1,1.2) [1.2,1.3) [1.3,1.5) [1.5,3) [3,10]) | controlled min {lo:+.3f} | max U/Ns {max(r['U']/r['Ns'] for r in v):.3g}")
open(os.path.join(HERE, 'cache', 'nf7_summary.txt'), 'w').write("\n".join(out) + "\n")
