"""u7: compact tables from the u4 scans (cache/recs_u4_<Dg>_Df<Df>.pkl) -> cache/u7_summary.txt.
Per (functional, Dg, skeleton, h): kappa (abs criterion, Cs = 0) per near-gap band, min controlled margin (F - E)/Ns and min kappa for
m >= 1.5 Lam, and the supported threshold n(kappa_req) = smallest grid m/Lam above which every state is positive (controlled) or has
kappa >= kappa_req, for kappa_req = 1 (a 100 % relative error on uncontrolled nodes tolerated) and 0.3; abs / coh and Cs = 0 / 40."""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from u4_scan import crit

NG = np.round(np.concatenate([np.arange(1.0, 1.3, 0.01), np.arange(1.3, 3.0, 0.05), [3.0]]), 4)
out = []
def pr(s):
    print(s, flush=True); out.append(s)
def nsup(v, req, Cs, mode):
    for n0 in NG:
        if all((x > 0) if c else (x >= req) for c, x in (crit(r, Cs, mode) for r in v if r['mfl'] >= n0 - 1e-9)):
            return n0
    return np.inf
BANDS = ((1.0, 1.1), (1.1, 1.2), (1.2, 1.3), (1.3, 1.5))
for name in ('M02', 'M01', 'Fg3'):
    pr(f"== {name}: kappa (abs, Cs=0) in bands [1,1.1) [1.1,1.2) [1.2,1.3) [1.3,1.5) | m>=1.5: min (F-E)/Ns, min kappa | n(kappa>=1) abs/coh, Cs=0/40 | n(kappa>=0.3) abs")
    for Df in ('Df3', 'Df2'):
        for Dg in (20, 50, 100, 300, 1000):
            f = os.path.join(HERE, 'cache', f'recs_u4_{Dg}_{Df}.pkl')
            if not os.path.exists(f):
                continue
            R = pickle.load(open(f, 'rb'))
            for h in ('++', '+-'):
                v = R[(name, h)]
                ks = []
                for b in BANDS:
                    c = [crit(r, 0.0, 'abs') for r in v if b[0] <= r['mfl'] < b[1]]
                    ks.append(min([x for ok, x in c if not ok], default=np.inf))
                hi = [crit(r, 0.0, 'abs') for r in v if r['mfl'] >= 1.5]
                lo = min([x for ok, x in hi if ok], default=np.inf); ka = min([x for ok, x in hi if not ok], default=np.inf)
                ns = [nsup(v, 1.0, Cs, mode) for mode in ('abs', 'coh') for Cs in (0.0, 40.0)]
                pr(f"  {Df} Dg={Dg:4d} {h}: " + " ".join(f"{k:7.3g}" for k in ks) + f" | {lo:+.3f} {ka:6.3g} | n1 = " + "/".join(f"{x:.2f}" for x in ns)
                   + f" | n0.3 = {nsup(v, 0.3, 0.0, 'abs'):.2f}")
open(os.path.join(HERE, 'cache', 'u7_summary.txt'), 'w').write("\n".join(out) + "\n")

# ---------------------------------------------------------------- after REFEREE_OMEGA_U.md M1: threshold as a function of the tolerance
# Criterion (abs, Cs = 0): every scanned state with m >= n Lam has  Fh - E - kappa_req U - delta Ns > 0,
#   kappa_req = tolerated relative AdS error on the uncontrolled nodes (per k, per node),
#   delta     = allowance for what E and U omit (graviton transfer I4/I5, odd-k prefactor I3: ~1-3 % Ns near 1.4 Lam).
# n is the worst case over helicities and the two skeletons (Df = 3, 2).  Also the conditional bound sqrt(t4^2+aT^2) at that n (M02).
import json
from fractions import Fraction as Fr
from scipy.integrate import quad
def n_of(recs_list, kreq, delta):
    for n0 in NG:
        if all(r['Fh'] - r['E'] - kreq*r['U'] - delta*r['Ns'] > 0 for v in recs_list for r in v if r['mfl'] >= n0 - 1e-9):
            return n0
    return np.inf
Jd = json.load(open(os.path.join(HERE, '..', 'robust', 'functionals_certified.json')))
def bound(name, Dg, n):
    if not np.isfinite(n) or name != 'M02':
        return np.nan
    c = [float(Fr(x)) for x in Jd['M02']['psi2_exact']]; k0 = 10884217/396000
    Lam = np.sqrt((Dg - 2.0)*(Dg - 1.0)); L = Lam/n
    X = n**8*16*quad(lambda t: sum(ci*t**i for i, ci in enumerate(c))/(t*t + 2.25/L**2), 0, 1, limit=400, points=[1.5/L])[0]/k0
    return 240*np.sqrt(X)/Lam**4
KREQ = (0.05, 0.1, 0.3, 0.5, 1.0, 2.0)
for name in ('M02', 'M01', 'Fg3'):
    pr(f"== {name}: n(kappa_req, delta) [worst over h and Df = 3, 2]; kappa_req = " + ", ".join(map(str, KREQ)) + " ; (M02: sqrt(t4^2+aT^2) bound at that n)")
    for Dg in (20, 50, 100, 300, 1000):
        L = []
        for Df in ('Df3', 'Df2'):
            f = os.path.join(HERE, 'cache', f'recs_u4_{Dg}_{Df}.pkl')
            if os.path.exists(f):
                R = pickle.load(open(f, 'rb')); L += [R[(name, '++')], R[(name, '+-')]]
        for delta in (0.0, 0.03):
            ns = [n_of(L, kq, delta) for kq in KREQ]
            pr(f"  Dg={Dg:4d} delta={delta:.2f}: n = " + " ".join(f"{x:5.2f}" for x in ns)
               + ("   bound: " + " ".join(f"{bound(name, Dg, x):.2g}" for x in ns) if name == 'M02' else ""))
    # binding states at the kappa_req = 1 threshold (delta = 0): spare margin (F - E - U)/Ns
    for Dg in (100,):
        L = []
        for Df in ('Df3', 'Df2'):
            R = pickle.load(open(os.path.join(HERE, 'cache', f'recs_u4_{Dg}_{Df}.pkl'), 'rb')); L += [R[(name, '++')], R[(name, '+-')]]
        n1 = n_of(L, 1.0, 0.0)
        b = min((r for v in L for r in v if r['mfl'] >= n1 - 1e-9 and r['U'] > 0), key=lambda r: (r['Fh'] - r['E'] - r['U'])/r['Ns'])
        pr(f"  Dg=100 binding state at n(kappa_req=1) = {n1}: m/Lam={b['mfl']}, J={b['J']}: F/Ns={b['Fh']/b['Ns']:.3f} E/Ns={b['E']/b['Ns']:.3f} "
           f"U/Ns={b['U']/b['Ns']:.3f} spare (F-E-U)/Ns={(b['Fh']-b['E']-b['U'])/b['Ns']:.3f}")
# negative kappa (some uncontrolled state with F - E < 0) anywhere below 1.5 Lam
pr("== states with F_hyb - E < 0 (negative kappa or negative controlled margin), any m/Lam:")
for name in ('M02', 'M01', 'Fg3'):
    for Df in ('Df3', 'Df2'):
        for Dg in (20, 50, 100, 300, 1000):
            R = pickle.load(open(os.path.join(HERE, 'cache', f'recs_u4_{Dg}_{Df}.pkl'), 'rb'))
            for h in ('++', '+-'):
                bad = [r for r in R[(name, h)] if r['Fh'] - r['E'] <= 0]
                if bad:
                    w = min(bad, key=lambda r: (r['Fh'] - r['E'])/r['Ns'])
                    pr(f"  {name} {Df} Dg={Dg} {h}: {len(bad)} states, worst (F-E)/Ns = {(w['Fh']-w['E'])/w['Ns']:+.4f} at m/Lam={w['mfl']}, J={w['J']} (U/Ns={w['U']/w['Ns']:.3g})")
open(os.path.join(HERE, 'cache', 'u7_summary.txt'), 'w').write("\n".join(out) + "\n")
