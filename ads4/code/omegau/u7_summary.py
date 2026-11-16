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
