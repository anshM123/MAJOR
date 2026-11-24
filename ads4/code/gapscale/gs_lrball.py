"""
gapscale/gs_lrball.py (GAP_SCALING.md §3.1): Lemma LR-inf UNIFORMLY in the gap.  lrtight/lrt3_rcert.py is re-run with Lam replaced by
an Arb BALL [Lam_a, Lam_b]; every bound it evaluates is then an enclosure valid for every Lam in the ball (ball arithmetic; the
min-selection in lrt_common.Weight.R compares upper endpoints, so it stays an upper bound).  Adjacent balls cover [Lam_lo, Lam_hi].
Per ball we record sup_{beta >= B} r(beta) for B in BETAS (certified cells of lrt3 to beta = 64 plus its tail majorant), and the
smallest certified beta0 with r < 1 / rtilde < 1.  The ball width adapts: it grows x2 (ratio up to RMAX) while sup r(beta >= 0.37) < 0.5.
Inputs copied read-only (Lam-independent): lrtight/cache/lrt1_constants.json, lrt2_m1.json.  Writes only gapscale/cache/lrball/.
Usage: python3 gs_lrball.py [Lam_hi=1e6] -> cache/lrball.txt, cache/lrball.json
"""
import os, sys, json, shutil, time
from fractions import Fraction as Fr
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
LRT = os.path.join(HERE, '..', 'lrtight')
sys.path.insert(0, LRT)
out = os.path.join(HERE, 'cache', 'lrball')
os.makedirs(out, exist_ok=True)
for f in ('lrt1_constants.json', 'lrt2_m1.json'):
    shutil.copy(os.path.join(LRT, 'cache', f), os.path.join(out, f))
import lrt_common as LC
LC.CACHE = out
import lrt3_rcert as L3
L3.CACHE = out
from flint import arb
L3.pr = lambda s='': None                      # silence the per-run report
BETAS = ('0.245', '0.3', '0.37', '0.5', '1')
LAM_LO = Fr(998*999)                           # Lam^2 at Delta_g = 1000
LAM_HI = float(sys.argv[1]) if len(sys.argv) > 1 else 1e6
RMAX = 2.0


def run_ball(L2a, L2b):
    lam = LC.A(L2a).sqrt().union(LC.A(L2b).sqrt())
    L3.LAM = lam
    L3.main()
    d = json.load(open(os.path.join(out, 'lrt3_rcert.json')))
    rows = d['rows']
    res = {}
    for B in BETAS:
        res['sup_r_' + B] = max([R['r'] for R in rows if Fr(R['ba']) >= Fr(B) - Fr(1, 10**9)] + [d['tail']])
        res['sup_rt_' + B] = max([R['rt'] for R in rows if Fr(R['ba']) >= Fr(B) - Fr(1, 10**9)] + [d['tail']])
    res['beta0_r'] = d['beta0_r_1.0']; res['beta0_rt'] = d['beta0_rt_1.0']
    return res


def main():
    t0 = time.time()
    lines = []

    def pr(s):
        print(s, flush=True); lines.append(s)
    pr("# gs_lrball: certified Lemma LR-inf remainder ratio r(beta) on Lam-BALLS covering [998.4999, %.3g]" % LAM_HI)
    pr("# Lam_a  Lam_b | sup r (beta >= 0.245 / 0.3 / 0.37 / 0.5 / 1) | beta0 (r<1) | beta0 (rtilde<1, F_LR>0)")
    L2a = LAM_LO; q = Fr(1001, 1000)
    allres = []
    worst = {B: 0.0 for B in BETAS}
    while True:
        La = float(L2a)**0.5
        Lb = min(La*float(q), LAM_HI)
        L2b = Fr(Lb*Lb).limit_denominator(10**6)
        if L2b <= L2a:
            break
        r = run_ball(L2a, L2b)
        r.update(Lam_a=La, Lam_b=float(L2b)**0.5)
        allres.append(r)
        for B in BETAS:
            worst[B] = max(worst[B], r['sup_r_' + B])
        pr(f"{La:12.2f} {float(L2b)**0.5:12.2f} | " + " ".join(f"{r['sup_r_' + B]:.3e}" for B in BETAS) +
           f" | {r['beta0_r']} | {r['beta0_rt']}")
        if float(L2b)**0.5 >= LAM_HI*(1 - 1e-12):
            break
        L2a = L2b
        if r['sup_r_0.37'] < 0.5 and q < RMAX:
            q = min(Fr(RMAX), q*q)
    pr("# worst sup r over all balls: " + ", ".join(f"beta>={B}: {worst[B]:.4e}" for B in BETAS))
    pr("# => Lemma LR-inf holds with r < 1 at beta >= 0.37 for EVERY Lam in [998.4999, %.3g] iff the beta>=0.37 entry is < 1. [%.0fs]"
       % (LAM_HI, time.time() - t0))
    json.dump(allres, open(os.path.join(HERE, 'cache', 'lrball.json'), 'w'), indent=0)
    open(os.path.join(HERE, 'cache', 'lrball.txt'), 'w').write("\n".join(lines) + "\n")


if __name__ == '__main__':
    main()
