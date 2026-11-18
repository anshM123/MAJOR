"""
f3: near-gap D-safe action (NF02 / Omega^{NF,SC}) on all exact-skeleton states at Delta_g = 100, 300, 1000 (e1 + e7 caches), for the
graviton objects of spf_common, in the bulk-point-factor model of objid/o2 [M]:
    F^X = F_ex + sum_k sum_i (R_X(k, x_i) - 1) I_k(p_i),     I_k = dens + curv + rem  (spf_common.NearGap),
so 'A' (route A, R = 1) is exactly the e7 action.  Then the H4 / H2 tolerance targets for route A:
  * rem-tolerance delta_rem(x0) = min over states of F_ex / sum_{k, x_i > x0} |rem_k(p_i)|: the graviton remainder may differ from the scalar
    remainder by any relative amount |eps| <= delta_rem pointwise (any sign) on x > x0 without changing n = 1 (sufficient condition);
  * dens-tolerance delta_dens(x0) = min F_ex / sum_{k, x_i > x0} |dens_k + curv_k|: same for a relative error of the leading (B' + curvature)
    graviton density on x > x0;
  * shares: fraction of sum |rem| and of sum |dens| from x > x0 (Df = 3 and 2 pooled; per state, then the max/median over states).
Validation: the 64-node sum with R = 1 against e7's F_ex (400-node flat parts).
Usage: python3 f3_ads.py [Dg ...]  -> cache/f3_ads_<Dg>.txt, cache/f3_<Dg>_recs.pkl
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
from spf_common import R, OBJS, NearGap, load_states
from common import KS

X0S = (0.0, 0.1, 0.3, 0.5, 0.7, 0.9)

def run(Dg):
    NG = NearGap(Dg)
    recs = []
    for h, DF, r0, ref in load_states(Dg):
        x, P = NG.pieces(h, r0)
        tot = {k: P[k]['dens'] + P[k]['curv'] + P[k]['rem'] for k in KS}
        rec = dict(h=h, Df=DF, mfl=r0['mfl'], J=r0['J'], Ns=ref['Ns'], Fex=ref['Fex'], base64=sum(tot[k].sum() for k in KS))
        for obj in OBJS:
            rec['dF_' + obj] = sum(((R(obj, k, x) - 1)*tot[k]).sum() for k in KS)
        for x0 in X0S:
            m = x > x0
            rec[('rem', x0)] = sum(np.abs(P[k]['rem'][m]).sum() for k in KS)
            rec[('dens', x0)] = sum(np.abs(P[k]['dens'][m] + P[k]['curv'][m]).sum() for k in KS)
        rec['Srem'] = sum(P[k]['rem'].sum() for k in KS)
        recs.append(rec)
    pickle.dump(recs, open(os.path.join(HERE, 'cache', 'f3_%d_recs.pkl' % Dg), 'wb'))
    return recs

def report(Dg, recs):
    out = []
    def pr(s=""):
        print(s, flush=True); out.append(s)
    pr(f"# f3: Delta_g = {Dg}: near-gap D-safe action, graviton objects (bulk-point-factor model on the exact skeleton); {len(recs)} records")
    d = np.array([abs(r['base64'] - r['Fex'])/r['Ns'] for r in recs])
    pr(f"  validation (R = 1 node sum vs e7 F_ex): |diff|/N_s median {np.median(d):.1e}, max {d.max():.1e}")
    grid = sorted(set(r['mfl'] for r in recs))
    for obj in OBJS:
        F = lambda r: r['Fex'] + r['dF_' + obj]
        for delta in (0.0, 0.03):
            n = next((n0 for n0 in grid if all(F(r) - delta*r['Ns'] > 0 for r in recs if r['mfl'] >= n0 - 1e-9)), None)
            w = min(recs, key=lambda r: F(r)/r['Ns'])
            bad = [r for r in recs if F(r) - delta*r['Ns'] <= 0]
            pr(f"   {obj:6s} delta={delta}: n = {n}; min F/N_s = {F(w)/w['Ns']:+.4f} ({w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}); "
               f"{len(bad)} non-positive")
    # largest object shifts
    for obj in ('PsiG', 'PsiG2'):
        w = max(recs, key=lambda r: abs(r['dF_' + obj])/r['Ns'])
        pr(f"   largest |dF_{obj}|/N_s = {abs(w['dF_' + obj])/w['Ns']:.4f} at {w['h']} m/Lam={w['mfl']} J={w['J']} Df={w['Df']}")
    pr("  route A (R = 1): exact scalar remainder S/N_s: max |S|/N_s = %.4f" % max(abs(r['Srem'])/r['Ns'] for r in recs))
    pr("  H4/H2 tolerance targets for route A (min over all states; sufficient conditions for n = 1):")
    for x0 in X0S:
        drem = min(recs, key=lambda r: r['Fex']/max(r[('rem', x0)], 1e-300))
        dden = min(recs, key=lambda r: r['Fex']/max(r[('dens', x0)], 1e-300))
        sh_rem = np.array([r[('rem', x0)]/max(r[('rem', 0.0)], 1e-300) for r in recs])
        sh_den = np.array([r[('dens', x0)]/max(r[('dens', 0.0)], 1e-300) for r in recs])
        pr(f"   x > {x0:.1f}: delta_rem = {drem['Fex']/drem[('rem', x0)]:.3g} [{drem['h']} m/Lam={drem['mfl']} J={drem['J']}]; "
           f"delta_dens = {dden['Fex']/dden[('dens', x0)]:.3g} [{dden['h']} m/Lam={dden['mfl']} J={dden['J']}]; "
           f"share of sum|rem| median {np.median(sh_rem):.3f} max {sh_rem.max():.3f}; of sum|dens+curv| median {np.median(sh_den):.3f}")
    nearg = [r for r in recs if r['mfl'] <= 1.1]
    if nearg:
        pr("  same restricted to m <= 1.1 Lam:")
        for x0 in (0.0, 0.5, 0.9):
            drem = min(nearg, key=lambda r: r['Fex']/max(r[('rem', x0)], 1e-300))
            dden = min(nearg, key=lambda r: r['Fex']/max(r[('dens', x0)], 1e-300))
            pr(f"   x > {x0:.1f}: delta_rem = {drem['Fex']/drem[('rem', x0)]:.3g}; delta_dens = {dden['Fex']/dden[('dens', x0)]:.3g}")
    open(os.path.join(HERE, 'cache', 'f3_ads_%d.txt' % Dg), 'w').write("\n".join(out) + "\n")

if __name__ == "__main__":
    for Dg in (map(int, sys.argv[1:]) if len(sys.argv) > 1 else (100, 300, 1000)):
        report(Dg, run(Dg))
