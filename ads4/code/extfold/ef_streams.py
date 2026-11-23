"""
ef_streams (EXT_FOLD.md §4.2): the J = 29 n3 small-nu nodes (nu = 10, 12.5, 15, 17.5), stream by stream.
Gamma (C2, C3H; c2, kappa subtracted, SMALLNU's baseline normalisation -- sn_core.small_nodes, unchanged) is evaluated with the n3
records of ONE distinct RQMC stream at a time (ef_merge's stream table: 4242..4249 = SMALLNU's job, 5001..5023 = the 16 ext jobs), the
rest of the row as in SMALLNU.  Then:
  * independent blocks: O = 4242-4249 (SMALLNU), A = 5001-5008 (ext job s5001), B = 5009-5016 (ext job s5009), C = 5017-5023 (the
    7 streams only ext jobs s5010-s5016 reach); block means, sigma = sd/sqrt(n), chi^2 of the four blocks against their weighted mean;
  * the 16 ext JOBS as submitted (overlapping streams): their per-job means, for the record (not independent: chi^2 meaningless);
  * pooled: equal-weight mean of the 31 streams, sigma = sd/sqrt(31), against ef_merge's 8-super-replica record (jackknife);
  * tails: max |z| of single streams (z = (Gamma_s - median)/(1.4826 MAD)), and the pooled value without the most extreme stream.
Usage: python3 ef_streams.py  -> cache/ef_streams.txt / .json
"""
import os, sys, json, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ef_core as C
import numpy as np
S, H = C.S, C.H
J, D = 29, 1029.0
NU3 = [10.0, 12.5, 15.0, 17.5]
ST = pickle.load(open(os.path.join(C.EFC, 'ef_streams_J29.pkl'), 'rb'))
KIND_OF = {'EE': 'C2', 'C3H': 'C3H'}
SEL = [None]          # list of stream ids to use for the n3 records (None: the set's records unchanged)
_p0 = C._pts


def _pts(kind, J_, D_):
    P = _p0(kind, J_, D_)
    if SEL[0] is None or J_ != J or D_ != D or S.BASE[0]:
        return P
    P = dict(P)
    for nu in NU3:
        for DD in (D, H.DREF):
            st = ST[kind][nu][DD]
            reps = np.array([st[s] for s in SEL[0]])
            if len(reps) < 8:
                reps = np.concatenate([reps]*int(np.ceil(8/len(reps))))[:8] if len(SEL[0]) == 1 else reps
            r = dict(P[(DD, nu)]); r['reps'] = reps; r['E'] = reps.mean(0); r['err'] = 0*np.abs(r['E']); r['R'] = len(reps)
            P[(DD, nu)] = r
    return P


C.HN.pts = _pts
C.S._pts = _pts


def gam(sel):
    SEL[0] = sel
    try:
        nodes, chk, da, m = S.small_nodes(J, D, numax=19.0, jkn=1)
    finally:
        SEL[0] = None
    return {round(r['nu'], 3): {k: r[k][0] for k in ('C2', 'C3H')} for r in nodes}, chk


if __name__ == '__main__':
    C.SET[0] = 'old'                     # rows exactly as SMALLNU; the n3 records are replaced stream by stream
    lines = []

    def pr(s=""):
        print(s, flush=True); lines.append(s)
    ids = sorted(set(ST['EE'][10.0][D]))
    per = {}
    for s in ids:
        per[s], chk = gam([s])
    pr(f"# ef_streams J=29 row 1029: Gamma - kappa (c2, baseline normalisation) per distinct n3 stream ({len(ids)} streams); transform chk {chk:.2e}")
    blocks = dict(O=list(range(4242, 4250)), A=list(range(5001, 5009)), B=list(range(5009, 5017)), C=list(range(5017, 5024)))
    res = dict(streams={str(s): per[s] for s in ids}, nu={})
    for nu in NU3:
        k = round(nu, 3)
        res['nu'][str(nu)] = {}
        for kd in ('C2', 'C3H'):
            v = np.array([per[s][k][kd] for s in ids])
            med = np.median(v); mad = 1.4826*np.median(np.abs(v - med)); z = (v - med)/mad
            bm, be = {}, {}
            for b, bl in blocks.items():
                vb = np.array([per[s][k][kd] for s in bl]); bm[b] = vb.mean(); be[b] = vb.std(ddof=1)/np.sqrt(len(vb))
            w = np.array([1/be[b]**2 for b in blocks]); mu = np.array([bm[b] for b in blocks])
            wm = (w*mu).sum()/w.sum(); chi2 = float((w*(mu - wm)**2).sum())
            from scipy.stats import chi2 as CH
            pv = float(CH.sf(chi2, len(blocks) - 1))
            # equal-variance version (one pooled sd) -- robust against block sd fluctuations
            sd = v.std(ddof=1); chi2e = float(sum(len(bl)*(bm[b] - v.mean())**2 for b, bl in blocks.items())/sd**2)
            pve = float(CH.sf(chi2e, len(blocks) - 1))
            imax = int(np.argmax(np.abs(z)))
            v_wo = np.delete(v, imax)
            o = dict(mean=float(v.mean()), sigma=float(sd/np.sqrt(len(v))), sd=float(sd), median=float(med), mad=float(mad),
                     zmax=float(z[imax]), smax=int(ids[imax]), mean_wo=float(v_wo.mean()), sigma_wo=float(v_wo.std(ddof=1)/np.sqrt(len(v_wo))),
                     blocks={b: (float(bm[b]), float(be[b])) for b in blocks}, chi2=chi2, p=pv, chi2_eqvar=chi2e, p_eqvar=pve)
            res['nu'][str(nu)][kd] = o
            pr(f"  nu={nu:5.1f} {kd:3s}: 31 streams mean {o['mean']:+7.1f} +- {o['sigma']:5.1f} (sd {sd:6.1f}; median {med:+7.1f}, MAD-sd {mad:6.1f}); "
               f"max |z| {abs(o['zmax']):4.1f} (stream {o['smax']}) -> without it {o['mean_wo']:+7.1f} +- {o['sigma_wo']:5.1f}")
            pr("        blocks " + "  ".join(f"{b}: {bm[b]:+7.1f} +- {be[b]:5.1f}" for b in blocks) +
               f" | chi2 = {chi2:.2f} (3 dof, p = {pv:.2f}); equal-variance chi2 = {chi2e:.2f} (p = {pve:.2f})")
    # the 16 ext jobs as submitted (overlapping streams)
    pr("#  ext jobs as submitted (Z2_SEED s uses streams s..s+7; jobs overlap, so these are NOT independent):")
    jobs = {}
    for s0 in range(5001, 5017):
        g, _ = gam(list(range(s0, s0 + 8)))
        jobs[s0] = g
    for nu in NU3:
        k = round(nu, 3)
        for kd in ('C2', 'C3H'):
            v = [jobs[s0][k][kd] for s0 in range(5001, 5017)]
            pr(f"  nu={nu:5.1f} {kd:3s}: job means " + " ".join(f"{x:+.0f}" for x in v) + f"  (spread sd {np.std(v, ddof=1):.0f})")
    res['jobs'] = {str(s0): jobs[s0] for s0 in jobs}
    # pooled with ALL 31 streams in one record (31 replicas; value = equal-weight ratio of means)
    g31, _ = gam(ids)
    res['pooled31'] = g31
    pr("#  pooled record, 31 streams equal weight (ratio of means): " + "; ".join(
        f"nu={nu}: C2 {g31[round(nu, 3)]['C2']:+.1f}, C3H {g31[round(nu, 3)]['C3H']:+.1f}" for nu in NU3))
    open(os.path.join(C.EFC, 'ef_streams.txt'), 'w').write("\n".join(lines) + "\n")
    json.dump(res, open(os.path.join(C.EFC, 'ef_streams.json'), 'w'), indent=1, default=float)
