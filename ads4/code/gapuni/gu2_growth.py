"""
gu2_growth (GAP_UNIFORM.md §5): chart-A behaviour of the graviton beyond-hybrid profile Gamma(x; J, m) (c2 units, K2 rows) between
Delta ~ 1000, 1600 and 2000, read from gapnear/cache/gn_hmono_*.json (READ ONLY; another agent's cache), and the growth allowance
against the chart-A tolerance of gu1_tol (cache/gu1_tol.json).

Per (J, kind, helicity) and node x in [0.05, 0.9]:
  * the two-term law  Gamma(m) = g_inf + g_1/m  fitted on the ~1000 and ~2000 rows, checked on the 1600 row (residual);
  * the K2 - K1 step on each row (the last Laplace term kept; its m-independence in c2 is the Laplace structure of §4);
  * the model majorant  sup_{m >= m0} |Gamma| <= |g_inf| + |g_1|/m0  (if the two-term law held exactly)   [NE/M]
Per J: sup|Gamma| on the 1000 row, the majorant, and  F_allow = G_cert(J, twist edge) / sup|Gamma|_1000  (and with G_hyb) = the factor by
which |Gamma| may grow (uniformly in x) at ANY larger gap without losing positivity (Lemma T: the tolerance itself grows like Lam^2).
Usage: python3 gu2_growth.py -> cache/gu2_growth.txt / .json
"""
import os, sys, json, glob
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gu_guard  # noqa: F401
import numpy as np

GN = os.path.join(HERE, '..', 'gapnear', 'cache')
C = os.path.join(HERE, 'cache')
OUT = []


def pr(s=''):
    print(s, flush=True); OUT.append(s)


def mphys(D, J):
    al = D - J - 2.0
    return np.sqrt(al*(al + 2*J + 1))


def main():
    H = {}
    for fn in sorted(glob.glob(os.path.join(GN, 'gn_hmono_*.json'))):
        H.update(json.load(open(fn)))
    tol = json.load(open(os.path.join(C, 'gu1_tol.json')))
    res = {}
    pr("# per node: Gamma (c2, K=2) on rows lo/1600/hi; fit g_inf + g1/m on lo, hi; resid = fit(1600) - Gamma_1600; step = K2-K1")
    for J in sorted(H, key=int):
        Ji = int(J)
        keys = H[J]
        Ds = sorted({D for k in keys for D in keys[k]}, key=float)
        lo = [D for D in Ds if 985 <= float(D) < 1100]
        hi = [D for D in Ds if 1900 < float(D) < 2100]
        mid = [D for D in Ds if float(D) == 1600]
        if not hi:
            continue
        sup_lo = sup_hi = sup_maj = 0.0; max_res = 0.0; max_g1 = 0.0; max_step = [0.0, 0.0]
        for k, rows in keys.items():
            if 'off' in k:
                continue
            vh = {round(x, 4): (v, s) for x, v, s in rows[hi[0]] if 0.05 - 1e-9 <= x <= 0.9 + 1e-9}
            mh = mphys(float(hi[0]), Ji)
            for x, (v, s) in vh.items():
                sup_hi = max(sup_hi, abs(v)); max_step[1] = max(max_step[1], s)
            if not lo:
                continue
            ml = mphys(float(lo[0]), Ji)
            vl = {round(x, 4): (v, s) for x, v, s in rows[lo[0]] if 0.05 - 1e-9 <= x <= 0.9 + 1e-9}
            vm = {round(x, 4): v for x, v, s in rows.get(mid[0], [])} if mid else {}
            mm = mphys(1600.0, Ji)
            for x in sorted(set(vl) & set(vh)):
                a, sa = vl[x]; b, sb = vh[x]
                sup_lo = max(sup_lo, abs(a)); max_step[0] = max(max_step[0], sa)
                g1 = (a - b)/(1/ml - 1/mh); ginf = b - g1/mh
                max_g1 = max(max_g1, abs(g1))
                sup_maj = max(sup_maj, abs(ginf) + abs(g1)/ml)
                if x in vm:
                    max_res = max(max_res, abs(ginf + g1/mm - vm[x]))
        # tolerance at the 1000 twist edge (worst helicity) and min over the tabulated masses
        tj = [t for t in tol if t['J'] == Ji]
        if tj:
            edge = min(t['mu'] for t in tj)
            te = [t for t in tj if t['mu'] == edge]
            Gc = min(t['1000']['Gcert'] for t in te); Gh = min(t['1000']['Ghyb'] for t in te)
            Gc_all = min(t['1000']['Gcert'] for t in tj); Gh_all = min(t['1000']['Ghyb'] for t in tj)
        else:
            edge = Gc = Gh = Gc_all = Gh_all = float('nan')
        r = dict(J=Ji, rows=lo + mid + hi, sup_lo=sup_lo, sup_hi=sup_hi, sup_maj=sup_maj, max_g1=max_g1, max_res=max_res,
                 step_lo=max_step[0], step_hi=max_step[1], mu_edge=edge, Gcert_edge=Gc, Ghyb_edge=Gh, Gcert_min=Gc_all, Ghyb_min=Gh_all)
        res[Ji] = r
    pr()
    pr(f"{'J':>3} {'rows':>16} | {'sup|G|lo':>8} {'sup|G|hi':>8} {'max|g1|':>9} {'|g1|/m0':>7} {'res1600':>7} {'majorant':>8} | "
       f"{'step lo':>7} {'step hi':>7} | {'Gcert(edge)':>11} {'F_cert':>7} {'Ghyb(edge)':>10} {'F_hyb':>7} {'F_cert,min':>10} {'F_hyb,min':>9}")
    for Ji in sorted(res):
        r = res[Ji]
        base = max(r['sup_lo'], r['sup_hi'])
        f = lambda G: (G/base if (G == G and G > 0 and base > 0) else float('nan'))
        pr(f"{Ji:>3} {'/'.join(r['rows']):>16} | {r['sup_lo']:8.1f} {r['sup_hi']:8.1f} {r['max_g1']:9.3g} {r['max_g1']/1000:7.1f} "
           f"{r['max_res']:7.1f} {r['sup_maj']:8.1f} | {r['step_lo']:7.0f} {r['step_hi']:7.0f} | {r['Gcert_edge']:11.3g} {f(r['Gcert_edge']):7.1f} "
           f"{r['Ghyb_edge']:10.3g} {f(r['Ghyb_edge']):7.1f} {f(r['Gcert_min']):10.1f} {f(r['Ghyb_min']):9.1f}")
    json.dump(res, open(os.path.join(C, 'gu2_growth.json'), 'w'), indent=1)
    open(os.path.join(C, 'gu2_growth.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


def incremental():
    """incremental (signed) transfer criterion per 1000 state of gapnear's base table (GAP_UNIFORM §3.2, Corollary T'):
    lambda(Lam) >= lam_hyb,min - (Lam0/Lam)^2 [ |dlam0| + E0 + W0 sup_x|Gamma_Lam - Gamma_Lam0| ],
    with the two-term law sup|Gamma_Lam - Gamma_Lam0| <= |g1| (1/m0 - 1/m)  (|g1| = max over the J rows, gu2 main) [NE/M]:
    margin(t) = lam_min t^2 - |dlam0| - E0 - W0 |g1|/m0 (1 - 1/t),  t = Lam/Lam0 >= 1;  report min_t margin / lam_min."""
    sys.path.insert(0, os.path.join(HERE, '..', 'poserr'))
    sys.path.insert(0, os.path.join(HERE, '..', 'spinfunc'))
    import large_m_map as LM
    base = json.load(open(os.path.join(GN, 'gn_bind_base1000.json')))
    g = json.load(open(os.path.join(C, 'gu2_growth.json')))
    pr()
    pr("# incremental transfer: E0 = sign-free 1000 budget (E_tot) / lambda-level (E_tot_lam); W0 = max_h W_h + W_off (1000);"
       " |g1| from the 1000/2000 rows")
    pr(f"{'J':>3} {'m/Lam':>8} | {'lam_hyb':>8} {'dlam0':>9} {'E0':>8} {'E0lam':>8} | {'W0*Lam0^2':>10} {'|g1|/m0':>8} {'W0|g1|/m0':>10} | "
       f"{'min-margin/lam (sf)':>19} {'(lam-level)':>11}")
    L0 = LM.E.lam_of(1000.0)
    for r in base:
        J = r['J']
        if str(J) not in g:
            continue
        s = LM.state(1000.0, J, r['mfl'])
        hs = [k for k in ('++', '+-') if k in s]
        W0 = max(s[h]['W_lo'] + s[h]['W_hi'] for h in hs) + s['W_off_lo'] + s['W_off_hi']
        m0 = r['mfl']*L0
        inc = W0*g[str(J)]['max_g1']/m0
        lam_min = min(r['lam_hybrid'], min(s[h]['lam'] for h in hs))
        dl = abs(r['lam_K2'] - r['lam_hybrid'])
        out = []
        for E0 in (r['E_tot'], r['E_tot_lam']):
            ts = np.linspace(1, 20, 4000)
            mg = lam_min*ts**2 - dl - E0 - inc*(1 - 1/ts)
            out.append(mg.min()/lam_min)
        pr(f"{J:>3} {r['mfl']:8.5f} | {r['lam_hybrid']:8.5f} {dl:9.2e} {r['E_tot']:8.2e} {r['E_tot_lam']:8.2e} | {W0*L0**2:10.4g} "
           f"{g[str(J)]['max_g1']/m0:8.1f} {inc:10.2e} | {out[0]:+19.3f} {out[1]:+11.3f}")
    open(os.path.join(C, 'gu2_growth.txt'), 'w').write("\n".join(OUT) + "\n")


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'inc':
    OUT[:] = open(os.path.join(C, 'gu2_growth.txt')).read().splitlines()
    incremental()
