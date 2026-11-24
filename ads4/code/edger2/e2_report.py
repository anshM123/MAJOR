"""
edger2/e2_report.py (EDGE_R2.md table): per state and entry, the (EDGE-R) cost_E = |cost_2| + |cost_3| split into
  determined part D_k (e2_gamma: bulk fit + measured K15 edge nodes) with its budget mc (3 sigma jackknife), repr (bulk degree),
  quad (3 |K15 - G7| of the measured edge sum),
  remainder, sign-free (S_E,t units):
    R_edge : the GK nodes beyond the last valid n1 node x_l (the n3-only nodes of the job list, never run), |Gamma| <= c A ((1 - x_l)/(1 - x))^q,
             A = max|Gamma| of the last three valid nodes of the key; main c = 3, q = 2 (HIGHSPIN's rule; e2_calib: the exact scalar needs
             c >= 1.2-2.1 pointwise at these J and x_l), variant q = the graviton's own measured growth exponent on [0.8, x_l];
    R_gap  : measured-range GK nodes without a usable value (reference-zero window, helicity off-diagonal > 2 %, rel. err > 0.2), covered
             by 3 x max|Gamma| of the 4 nearest valid edge nodes (HIGHSPIN's zero-node rule);
  upper = sum_k (|D_k| + mc + repr + quad + R_edge + R_gap);  lower = sum_k max(0, |D_k| - (same budget)).
Verdict: CLOSED if upper <= tol (grade [NE + M]: the envelope is a model), VIOLATED if lower > tol, else OPEN.
-> cache/report.json / .txt
"""
import e2_guard  # noqa: F401
import os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, 'cache')
STATES = [26, 30, 34, 33, 38, 42, 46, 50]


def growth_exponent(x, g, xl):
    s = (x >= 0.8) & (x <= xl + 1e-9)
    if s.sum() < 6:
        return np.nan
    xs, gs = x[s], np.abs(g[s])
    e = np.linspace(0.8, xl + 1e-9, 6)
    mids, mx = [], []
    for u, v in zip(e[:-1], e[1:]):
        t = (xs >= u) & (xs <= v)
        if t.any():
            mids.append(0.5*(u + v)); mx.append(gs[t].max())
    if len(mids) < 3:
        return np.nan
    return float(-np.polyfit(np.log(1 - np.array(mids)), np.log(mx), 1)[0])


def main():
    OUT = []
    pr = lambda s='': (print(s, flush=True), OUT.append(s))
    calib = json.load(open(os.path.join(CACHE, 'calib.json')))
    res = {}
    for J in STATES:
        d = json.load(open(os.path.join(CACHE, 'gamma_J%d.json' % J)))
        SEt, tol_tr = d['SEt'], d['tol_tr']
        rs = {}
        for h, ent in d['entries'].items():
            rk = {}
            for kd in ('C2', 'C3H'):
                e = ent[kd]
                node = d['nodes'][str((kd, h))]
                x = np.array(node['x']); g = np.array(node['G'])
                xl, A = e['xl'], e['A']
                q_meas = growth_exponent(x, g, xl)
                edge_x = x[x >= 0.8]
                Redge = {2.0: 0.0, 'meas': 0.0}; Rgap = 0.0; n_edge = n_gap = 0
                for (xm, wk, wg) in e['miss']:
                    if xm > xl + 1e-9:
                        n_edge += 1
                        Redge[2.0] += wk*3*A*max(1.0, ((1 - xl)/(1 - xm))**2)
                        qq = max(2.0, q_meas) if np.isfinite(q_meas) else 3.0
                        Redge['meas'] += wk*3*A*max(1.0, ((1 - xl)/(1 - xm))**qq)
                    else:
                        n_gap += 1
                        near = np.argsort(np.abs(edge_x - xm))[:4]
                        gz = 3*np.max(np.abs(g[x >= 0.8][near]))
                        Rgap += wk*gz
                rk[kd] = dict(D=e['det']/SEt, mc=e['mc']/SEt, repr=e['repr']/SEt, quad=e['quad']/SEt, Redge=Redge[2.0]/SEt,
                              Redge_q=Redge['meas']/SEt, Rgap=Rgap/SEt, q_meas=q_meas, n_edge=n_edge, n_gap=n_gap, xl=xl, A=A)
            bud = {kd: r['mc'] + r['repr'] + r['quad'] for kd, r in rk.items()}
            rem = {kd: r['Redge'] + r['Rgap'] for kd, r in rk.items()}
            det = sum(abs(r['D']) for r in rk.values())
            up = sum(abs(r['D']) + bud[kd] + rem[kd] for kd, r in rk.items())
            lo = sum(max(0.0, abs(r['D']) - bud[kd] - rem[kd]) for kd, r in rk.items())
            up_q = sum(abs(r['D']) + bud[kd] + r['Redge_q'] + r['Rgap'] for kd, r in rk.items())
            verdict = 'CLOSED' if up <= tol_tr else ('VIOLATED' if lo > tol_tr else 'OPEN')
            rs[h] = dict(k=rk, det=det, budget=sum(bud.values()), rem=sum(rem.values()), upper=up, upper_q=up_q, lower=lo,
                         verdict=verdict, comb=(rk['C2']['D'] + rk['C3H']['D']))
        res[J] = dict(tol=tol_tr, mfl=d['mfl'], x1=d['x1'], entries=rs)
        for h, r in rs.items():
            k2, k3 = r['k']['C2'], r['k']['C3H']
            pr(f"J = {J} (m/L {d['mfl']:.3f}, x1 {d['x1']:.4f}) entry {h}: tol {tol_tr:.3f} | D_2 {k2['D']:+.3f}, D_3 {k3['D']:+.3f} "
               f"(|D_2|+|D_3| = {r['det']:.3f}; D_2 + D_3 = {r['comb']:+.3f}) | budget of D: mc {k2['mc'] + k3['mc']:.3f} repr "
               f"{k2['repr'] + k3['repr']:.3f} quad {k2['quad'] + k3['quad']:.3f} | remainder: edge (x > {k2['xl']:.4f}, {k2['n_edge']} "
               f"nodes) {k2['Redge'] + k3['Redge']:.3f} [q = meas ({k2['q_meas']:.2f}/{k3['q_meas']:.2f}): "
               f"{k2['Redge_q'] + k3['Redge_q']:.3f}], gap ({k2['n_gap']} nodes) {k2['Rgap'] + k3['Rgap']:.3f} | upper {r['upper']:.3f} "
               f"(q meas: {r['upper_q']:.3f}), lower {r['lower']:.3f} -> {r['verdict']}   [S_E,t]")
    pr("\nscalar calibration (e2_calib): " + "; ".join(
        f"J = {J}: c_pt {v['c_pt']:.2f}, c_rem {v['c_rem']:.2f}, E_det {v['E_det']:+.2f} / E_rem {v['E_rem']:+.2f} / E_full {v['E_full']:+.3f} S_E2, "
        f"growth q {v['expo_meas']:.2f} (beyond x_l {v['expo_beyond']:.2f})" for J, v in calib.items()))
    json.dump({str(k): v for k, v in res.items()}, open(os.path.join(CACHE, 'report.json'), 'w'), default=float, indent=1)
    open(os.path.join(CACHE, 'report.txt'), 'w').write('\n'.join(OUT) + '\n')


if __name__ == '__main__':
    main()
