"""
Approach to 1 of the corner-completed k = 2, 3 (and 4, 6) ratios: m^2 (R - 1) sequences, O(1/m) and O(1/m^2) coefficients, and a test for
non-analytic terms.  Reads cache/data_t1_l{1,2}.pkl (t1_lowk_e2e.py).

Ansaetze for R - 1 (three largest m, solved exactly; the fourth largest m is predicted and the misfit printed):
  (A)  c1/m + c2/m^2 + c4/m^4                (O(1/m) test, as in e2e/t7)
  (B)  c2/m^2 + c4/m^4 + c6/m^6              (analytic in 1/m^2: the expectation)
  (N)  c2/m^2 + c_{k+1/2}/m^{k+1/2} + c4/m^4 (a non-analytic term at the order suggested by the heuristic of Assumption T,
                                               u_0 ~ 1/m:  int_0^{1/m} u^{k-1/2} du ~ m^{-k-1/2})
  (N2) c2/m^2 + c_{2k+1}/m^{2k+1} + c4/m^4    (non-analytic term for a breakdown scale u_0 ~ 1/m^2)
The coefficient errors are propagated from the per-point errors of R (worst case over the sign pattern).
Usage: python3 t2_corrections.py        Output: out_t2_corrections.txt
"""
import os, sys, pickle, itertools
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from corner import CACHE

def solve(ms, R, pows):
    A = np.stack([ms**(-p) for p in pows], 1)
    return np.linalg.solve(A, R)

def solve_err(ms, R, dR, pows):
    c0 = solve(ms, R, pows)
    worst = np.zeros_like(c0)
    for sg in itertools.product((-1, 1), repeat=len(ms)):
        worst = np.maximum(worst, np.abs(solve(ms, R + np.array(sg)*dR, pows) - c0))
    return c0, worst

if __name__ == "__main__":
    out = []
    def pr(*a):
        s = " ".join(str(x) for x in a); print(s, flush=True); out.append(s)
    summary = []
    for l in (1, 2):
        fn = os.path.join(CACHE, f'data_t1_l{l}.pkl')
        if not os.path.exists(fn):
            continue
        data = pickle.load(open(fn, 'rb'))
        for beta in (0.4, 1.0, 2.0):
            for k in (2, 3, 4, 6):
                for nu in (0.0, 1.0, 3.0):
                    rows = sorted([v for key, v in data.items() if key[0] == beta and key[1] == k and key[2] == nu], key=lambda v: v['m'])
                    if not rows:
                        continue
                    pr(f"# l={l} beta~{beta} k={k} nu={nu}")
                    pr("#      m       R_pres  (err)            m^2(R_pres-1)      R_flip  (err)            m^2(R_flip-1)")
                    for v in rows:
                        pr(f"  {v['m']:7.2f}  {v['Rp']:.9f} ({v['dRp']:.0e})  {v['m']**2*(v['Rp']-1):12.4f}    {v['Rf']:.9f} ({v['dRf']:.0e})  {v['m']**2*(v['Rf']-1):12.4f}")
                    if len(rows) < 4:
                        continue
                    ms = np.array([v['m'] for v in rows])
                    for comp, key, ekey in (("pres", 'Rp', 'dRp'), ("flip", 'Rf', 'dRf')):
                        R = np.array([v[key] for v in rows]) - 1; dR = np.array([v[ekey] for v in rows])
                        line = []
                        res = {}
                        for name, pows in (("A", (1, 2, 4)), ("B", (2, 4, 6)), ("N", (2, k + 0.5, 4)), ("N2", (2, 2*k + 1, 4))):
                            c, e = solve_err(ms[-3:], R[-3:], dR[-3:], pows)
                            pred4 = sum(c[i]*ms[-4]**(-p) for i, p in enumerate(pows))
                            mis = (pred4 - R[-4])/max(abs(R[-4]), 1e-300)
                            res[name] = (c, e, mis)
                            line.append(f"({name}) " + ", ".join(f"c[{p:g}]={c[i]:+.4g}±{e[i]:.1g}" for i, p in enumerate(pows)) + f"; misfit@m={ms[-4]:.0f}: {mis:+.1e}")
                        pr(f"#   {comp}: " + "   ".join(line))
                        summary.append((l, beta, k, nu, comp, res['A'][0][0], res['A'][1][0], res['B'][0][0], res['B'][1][0],
                                        res['N'][0][1], res['N'][1][1], res['B'][2], res['N'][2], res['N2'][0][1], res['N2'][1][1]))
    pr("")
    pr("# ---- summary: c1 of (A); c2 of (B); coefficient of m^{-k-1/2} in (N) and of m^{-2k-1} in (N2), with propagated errors;")
    pr("#      misfit of (B) and (N) at the fourth largest m (relative to R-1)")
    pr("#  l  beta  k  nu  comp     c1(A)             c2(B)                 c_{k+1/2}(N)            misfit B    misfit N     c_{2k+1}(N2)")
    for s in summary:
        pr(f"  {s[0]}  {s[1]:3.1f}  {s[2]}  {s[3]:3.1f}  {s[4]:4s}  {s[5]:+8.4f}±{s[6]:.1e}  {s[7]:10.3f}±{s[8]:.1e}  {s[9]:+11.3f}±{s[10]:.1e}  {s[11]:+.1e}  {s[12]:+.1e}   {s[13]:+.3e}±{s[14]:.1e}")
    open(os.path.join(HERE, 'out_t2_corrections.txt'), 'w').write("\n".join(out) + "\n")
