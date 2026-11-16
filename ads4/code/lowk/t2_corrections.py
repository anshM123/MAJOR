"""
Approach to 1 of the corner-completed ratios (k = 2, 3; k = 4, 6 for comparison): m^2 (R - 1), O(1/m), O(1/m^2), non-analytic terms.
Reads cache/data_t1_l{1,2}.pkl (t1_lowk_e2e.py).  FITTED coefficients (labelled as such in SPINNING_LOWK.md).

Weighted least squares (weights 1/err^2 with err = max(propagated error of R, 1e-9)) on the masses m >= 20, for R - 1:
  (B2) c2/m^2 + c4/m^4                        analytic in 1/m^2, two terms
  (B3) c2/m^2 + c4/m^4 + c6/m^6               analytic in 1/m^2, three terms
  (A)  c1/m + c2/m^2 + c4/m^4                 O(1/m) test (as e2e/t7)
  (N)  c2/m^2 + cN/m^{k+1/2} + c4/m^4         non-analytic term at the order of the E2E Assumption-T heuristic (u_0 ~ 1/m)
With four masses >= 20 (beta ~ 0.4, 1) the three-parameter fits have one degree of freedom; chi2 is printed.  With three masses
(beta ~ 2) they are exact.  Quoted parameter errors: statistical (from the propagated numerical errors) only; the truncation of the
ansatz is visible in the spread between (B2) and (B3) and in chi2.
Usage: python3 t2_corrections.py        Output: out_t2_corrections.txt
"""
import os, sys, pickle
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'e2e'))
import numpy as np
from corner import CACHE

def wls(ms, y, e, pows):
    A = np.stack([ms**(-p) for p in pows], 1)
    w = 1/np.maximum(e, 1e-9)
    Aw = A*w[:, None]; yw = y*w
    c, *_ = np.linalg.lstsq(Aw, yw, rcond=None)
    cov = np.linalg.pinv(Aw.T @ Aw)
    chi2 = float(np.sum((Aw @ c - yw)**2))
    return c, np.sqrt(np.diag(cov)), chi2

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
                    pr("#      m       R_pres  (err)            m^2(R_pres-1)        R_flip  (err)            m^2(R_flip-1)")
                    for v in rows:
                        pr(f"  {v['m']:7.2f}  {v['Rp']:.9f} ({v['dRp']:.0e})  {v['m']**2*(v['Rp']-1):10.3f}±{v['m']**2*v['dRp']:.1g}    "
                           f"{v['Rf']:.9f} ({v['dRf']:.0e})  {v['m']**2*(v['Rf']-1):10.3f}±{v['m']**2*v['dRf']:.1g}")
                    sel = [v for v in rows if v['m'] >= 19]
                    if len(sel) < 3:
                        continue
                    ms = np.array([v['m'] for v in sel])
                    for comp, key, ekey in (("pres", 'Rp', 'dRp'), ("flip", 'Rf', 'dRf')):
                        y = np.array([v[key] for v in sel]) - 1; e = np.array([v[ekey] for v in sel])
                        res = {}
                        for name, pows in (("B2", (2, 4)), ("B3", (2, 4, 6)), ("A", (1, 2, 4)), ("N", (2, k + 0.5, 4))):
                            if len(pows) > len(ms):
                                continue
                            res[name] = wls(ms, y, e, pows)
                        line = "   ".join(f"({n}) " + ", ".join(f"c[{p:g}]={c[i]:+.4g}±{s[i]:.1g}" for i, p in enumerate(pw)) + f" chi2={x2:.1g}"
                                           for n, pw in (("B2", (2, 4)), ("B3", (2, 4, 6)), ("A", (1, 2, 4)), ("N", (2, k + 0.5, 4)))
                                           for (c, s, x2) in [res[n]])
                        pr(f"#   {comp}: " + line)
                        # adopted c2: mid-range of the c2 values of the fits that are acceptable (chi2/dof < 4), error = half the range
                        # (statistical errors are not meaningful where chi2 >> 1; referee M2)
                        fits = {n: res[n] for n in res}
                        dof = {n: len(ms) - len(pw) for n, pw in (("B2", (2, 4)), ("B3", (2, 4, 6)), ("A", (1, 2, 4)), ("N", (2, k + 0.5, 4)))}
                        c2s = {"B2": fits["B2"][0][0], "B3": fits["B3"][0][0], "A": fits["A"][0][1], "N": fits["N"][0][0]}
                        okf = [n for n in c2s if (dof[n] == 0 or fits[n][2]/max(dof[n], 1) < 4)]
                        vals = [c2s[n] for n in okf] or list(c2s.values())
                        c2_mid, c2_err = 0.5*(max(vals) + min(vals)), 0.5*(max(vals) - min(vals))
                        summary.append((l, beta, k, nu, comp, c2s, {n: fits[n][2] for n in fits}, okf, c2_mid, c2_err,
                                        res['A'][0][0], res['A'][1][0], res['N'][0][1], res['N'][1][1]))
    pr("")
    pr("# ---- summary (fits on m >= 20): c2 of each ansatz with chi2 in brackets; adopted c2 = mid-range of the acceptable fits (chi2/dof < 4)")
    pr("#      +- half the range; c1 of (A) and the coefficient cN of m^{-k-1/2} in (N) with their statistical errors")
    pr("#  l  beta  k  nu  comp   c2 B2 [chi2]      c2 B3 [chi2]      c2 A [chi2]       c2 N [chi2]      adopted c2        c1(A)              cN(N)")
    for s_ in summary:
        (l, beta, k, nu, comp, c2s, chis, okf, mid, err, c1, c1e, cN, cNe) = s_
        cols = "  ".join(f"{c2s[n]:9.3f} [{chis[n]:.0e}]" for n in ("B2", "B3", "A", "N"))
        pr(f"  {l}  {beta:3.1f}  {k}  {nu:3.1f}  {comp:4s}  {cols}  {mid:9.3f}±{err:.2g}  {c1:+8.4f}±{c1e:.1e}  {cN:+10.3f}±{cNe:.1e}")
    open(os.path.join(HERE, 'out_t2_corrections.txt'), 'w').write("\n".join(out) + "\n")
