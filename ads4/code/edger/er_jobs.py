"""
edger/er_jobs.py (EDGE_R.md §6): the external job list ads4/ext/jobs_batch3.txt (format of ext/jobs.txt: TAG METHOD J DELTA SEED
nu1,nu2,...) for the graviton (EDGE-R) quantity at Theorem W's sampled states with J <= 51 (AG_PROOF chi.json rows).

Per state (J, m = mfl Lam; row Delta from alpha(alpha + 2J + 1) = m^2, alpha = Delta - J - 2):
  * support of chi: p >= max(p_w, p_e - 3 sigma) (chi(p_e - 3 sigma) = 1e-5), i.e. x >= x_lo = that / mfl;
  * bulk x in [x_lo, 0.8): the 0.025 grid (HIGHSPIN's direct-row layout; degree-6 shape fit), method n1;
  * edge x in [0.8, x1 = 1/mfl]: two GK15 windows [0.8, x_mid], [x_mid, x1] (x_mid = (0.8 + x1)/2; K15 values used directly, G7 for
    the quadrature error) -- n1 where nu/Delta < 0.93, n3 where nu/Delta >= 0.93 (n1 fails above ~0.94: EDGE_R §2);
  * validation: the two highest nodes with nu/Delta in [0.89, 0.93) also with n3 (n1/n3 overlap), and every n3 node with a second,
    independent seed (seeds spaced by 8: README_EXT note) for the n3 jackknife;
  * reference-zero guard: a node within 0.004 of x0 (the Delta_ref = 2e4 reference's first partial-wave zero, HIGHSPIN §2) is moved
    to x0 +- 0.004 (grid) or the GK window is shifted by 0.006 (edge);
  * states with x1 <= 0.93 Delta/m need no n3 at all (every m >= 1.13 Lam state).
Priority tags: P1 = the three thinnest gap states (26, 1.026), (30, 1.030), (34, 1.033); P2 = the other gap-mass states; P3 = m >= 1.13 Lam.
Usage: python3 er_jobs.py   -> writes ads4/ext/jobs_batch3.txt (the only file outside edger/cache this note writes; created, not edited)
                               and cache/jobs_summary.txt
"""
import er_core as R
import os, sys
import numpy as np
import hs_core as HC
import builtins

EXT = os.path.join(R.CODE, '..', 'ext', 'jobs_batch3.txt')


def delta_of(J, m):
    a = (-(2*J + 1) + np.sqrt((2*J + 1)**2 + 4*m*m))/2
    return a + J + 2


def ref_zero(J, D):
    h = '+-' if J % 2 else '++'
    from scipy.optimize import brentq
    f = lambda t: HC.E.Dc_fixed(h, J, t)
    ts = np.linspace(1e-5, 0.3, 6000); v = f(ts)
    i = int(np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0][0])
    return brentq(f, ts[i], ts[i + 1])*HC.HN.mphys(20000.0, J)/HC.HN.mphys(D, J)


def gk(a, b):
    return list(HC.gk15_nodes(a, b)[0])


def state_nodes(J, mfl, pe, sg, pw):
    L_ = R.lam(); m = mfl*L_; D = delta_of(J, m); x1 = 1/mfl
    x0 = ref_zero(J, D)
    xlo = max(pw, pe - 3*sg)/mfl
    g = np.arange(np.floor(xlo/0.025)*0.025, 0.8 - 1e-9, 0.025)
    g = g[g >= max(xlo - 0.025, 0.025)]
    g = np.array([x0 + np.sign(x - x0)*0.004 if abs(x - x0) < 0.004 else x for x in g])
    a = 0.8
    xm = 0.5*(a + x1)
    for _ in range(4):
        e = gk(a, xm) + gk(xm, x1)
        if min(abs(np.array(e) - x0)) >= 0.004:
            break
        a += 0.006; xm = 0.5*(a + x1)
    e = np.array(sorted(set(np.round(e, 9))))
    if a > 0.8:
        g = np.concatenate([g, np.arange(0.8, a - 1e-9, 0.025)])
    xs = np.concatenate([g, e])
    nud = xs*m/D
    n1 = xs[nud < 0.93]; n3 = xs[nud >= 0.93]
    ov = xs[(nud >= 0.89) & (nud < 0.93)][-2:]
    return dict(J=J, mfl=mfl, m=m, D=D, x1=x1, x0=x0, xlo=xlo, n1=n1, n3=n3, ov=ov, edge_a=a)


def main():
    rows = [r for r in R.chi_rows() if r['J'] <= 51]
    lines = ['# ext/jobs_batch3.txt -- EDGE_R.md §6: real-axis graviton rows for (EDGE-R) at Theorem W\'s sampled states with J <= 51',
             '# format as jobs.txt: TAG METHOD J DELTA SEED nu1,nu2,...  (one job per nu).  n1 seeds are irrelevant (fixed Sobol seeds);',
             '# n3 seeds: every n3 node twice (seeds 6000 + 8k and 6000 + 8k + 4096), independent streams (README_EXT note: spacing >= 8).',
             '# Priority: P1 (thinnest tolerances) first, then P2, then P3.  Measured cost (one core, EDGE_R.md §6): n1 (M = 11) 2.5 min/node',
             '# at J = 26; n3 (M = 10) 8-10 min/node.  P3 (m >= 1.13 Lam, x1 <= 0.887) is n1-only: it needs no edge method.',
             '# CAUTION (EDGE_R.md §2.4): n3 is validated for the scalar and the EE (k = 2) kind only; its C3H (k = 3) STATE values are wrong',
             '# (x200-1000 too small, 30-600 % errors) at nu/Delta >= 0.93 for both parities.  The *_n3 jobs give the k = 2 edge; the k = 3',
             '# edge needs a fixed C3H evaluator first (run P3 and the *_n1 jobs of P1/P2 now; the *_n3 jobs only after that fix, or for k = 2).']
    summ = []
    seed = 6000
    order = sorted(rows, key=lambda r: (0 if (r['J'], round(r['mfl'], 2)) in ((26, 1.03), (30, 1.03), (34, 1.03)) else
                                        (1 if r['mfl'] < 1.06 else 2), r['J'], r['mfl']))
    nn1 = nn3 = 0
    for r in order:
        s = state_nodes(r['J'], r['mfl'], r['pe'], r['sigma'], r['pw'])
        pri = 'P1' if (r['J'], round(r['mfl'], 2)) in ((26, 1.03), (30, 1.03), (34, 1.03)) else ('P2' if r['mfl'] < 1.06 else 'P3')
        tag = f"{pri}_J{r['J']}_m{r['mfl']:.3f}"
        D = s['D']; m = s['m']
        nus1 = ','.join(f'{x*m:.6f}' for x in s['n1'])
        lines.append(f"# {tag}: Delta = {D:.4f}, x1 = {s['x1']:.4f}, x_lo = {s['xlo']:.3f}, reference zero x0 = {s['x0']:.4f}, edge from "
                     f"{s['edge_a']:.3f}; tolerance {r['tolE_tr']:.3g} S_E,t")
        lines.append(f"{tag}_n1  n1  {r['J']}  {D:.6f}  4242  {nus1}")
        nn1 += len(s['n1'])
        n3x = list(s['n3']) + list(s['ov'])
        if n3x:
            nus3 = ','.join(f'{x*m:.6f}' for x in n3x)
            lines.append(f"{tag}_n3a  n3  {r['J']}  {D:.6f}  {seed}  {nus3}")
            lines.append(f"{tag}_n3b  n3  {r['J']}  {D:.6f}  {seed + 4096}  {nus3}")
            nn3 += 2*len(n3x)
            seed += 8
        summ.append(f"{tag}: Delta {D:.3f}, n1 nodes {len(s['n1'])}, n3 nodes {len(s['n3'])} (+ {len(s['ov'])} overlap) x 2 seeds")
    lines.append(f"# totals: {nn1} n1 jobs, {nn3} n3 jobs")
    open(EXT, 'w').write('\n'.join(lines) + '\n')
    summ.append(f"totals: {nn1} n1 jobs, {nn3} n3 jobs")
    open(os.path.join(R.CACHE, 'jobs_summary.txt'), 'w').write('\n'.join(summ) + '\n')
    print('\n'.join(summ))


if __name__ == '__main__':
    # the one write outside edger/cache: the NEW file ads4/ext/jobs_batch3.txt (refuse to overwrite anything else)
    target = os.path.realpath(EXT)
    _o = builtins.open

    def _open(f, mode='r', *a, **k):
        if any(c in mode for c in 'wax+') and os.path.realpath(os.fspath(f)) == target:
            return R.er_guard._open(f, mode, *a, **k)
        return _o(f, mode, *a, **k)
    builtins.open = _open
    main()
