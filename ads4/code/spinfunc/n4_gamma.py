"""
n4_gamma (NEARGAP100_SMALLNU.md): common-random-number (CRN) errors of the graviton state/reference ratios at the n4 production points
(cache/n4_raw.pkl): for every kind and matrix element, the ratio R = <E_state>/<E_ref> (thimble sums A+ + A-, each with its own logN)
with a jackknife error over the 8 common-seed replicas, compared with the independent-error estimate that n1_analysis / p8_analysis
use (state and reference spreads added in quadrature).  Delta Gamma ~ m^2 * relative error (c2 units).
Usage: python3 n4_gamma.py -> cache/n4_gamma.txt
"""
import os, sys, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np


def reps_sum(raw, k):
    a, b = raw['th']['A+'][k], raw['th']['A-'][k]
    return a['reps'] + b['reps']*np.exp(b['logN'] - a['logN']), a['logN']


def jk_ratio(xs, xr):
    n = len(xs)
    R = xs.mean()/xr.mean()
    jk = np.array([np.delete(xs, i).mean()/np.delete(xr, i).mean() for i in range(n)])
    lj = np.log(jk/R)
    e = np.sqrt((n - 1)/n*np.sum(np.abs(lj - lj.mean())**2))
    ind = np.hypot(xs.std(ddof=1)/np.sqrt(n)/abs(xs.mean()), xr.std(ddof=1)/np.sqrt(n)/abs(xr.mean()))
    return R, e, ind


def main(pr):
    db = pickle.load(open(os.path.join(HERE, 'cache', 'n4_raw.pkl'), 'rb'))
    keys = sorted(set((k[0], k[2], k[3]) for k in db))
    for (mode, J, nu) in keys:
        Ds = sorted(k[1] for k in db if k[0] == mode and k[2] == J and k[3] == nu)
        if 20000.0 not in Ds or len(Ds) < 2:
            continue
        st = [D for D in Ds if D != 20000.0][0]
        rs, rr = db[(mode, st, J, nu)], db[(mode, 20000.0, J, nu)]
        m = rs['m']
        pr(f"## {mode} J={J} nu={nu:.4f} state Delta={st:.4f} (m={m:.2f}) / ref 2e4")
        for k in rs['kinds']:
            Xs, _ = reps_sum(rs, k); Xr, _ = reps_sum(rr, k)
            na = Xs.shape[1]
            for i in range(na):
                for j in range(Xs.shape[2]):
                    if i != j and na > 1:
                        continue
                    R, e, ind = jk_ratio(Xs[:, i, j], Xr[:, i, j])
                    pr(f"   {k:8s} [{i},{j}]: state/ref rel. error CRN {e:.3f} (indep. {ind:.3f});  m^2 x CRN error = {m*m*e:7.0f} c2")


if __name__ == '__main__':
    lines = []
    def pr(s=""):
        print(s, flush=True); lines.append(s)
    main(pr)
    open(os.path.join(HERE, 'cache', 'n4_gamma.txt'), 'w').write("\n".join(lines) + "\n")
