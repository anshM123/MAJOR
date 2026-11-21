"""
bg_validate (BELOW_FOLD2.md §3.1): scalar validation of the analytic-radial expansion (bg_scalar) against the exact Mellin skeleton
(n3_exact.X = g C(nu+i) + g~ C(nu-i), records via bf_exact), in the convention of bf_validate:
   m^2 eps(nu) = m^2 [log(R_bg(nu)/R_exact(nu)) - const],  R = state/ref (Delta_ref = 2e4, same J),  const fitted on nu < 1.4.
Also the size of each order: m^2 [log R_K - log R_{K-1}] (c2 units), i.e. the truncation table.
Usage: python3 bg_validate.py [K]  -> cache/bg_validate.txt
"""
import os, sys, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import warnings
warnings.filterwarnings('ignore')
import numpy as np
import bg_scalar as BS, bf_exact as BE, n3_exact as NX

lines = []


def pr(s=''):
    print(s, flush=True); lines.append(s)


STATES = [(100.0, 0), (100.0, 2), (100.0, 4), (100.0, 24), (300.0, 0), (300.0, 2), (300.0, 4), (300.0, 24),
          (1000.0, 0), (1000.0, 2), (1000.0, 4), (1000.0, 24), (336.0, 24), (336.0, 28), (336.0, 32)]
NUSV = np.array([0.1, 0.3, 0.6, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 14.0, 16.0, 18.0])


def check(r, D, J, nus, Kt=None):
    rs, rr = BE.get(D, J), BE.get(2e4, J)
    Rex = np.array([NX.X(nu, rs)[0]/NX.X(nu, rr)[0] for nu in nus])
    m = rs['m']
    R, ms, mr = BS.ratio(r, D, J, Kt=Kt)
    L = np.log(R/Rex); c = np.mean(L[nus < 1.4])
    eps = m*m*(L - c)
    # per-order sizes: m^2 [log R_{<=k} - log R_{<=k-1}] - (same at nu_0), k = 1..K
    cum_s = np.cumsum(ms, 0); cum_r = np.cumsum(mr, 0)
    lR = np.log(cum_s/cum_r)
    dk = m*m*np.diff(lR, axis=0)
    dk = dk - dk[:, :1]
    return eps, dk, m


if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    t0 = time.time()
    r = BS.build(K, NUSV, verbose=False)
    pr(f'# bg_validate: analytic-radial 1/Delta expansion (K = {K}) vs exact Mellin skeleton; m^2 eps in c2 units (relative x m^2)')
    pr('# columns: nu | m^2 eps (Re, Im) | order-k contributions m^2[log R_k - log R_(k-1)] relative to nu_0, k = 1..K')
    for (D, J) in STATES:
        eps, dk, m = check(r, D, J, NUSV)
        pr(f'## Delta = {D}, J = {J} (m = {m:.2f});  max |Re m^2 eps| on nu<=3 / 3-10 / 10-18: '
           f'{np.max(np.abs(eps.real[NUSV <= 3])):.3f} / {np.max(np.abs(eps.real[(NUSV > 3) & (NUSV <= 10)])):.3f} / '
           f'{np.max(np.abs(eps.real[NUSV > 10])):.3f}')
        for i, v in enumerate(NUSV):
            pr(f'   nu = {v:5.2f}: {eps[i].real:+10.3f} {eps[i].imag:+9.3f}   | ' + ' '.join(f'{x.real:+9.2f}' for x in dk[:, i]))
    pr(f'# time {time.time() - t0:.0f}s')
    open(os.path.join(HERE, 'cache', 'bg_validate.txt'), 'w').write('\n'.join(lines) + '\n')
