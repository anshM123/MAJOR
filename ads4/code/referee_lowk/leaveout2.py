"""signed comparison of the power-law estimate with the true values at left-out accepted points (see leaveout.py)"""
import os, sys, inspect, textwrap
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'e2e')); sys.path.insert(0, os.path.join(HERE, '..', 'lowk'))
import numpy as np
import corner
from corner import TensorMomentsLC
src = textwrap.dedent(inspect.getsource(TensorMomentsLC.moment_err)).replace("self.last_est_share =", "self.last_est = est; self.last_est_share =")
ns = {}; exec(src, corner.__dict__, ns); TensorMomentsLC.moment_err = ns['moment_err']
mom = TensorMomentsLC.load(sys.argv[1]); l = mom.l; g = mom.g
acc0 = mom.conv & g['ok']
print(f"# {os.path.basename(sys.argv[1])}")
for mode in ('deep', 'tiny'):
    for n in (1, 2):
        removed = np.zeros_like(mom.conv)
        if mode == 'deep':
            for j in range(removed.shape[1]):
                rows = np.where(acc0[:, j])[0]
                if len(rows) > n + 2 and (g['ok'][:, j] & ~mom.conv[:, j]).any(): removed[rows[:n], j] = True
        else:
            for i in range(removed.shape[0]):
                cols = np.where(acc0[i, :])[0]
                if len(cols) > n + 2 and (g['ok'][i, :] & ~mom.conv[i, :]).any(): removed[i, cols[:n]] = True
        m2 = TensorMomentsLC.__new__(TensorMomentsLC); m2.__dict__.update(mom.__dict__); m2.conv = mom.conv & ~removed
        for k in (2, 3):
            for nu in (0.0,):
                for (q1, q3) in ((l, -l), (l, l)):
                    (v0, s0), _, ed0, _ = mom.moment_err(k, nu, q1, q3)
                    (v1, s1), _, ed1, _ = m2.moment_err(k, nu, q1, q3)
                    dens, sh = mom._density(k, nu, q1, q3)
                    tr = dens[removed].sum(); es = m2.last_est[removed].sum()   # same scale sh (dens scale identical)
                    tot = dens.sum()
                    print(f"  {mode} n={n} k={k} ({q1},{q3}): true(removed)/total = {complex(tr/tot).real:+.2e}{complex(tr/tot).imag:+.1e}j   est/true = {complex(es/tr).real:+.2f}{complex(es/tr).imag:+.2f}j   (quoted err_drop of full moment {ed0:.1e})")
