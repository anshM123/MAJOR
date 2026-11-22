"""
rd_corner (RD_CLOSE.md section 3.3): EE weight on the crystal-thimble flow lines beyond the n3/n4 action cut (s > 20), per terminal
divisor, from the rd_flow caches.  EE weight of a thimble point relative to the near region = exp(-s) |Q|/Q_eff.
Usage: python3 rd_corner.py -> cache/rd_corner.txt
"""
import os, sys, glob, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
np.seterr(all='ignore')
LAB = ['m45=0', 'm35=0', 'X1->0', 'X2->0', 'Y1->0', 'Y2->0', 's->+-1', 'y^2=0 corner', 'sin psi->0']
lines = ['# rd_corner: EE weight log(e^{-s}|Q|/Q_eff) at the action checkpoints s = 10, 20, 40, 80, ... along the crystal-thimble lines',
         '# columns: file saddle | terminal divisor: #lines, max EE log-weight at s >= 20, fraction of those lines whose EE weight rises between the last two checkpoints']
for f in sorted(glob.glob(os.path.join(HERE, 'cache', 'rd_flow_D*_nu*.pkl'))):
    R = pickle.load(open(f, 'rb'))
    for k, o in R.items():
        if 'sck' not in o:
            continue
        lQe = np.log(o['qnear'])
        ew = o['lqck'] - o['sck'] - lQe                     # EE log-weight at checkpoints (relative to the near region)
        L = o['last']; okl = np.all(np.isfinite(L), axis=0)
        am = np.full(L.shape[1], -1); am[okl] = np.argmin(np.log(L[:, okl] + 1e-300), axis=0)
        parts = []
        for j in range(9):
            sel = am == j
            if sel.sum() == 0:
                continue
            e2 = ew[sel][:, 1:]
            mx = np.nanmax(e2) if np.isfinite(e2).any() else np.nan
            # rising between last two finite checkpoints
            rise = 0; tot = 0
            for row in ew[sel]:
                fin = np.where(np.isfinite(row))[0]
                if fin.size >= 2:
                    tot += 1; rise += row[fin[-1]] > row[fin[-2]]
            parts.append(f'{LAB[j]}: {sel.sum()}, {mx:+.1f}, {rise}/{tot}')
        lines.append(f'{os.path.basename(f)[8:-4]:>14s} {k:2s} | ' + ' ; '.join(parts))
open(os.path.join(HERE, 'cache', 'rd_corner.txt'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(lines))
