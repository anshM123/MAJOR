"""
Partial completion of the cached <JJJJ> mode grids (fourpt/cache/modes4) at the points dropped by the radial-series convergence
mask (tiny t, r -> 1): the dropped points that converge with N2 radial levels (default 30000, vs N <= 10^4 in the production runs) are
evaluated with the unchanged fourpt machinery (SpinBlock4.eval_coeffs: words, jets of order 4, generated code), the seed blocks being
built one at a time (memory).  The remaining dropped points are estimated by corner_est.  Output: cache/modes4c_*.npz (same format as
TensorMoments4.save, with conv extended and the N2 points filled).  The lowk corner method (ads4/code/lowk) does not apply here:
it needs c1 + c2 > 0 for every seed, while <JJJJ> has seeds with c1 + c2 = 0.
Usage: python3 complete_n2.py <modes4 file name> [N2]
"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, os.path.join(HERE, '..', '..', 'e2e'))
from blocks import CBlock
from moments import frame_data, frame_positions, series_converged
from spinblock2 import SpinBlock4, invariants4, monomial_values
from moments4 import pol_grid4, fourier4
from tmoments4 import TensorMoments4, DATADIR, basis_coeffs
CACHE = os.path.join(HERE, 'cache')

class LazyBlocks(dict):
    """dict (a,b) -> (block, vpow) whose blocks are built on access and released afterwards (one at a time)."""
    def __init__(self, specs, Delta, J, N):
        super().__init__()
        self.specs, self.Delta, self.J, self.N = specs, Delta, J, N
        for k in specs: dict.__setitem__(self, k, None)
    def items(self):
        for ab, (c1, c2, vpow) in self.specs.items():
            blk = CBlock(self.Delta, self.J, c1, c2, self.N)
            yield ab, (blk, vpow)
            del blk

class SpinBlock4N2(SpinBlock4):
    def __init__(self, l, Delta, J, N):
        self.l, self.Delta, self.J, self.N, self.prec = l, float(Delta), int(J), int(N), 'd'
        from structures2 import words
        from spinblock2 import phi_funcs
        from embed2 import shifts
        self.words = words(l); self.phi = phi_funcs(l); self.ab = {}
        specs = {}
        Dx = l + 1
        for wL in self.words:
            for wR in self.words:
                s1, s2 = shifts(wL); s3, s4 = shifts(wR)
                a = -((Dx + s1) - (Dx + s2))/2; b = ((Dx + s3) - (Dx + s4))/2
                self.ab[(wL, wR)] = (a, b)
                if (a, b) not in specs:
                    c1, c2 = -a, b
                    specs[(a, b)] = (c1, c2, 0.0) if c1 + c2 >= 0 else (-c1, -c2, c1 + c2)
        self.blocks = LazyBlocks(specs, Delta, J, N)

def complete(fname, N2=30000, chunk=300):
    mom = TensorMoments4.load(os.path.join(DATADIR, fname))
    g = mom.g; ok = g['ok']
    fdall = frame_data(g['U'][ok], g['S'][ok])
    conv_new = series_converged(fdall['rs'], fdall['chi'], mom.J, N2, cut=25.0)
    conv_old = mom.conv[ok]
    sel = np.where(conv_new & ~conv_old)[0]
    print(f"# {fname}: {ok.sum()} grid points, {conv_old.sum()} converged at N <= 10^4, {len(sel)} more at N2 = {N2}, "
          f"{(~conv_new).sum()} still dropped", flush=True)
    fd = {k: v[sel] for k, v in fdall.items()}
    ws, vecs = basis_coeffs(1, mom.Delta, mom.J)
    combos = {(a, b): (vecs[a], vecs[b]) for a in range(2) for b in range(2)}
    t0 = time.time()
    sb = SpinBlock4N2(1, mom.Delta, mom.J, N2)
    fs, monos = sb.phi
    new = {ab: {} for ab in combos}; ls_new = np.zeros(len(sel))
    for c0 in range(0, len(sel), chunk):
        sl = slice(c0, min(len(sel), c0 + chunk))
        x1, x2, x3, x4 = frame_positions(fd['ru'][sl], fd['xi'][sl])
        z1, z3, z2, z4 = pol_grid4(fd['xi'][sl], 1)
        s, z = invariants4({1: x1, 2: x2, 3: x3, 4: x4}, {1: z1, 2: z2, 3: z3, 4: z4})
        lsc, C = sb.eval_coeffs(s, fd['rs'][sl], fd['chi'][sl], fd['A0'][sl], fd['E0'][sl], combos)
        ls_new[sl] = lsc
        mv = monomial_values(monos, z)
        for ab in combos:
            M = np.einsum('mp,mp...->p...', C[ab], mv)
            for q, v in fourier4(M, 1).items():
                new[ab].setdefault(q, []).append(v)
        print(f"  points {sl.start}-{sl.stop}: {time.time() - t0:.0f} s", flush=True)
    # merge
    idx_ok = np.flatnonzero(ok.ravel())
    flat = idx_ok[sel]
    for ab in combos:
        for q, lst in new[ab].items():
            arr = mom.modes[ab][q].copy().ravel(); arr[flat] = np.concatenate(lst); mom.modes[ab][q] = arr.reshape(g['S'].shape)
    lsf = mom.logscale.copy().ravel(); lsf[flat] = ls_new; mom.logscale = lsf.reshape(g['S'].shape)
    cf = mom.conv.copy().ravel(); cf[flat] = True; mom.conv = cf.reshape(g['S'].shape)
    os.makedirs(CACHE, exist_ok=True)
    out = os.path.join(CACHE, fname.replace('modes4_', 'modes4c_'))
    mom.save(out)
    print(f"  saved {out} ({time.time() - t0:.0f} s)", flush=True)
    return out

if __name__ == "__main__":
    complete(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 30000)
