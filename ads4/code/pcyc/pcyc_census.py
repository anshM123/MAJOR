"""
pcyc_census (PCYC_CLOSE.md section 4): extension of RD_CLOSE's pilot census of critical points of the scalar flow function Phi_w
(rd/rd_census.run, imported read-only) to more (Delta, nu) and more starts.  rd_saddle's saddle cache is redirected to pcyc/cache (the
RD cache is not touched).  Output: cache/census_D<Delta>_nu<nu>_s<seed>.txt / .pkl, plus a one-line summary:
  counts by raw radius r: Regge region (r < 0.1), FAR region (0.1 <= r < 10: the RD far region |raw x| ~ 1 and beyond the B zeros up
  to 10), raw infinity (r >= 10); the max height dH (sheet-independent part) above the crystal in each class.
Usage: python3 pcyc_census.py Delta nu nstart seed
"""
import os, sys, pickle, time
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
RDD = os.path.join(HERE, '..', 'rd')
sys.path.insert(0, RDD)
import numpy as np
import rd_saddle as SD
SD.HERE = HERE                       # saddle cache -> pcyc/cache (read-only use of the RD module)
import rd_census as RC


def main(Delta, nu, nstart, seed):
    t = time.time()
    rng_patch = np.random.default_rng
    np.random.default_rng = lambda s=None: rng_patch(seed)        # rd_census.run uses default_rng(seed=3): give it ours
    try:
        out, L = RC.run(Delta, nu, nstart)
    finally:
        np.random.default_rng = rng_patch
    r = out['r']; dH = out['Hs'] - out['Hc'][0]
    cls = {'Regge r<0.1': r < 0.1, 'mid 0.1<=r<0.3': (r >= 0.1) & (r < 0.3), 'FAR 0.3<=r<1.5': (r >= 0.3) & (r < 1.5),
           'beyond 1.5<=r<10': (r >= 1.5) & (r < 10), 'raw inf r>=10': r >= 10}
    lB = out['lB']
    S = [f'# pcyc_census Delta={Delta:g} nu={nu:g} starts={nstart} seed={seed}: {out["nconv"]} converged, {r.size} distinct; {time.time() - t:.0f}s']
    for k, m in cls.items():
        S.append(f'#   {k:17s}: {int(m.sum()):4d} points; max dH {dH[m].max() if m.any() else float("nan"):+8.2f}; '
                 f'(Delta/2)log|B| in [{lB[m].min() if m.any() else float("nan"):+.2f}, {lB[m].max() if m.any() else float("nan"):+.2f}] '
                 f'(crystal {out["lBc"][0]:+.2f}); r range [{r[m].min() if m.any() else float("nan"):.3g}, {r[m].max() if m.any() else float("nan"):.3g}]')
    L = S + L
    tag = f'D{Delta:g}_nu{nu:g}_s{seed}'
    open(os.path.join(HERE, 'cache', f'census_{tag}.txt'), 'w').write('\n'.join(L) + '\n')
    pickle.dump(out, open(os.path.join(HERE, 'cache', f'census_{tag}.pkl'), 'wb'))
    print('\n'.join(S), flush=True)


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
