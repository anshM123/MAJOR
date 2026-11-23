"""dc_tailan (DIRCV.md §3.7): heavy-tail / CV diagnostics from the per-point samples of dc_tail (nu = 13.3, two seeds pooled)."""
import dc_guard, os, glob, pickle
import numpy as np
lines = []
def pr(s=''):
    print(s, flush=True); lines.append(s)
fs = sorted(glob.glob(os.path.join(dc_guard.CACHE, 'tail_J100_D1100_nu13.3_s*.pkl')))
ds = [pickle.load(open(f, 'rb')) for f in fs]
def cat(mode, key, kind=None):
    if kind is None:
        return np.concatenate([d['data'][mode][key] for d in ds])
    return np.concatenate([d['data'][mode][key][kind] for d in ds], axis=-1)

def hill(x, k):
    a = np.sort(np.abs(x))[::-1]
    return 1.0/np.mean(np.log(a[:k]/a[k]))

pr(f"# per-point diagnostics, (1100,100) nu=13.3, files {[os.path.basename(f) for f in fs]}")
rng = np.random.default_rng(1)
for mode in (1, -1, 0):
    S = cat(mode, 'S'); G = cat(mode, 'G', 'EE'); T = G[0, 0] + G[1, 1]
    n = len(S)
    for name, x in (('scalar', S), ('EE trace', T)):
        mu = x.mean()
        # contribution of the largest |x| points to the real part of the mean (relative to the mean)
        o = np.argsort(np.abs(x))[::-1]
        top = [np.sum(x[o[:k]]).real/np.sum(x).real for k in (1, 10, 100)]
        boot = np.array([x[rng.integers(0, n, n)].mean() for _ in range(400)])
        pr(f"  mode {mode:+d} {name:8s}: n={n}  rel err {np.std(x)/np.sqrt(n)/abs(mu):.4f};  Hill alpha(k=50) {hill(x, 50):.2f}, (k=200) {hill(x, 200):.2f};  "
           f"share of Re sum from top 1/10/100 |x|: {top[0]:+.3f}/{top[1]:+.3f}/{top[2]:+.3f};  boot/naive err {np.std(boot)/(np.std(x)/np.sqrt(n)):.2f}")
    # control-variate structure: relative fluctuations g = T/E[T], s = S/E[S]
    g = T/T.mean(); s = S/S.mean()
    beta = np.real(np.vdot(s - 1, g - 1))/np.real(np.vdot(s - 1, s - 1))
    corr = np.corrcoef(np.real(g), np.real(s))[0, 1]
    ratio = T/S
    pr(f"  mode {mode:+d} CV: beta_opt (Re) = {beta:.3f}, corr(Re g, Re s) = {corr:.3f};  point-wise T/S: median {np.median(ratio.real):.4g}, "
       f"IQR/median {np.subtract(*np.percentile(ratio.real, [75, 25]))/abs(np.median(ratio.real)):.3f};  E[T]/E[S] = {(T.mean()/S.mean()).real:.4g}")
    # where does the graviton/scalar ratio deviate: weight by |S| in the top-1% tail vs bulk
    o = np.argsort(np.abs(S))[::-1]; k = n//100
    pr(f"           mean T/S over the top 1% |S| points: {np.mean(ratio[o[:k]]).real:.4g};  over the rest: {np.mean(ratio[o[k:]]).real:.4g}")
open(os.path.join(dc_guard.CACHE, 'dc_tailan.txt'), 'w').write('\n'.join(lines) + '\n')
