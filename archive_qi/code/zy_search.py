"""Random-start search for quantum violations of Zhang-Yeung: maximise 2I(C:D)/RHS (>1 would violate)."""
import numpy as np, sys, time
from multiprocessing import Pool
from entropy_opt import *
S="/tmp/claude-0/-home-user-MAJOR/60eea3a7-07d6-5453-af55-c4346d86f9b5/scratchpad/"
def job(args):
    dims, r, seed, trials = args
    D = int(np.prod(dims)); lhs, rhs = zy_parts()
    num = EntropyFunctional(dims, comb((-1, lhs))); den = EntropyFunctional(dims, rhs)
    rng = np.random.default_rng(seed); best = (0, None)
    for t in range(trials):
        x0 = rng.normal(size=2*D*r)
        if t % 3 == 1:   # near-classical start: dominant diagonal-ish structure
            x0 *= 0.05; k = rng.integers(0, D, size=r); x0[(k*r + np.arange(r)) % (D*r)] += 1
        f, x = optimise(num, den, D, r, x0, "ratio", maxiter=1500)
        if f < best[0]: best = (f, x)
    return dims, r, -best[0], best[1]
if __name__ == "__main__":
    cfgs = []
    for dims in [(2,2,2,2)]:
        for r in [1, 2, 3, 4, 8, 16]: cfgs += [(dims, r, 100+s, 12) for s in range(4)]
    for dims in [(3,2,2,2), (2,3,2,2), (2,2,3,2), (2,2,3,3), (3,3,2,2)]:
        for r in [1, 2, 4]: cfgs += [(dims, r, 200+s, 6) for s in range(2)]
    t0 = time.time(); res = {}
    with Pool(4) as p:
        for dims, r, val, x in p.imap_unordered(job, cfgs):
            key = (dims, r)
            if key not in res or val > res[key][0]: res[key] = (val, x)
            print(f"{time.time()-t0:7.0f}s dims={dims} rank={r}: max ratio {val:.8f}", flush=True)
    print("SUMMARY")
    for key in sorted(res): print(key, f"{res[key][0]:.8f}")
    best = max(res.items(), key=lambda kv: kv[1][0]); np.save(S+"zy_best.npy", best[1][1]); print("best", best[0], best[1][0])
