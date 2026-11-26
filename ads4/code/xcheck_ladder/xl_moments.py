"""
xcheck_ladder/xl_moments.py: power moments S_j = sum_i w_i g(nu_i) nu_i^{2j} of a discrete realisation of the k = 2 weight,
in Arb, split into nu-ranges, written to out/pmom_<tag>.json (arb strings).

  python3 xl_moments.py cont 1000 NUP JMAX PREC   # Gauss-Legendre (NUP nodes/panel) on the continuum (repaired) weight x theta
  python3 xl_moments.py grid 1000 0   JMAX PREC   # the published midpoint grid (read as zero beyond Lam + 10)
  python3 xl_moments.py grid 2000 0   JMAX PREC
"""
import sys, os, json, time, math
sys.dont_write_bytecode = True
import numpy as np
import xl_common as X
from flint import ctx


def main():
    kind, Dg, nup, jmax, prec = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5])
    Lam = math.sqrt(X.lam_int(Dg))
    t0 = time.time()
    if kind == 'cont':
        W = X.load_cont(Dg)
        nu, w = X.cont_nodes(Lam, nup)
        chi = X.chi_re(W, nu)
        g = X.g_of_chi(chi, nu, Lam)*X.theta(nu, Lam)
        tag = f'cont_{Dg}_n{nup}_p{prec}'
    else:
        nu, w, chi = X.grid_measure(Dg)
        g = X.g_of_chi(chi, nu, Lam)
        tag = f'grid_{Dg}_p{prec}'
    wg = w*g
    print(f'{tag}: {len(nu)} nodes, nu_max = {nu.max():.3f}, sum|w g| = {np.abs(wg).sum():.6e}, eval {time.time()-t0:.1f}s', flush=True)
    edges = [0, Lam, Lam + 5, Lam + 10, Lam + 20, Lam + 30, Lam + 40, 1e30]
    ranges = list(zip(edges[:-1], edges[1:]))
    ranges[0] = (-1, Lam)
    t0 = time.time()
    S = X.power_moments(nu, wg, jmax, prec, ranges)
    print(f'  power moments j <= {jmax}: {time.time()-t0:.1f}s', flush=True)
    absw = [float(np.abs(wg[(nu > lo) & (nu <= hi)]).sum()) for lo, hi in ranges]
    out = dict(tag=tag, Dg=Dg, prec=prec, jmax=jmax, n_nodes=int(len(nu)), ranges=[[float(a), float(b)] for a, b in ranges],
               abs_wg=absw, S=[[s.str(int(prec*0.30103) + 5, radius=True) for s in Sr] for Sr in S])
    json.dump(out, open(os.path.join(X.OUT, f'pmom_{tag}.json'), 'w'))
    print('  written', f'out/pmom_{tag}.json', flush=True)


if __name__ == '__main__':
    main()
