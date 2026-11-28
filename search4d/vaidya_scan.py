"""Physical control scan: concentric annuli in Vaidya-AdS4 (M=1).

For a pool of boundary edges (r_p, t_p) on a spacelike boundary slice
t = t0 + kappa (r - r0), compute every cap branch and every tube branch
(tubes found from BOTH ends as a branch-completeness check), take the minimal
branch, then:
  * K-slacks (SSA/WM instances of Lemma B) on the physical data,
  * SSA/WM/MMI on all unions of intervals for random sub-pools,
  * all 372 facets for random assignments of annuli to 5 parties,
reporting min margin per facet orbit.
Usage: python3 -I vaidya_scan.py dlt t0 kappa npool seed out.json
"""
import sys, os, json, time, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from vaidya4 import Vaidya
from hrt_pieces import Pieces
from abstract_axis import F, ORB, entropy_vector, entropy, random_labels

INF = 1e300


def main():
    dlt = float(sys.argv[1]); t0 = float(sys.argv[2]); kap = float(sys.argv[3])
    npool = int(sys.argv[4]); seed = int(sys.argv[5]); outp = sys.argv[6]
    rng = np.random.default_rng(seed)
    V = Vaidya(M=1.0, dlt=dlt)
    P = Pieces(V, cap_grid=(220, 200), tube_grid=(181, 181))
    T = time.time()
    radii = np.sort(np.exp(rng.uniform(np.log(0.25), np.log(2.4), npool)))
    r0 = radii[0]
    times = t0 + kap * (radii - r0)
    P.cap_scan(vrange=(min(times) - 7.0, max(times) + 0.5), zrange=(0.02, 6.0))
    c = np.full(npool, INF); ncap = []
    for i in range(npool):
        sols = P.caps(radii[i], times[i])
        ncap.append(len(sols))
        if sols:
            c[i] = min(r[2] for x, r in sols)
    L = np.full((npool, npool), INF); nbr = np.zeros((npool, npool), int); mism = 0.0
    for i in range(npool):
        for j in range(npool):
            if i == j:
                continue
            sols = P.tubes(radii[i], times[i], radii[j], times[j])
            nbr[i, j] = len(sols)
            a = min((r[2] for x, r in sols), default=INF)
            if i < j:
                L[i, j] = a
            else:
                b = L[j, i]
                if (a < INF) != (b < INF):
                    # a branch seen from one end only: record (could be a missed branch)
                    mism = max(mism, 1.0 if min(a, b) < c[i] + c[j] else 0.5)
                elif a < INF:
                    mism = max(mism, abs(a - b))
                L[j, i] = L[i, j] = min(a, b)
    tpieces = time.time() - T
    # K-slacks with L' = min(L, c+c)
    Lp = np.minimum(L, c[:, None] + c[None, :])
    ks = []
    for i, j, k in itertools.combinations(range(npool), 3):
        ks.append(Lp[i, k] + c[j] - Lp[i, j] - c[k]); ks.append(Lp[i, k] + c[j] - c[i] - Lp[j, k])
    for i, j, k, l in itertools.combinations(range(npool), 4):
        ks.append(Lp[i, k] + Lp[j, l] - Lp[i, j] - Lp[k, l]); ks.append(Lp[i, k] + Lp[j, l] - Lp[i, l] - Lp[j, k])
    kmin = float(min(ks)) if ks else None
    # facets on random sub-pools / labellings
    worst = np.full(len(ORB), np.inf); worst_info = [None] * len(ORB)
    nvec = 0
    for trial in range(400):
        m = int(rng.integers(5, min(9, npool) + 1))
        sub = np.sort(rng.choice(npool, m, replace=False))
        Ls = L[np.ix_(sub, sub)]; cs = c[sub]
        for _ in range(10):
            lab = random_labels(m, rng)
            v = entropy_vector(lab, Ls, cs)
            if np.any(v >= INF / 2):
                continue
            mg = F @ v; nvec += 1
            for o in range(len(ORB)):
                if o == 2:
                    continue
                k = int(np.argmin(mg[ORB[o]]))
                if mg[ORB[o]][k] < worst[o]:
                    worst[o] = mg[ORB[o]][k]
                    worst_info[o] = dict(sub=sub.tolist(), lab=lab, facet=int(ORB[o][k]))
    res = dict(dlt=dlt, t0=t0, kappa=kap, radii=radii.tolist(), times=times.tolist(), c=c.tolist(),
               ncap_branches=ncap, L=np.where(L < INF, L, None).tolist(), tube_branches=nbr.tolist(),
               two_end_mismatch=mism, K_min_slack=kmin, n_vectors=nvec,
               min_margin_per_orbit=[None if not np.isfinite(w) else float(w) for w in worst],
               worst_info=worst_info, seconds_pieces=tpieces, seconds_total=time.time() - T)
    json.dump(res, open(outp, 'w'), indent=1)
    print(json.dumps({k: res[k] for k in ('dlt', 't0', 'kappa', 'ncap_branches', 'two_end_mismatch', 'K_min_slack',
                                          'n_vectors', 'min_margin_per_orbit', 'seconds_total')}))


if __name__ == '__main__':
    main()
